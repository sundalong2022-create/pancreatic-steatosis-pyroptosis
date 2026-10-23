#!/usr/bin/env Rscript
# 12-week PSPS robustness analyses.
# Leave-one-core-gene-out, random six-gene nulls, label permutations,
# bootstrap CIs, GSVA ssGSEA scores, non-core scores, and forest plots.

suppressPackageStartupMessages({
  library(limma); library(GSVA); library(ggplot2); library(dplyr)
  library(tidyr); library(tibble)
})
options(stringsAsFactors = FALSE)
set.seed(292815)

root <- normalizePath('.', mustWork=TRUE)
base <- file.path(root, 'results/12wk_primary')
outdir <- file.path(base, 'robustness')
figdir <- file.path(outdir, 'figures')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr
meta <- obj$metadata
stopifnot(ncol(expr) == nrow(meta), all(colnames(expr) == meta$matrix_column), ncol(expr) == 21)
logexpr <- log2(expr + 1)

gs_obj <- readRDS(file.path(root, 'scripts/results/gene_sets_mouse/gene_sets.rds'))
core <- unique(gs_obj$core_axis)
pathways <- gs_obj$gene_sets
noncore_inflammation <- c('Tnf','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl2','Ccl3','Ccl4','Ccl5',
                          'Nfkb1','Rela','Tlr4','Tgfb1','Ptgs2','S100a8','S100a9','Nfkbia','Tnfaip3')
sets <- list(
  PSPS_core = core,
  pyroptosis_noncore = setdiff(pathways$pyroptosis, core),
  inflammasome_noncore = setdiff(pathways$inflammasome, core),
  inflammation_noncore = noncore_inflammation
)
sets_present <- lapply(sets, function(x) intersect(unique(x), rownames(logexpr)))
coverage <- bind_rows(lapply(names(sets), function(nm) {
  tibble(set=nm, requested=length(unique(sets[[nm]])), core_overlap=sum(unique(sets[[nm]]) %in% core),
         n_present=length(sets_present[[nm]]), coverage=length(sets_present[[nm]])/length(unique(sets[[nm]])),
         genes_present=paste(sets_present[[nm]], collapse=';'))
}))
write.csv(coverage, file.path(outdir, '12wk_robustness_gene_set_coverage.csv'), row.names=FALSE)

score_z <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[genes, , drop=FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm=TRUE)
}
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(unname(cor(x[ok], y[ok], method='spearman')))
}

planned <- c(
  primary__HO_Safflower_vs_Safflower='HO_Safflower_12wk-Safflower_12wk',
  primary__HO_Sunflower_vs_Sunflower='HO_Sunflower_12wk-Sunflower_12wk'
)
design <- model.matrix(~0 + group, data=meta)
colnames(design) <- sub('^group', '', colnames(design))
group_indices <- list(
  Safflower=list(treated=which(meta$group == 'HO_Safflower_12wk'), control=which(meta$group == 'Safflower_12wk')),
  Sunflower=list(treated=which(meta$group == 'HO_Sunflower_12wk'), control=which(meta$group == 'Sunflower_12wk'))
)

contrast_table <- function(score_df, method_prefix) {
  score_cols <- setdiff(colnames(score_df), c('sample','group'))
  score_mat <- t(as.matrix(score_df[, score_cols, drop=FALSE]))
  rownames(score_mat) <- score_cols
  fit0 <- lmFit(score_mat, design)
  out <- bind_rows(lapply(names(planned), function(nm) {
    ff <- eBayes(contrasts.fit(fit0, makeContrasts(contrasts=planned[[nm]], levels=design)), robust=TRUE)
    tt <- topTable(ff, number=Inf, sort.by='none') %>% rownames_to_column('score')
    se <- ff$coefficients[,1] / ff$t[,1]
    tt$SE <- unname(se[tt$score])
    tt$CI95_lower <- tt$logFC - qt(0.975, df=ff$df.total[1]) * tt$SE
    tt$CI95_upper <- tt$logFC + qt(0.975, df=ff$df.total[1]) * tt$SE
    tt$contrast <- nm
    tt$method <- method_prefix
    tt
  }))
  out$padj_BH <- ave(out$P.Value, out$contrast, FUN=function(x) p.adjust(x, 'BH'))
  out
}

# Full PSPS and non-core scores.
full_psps <- score_z(logexpr, core)
noncore_scores <- data.frame(
  sample=colnames(expr),
  pyroptosis_noncore_score=score_z(logexpr, sets$pyroptosis_noncore),
  inflammasome_noncore_score=score_z(logexpr, sets$inflammasome_noncore),
  inflammation_noncore_score=score_z(logexpr, sets$inflammation_noncore),
  group=meta$group, check.names=FALSE
)
write.csv(noncore_scores, file.path(outdir, '12wk_noncore_scores.csv'), row.names=FALSE)
noncore_contrasts <- contrast_table(noncore_scores %>% select(sample, pyroptosis_noncore_score, inflammasome_noncore_score, inflammation_noncore_score, group), 'mean_z_noncore')
write.csv(noncore_contrasts, file.path(outdir, '12wk_noncore_score_contrasts.csv'), row.names=FALSE)

# Leave one of the six core genes out, recalculate PSPS from the remaining five,
# and repeat the planned contrasts and independent inflammation correlation.
loo_scores <- data.frame(sample=colnames(expr), group=meta$group, PSPS_core=full_psps, check.names=FALSE)
for (g in core) loo_scores[[paste0('PSPS_without_', g)]] <- score_z(logexpr, setdiff(core, g))
write.csv(loo_scores, file.path(outdir, '12wk_psps_leave_one_gene_out_scores.csv'), row.names=FALSE)
loo_contrasts <- contrast_table(loo_scores, 'leave_one_gene_out')
write.csv(loo_contrasts, file.path(outdir, '12wk_psps_leave_one_gene_out_contrasts.csv'), row.names=FALSE)
loo_cor <- bind_rows(lapply(setdiff(colnames(loo_scores), c('sample','group')), function(sc) {
  ok <- is.finite(loo_scores[[sc]]) & is.finite(noncore_scores$inflammation_noncore_score)
  ct <- if (sum(ok) >= 4) suppressWarnings(cor.test(loo_scores[[sc]][ok], noncore_scores$inflammation_noncore_score[ok], method='spearman', exact=FALSE)) else NULL
  tibble(score=sc, n=sum(ok), spearman_rho=if(is.null(ct)) NA_real_ else unname(ct$estimate), p_value=if(is.null(ct)) NA_real_ else ct$p.value)
}))
loo_cor$padj_BH <- p.adjust(loo_cor$p_value, 'BH')
write.csv(loo_cor, file.path(outdir, '12wk_psps_leave_one_gene_out_correlations.csv'), row.names=FALSE)

# Standard GSVA ssGSEA implementation (GSVA 2.x).
ssgsea_sets <- lapply(sets[c('PSPS_core','pyroptosis_noncore','inflammasome_noncore')], function(x) intersect(unique(x), rownames(logexpr)))
ssgsea_param <- GSVA::ssgseaParam(logexpr, ssgsea_sets, alpha=0.25, normalize=TRUE)
ssgsea_mat <- as.matrix(GSVA::gsva(ssgsea_param, verbose=FALSE))
colnames(ssgsea_mat) <- colnames(expr)
ssgsea_df <- as.data.frame(t(ssgsea_mat)) %>% rownames_to_column('sample') %>%
  rename(PSPS_core_ssGSEA=PSPS_core, pyroptosis_noncore_ssGSEA=pyroptosis_noncore, inflammasome_noncore_ssGSEA=inflammasome_noncore) %>%
  left_join(meta %>% select(matrix_column, group), by=c('sample'='matrix_column'))
write.csv(ssgsea_df, file.path(outdir, '12wk_ssgsea_scores.csv'), row.names=FALSE)
ssgsea_contrasts <- contrast_table(ssgsea_df %>% select(sample, PSPS_core_ssGSEA, pyroptosis_noncore_ssGSEA, inflammasome_noncore_ssGSEA, group), 'ssGSEA')
write.csv(ssgsea_contrasts, file.path(outdir, '12wk_ssgsea_contrasts.csv'), row.names=FALSE)

# Score correlations against the independent non-core inflammation module.
cor_score_df <- bind_rows(
  tibble(score='PSPS_core', value=full_psps),
  tibble(score='pyroptosis_noncore_score', value=noncore_scores$pyroptosis_noncore_score),
  tibble(score='inflammasome_noncore_score', value=noncore_scores$inflammasome_noncore_score),
  tibble(score='PSPS_core_ssGSEA', value=ssgsea_df$PSPS_core_ssGSEA),
  tibble(score='pyroptosis_noncore_ssGSEA', value=ssgsea_df$pyroptosis_noncore_ssGSEA),
  tibble(score='inflammasome_noncore_ssGSEA', value=ssgsea_df$inflammasome_noncore_ssGSEA)
)
cor_tbl <- cor_score_df %>% group_by(score) %>% summarise(n=n(), spearman_rho=safe_spearman(value, noncore_scores$inflammation_noncore_score), .groups='drop')
cor_tbl$p_value <- vapply(cor_tbl$score, function(sc) {
  v <- cor_score_df$value[cor_score_df$score == sc]
  suppressWarnings(cor.test(v, noncore_scores$inflammation_noncore_score, method='spearman', exact=FALSE)$p.value)
}, numeric(1))
cor_tbl$padj_BH <- p.adjust(cor_tbl$p_value, 'BH')
write.csv(cor_tbl, file.path(outdir, '12wk_robustness_score_correlations.csv'), row.names=FALSE)
psps_noncore_cor <- bind_rows(lapply(c('pyroptosis_noncore_score','inflammasome_noncore_score','inflammation_noncore_score'), function(nm) {
  ct <- suppressWarnings(cor.test(full_psps, noncore_scores[[nm]], method='spearman', exact=FALSE))
  tibble(score=nm, n=length(full_psps), spearman_rho=unname(ct$estimate), p_value=ct$p.value)
}))
psps_noncore_cor$padj_BH <- p.adjust(psps_noncore_cor$p_value, 'BH')
write.csv(psps_noncore_cor, file.path(outdir, '12wk_psps_noncore_correlations.csv'), row.names=FALSE)

# Random six-gene score null distribution.
n_perm <- as.integer(Sys.getenv('PSPS_PERM_N', '5000'))
universe <- rownames(logexpr)[apply(logexpr, 1, sd) > 0]
universe <- setdiff(universe, core)
if (length(universe) < 6) stop('Random-gene universe has fewer than six genes.')
random_rows <- vector('list', n_perm)
for (i in seq_len(n_perm)) {
  genes_i <- sample(universe, length(core), replace=FALSE)
  sc <- score_z(logexpr, genes_i)
  a <- group_indices$Safflower; b <- group_indices$Sunflower
  random_rows[[i]] <- tibble(iteration=i, gene1=genes_i[1], gene2=genes_i[2], gene3=genes_i[3], gene4=genes_i[4], gene5=genes_i[5], gene6=genes_i[6],
    effect_HO_Safflower_vs_Safflower=mean(sc[a$treated])-mean(sc[a$control]),
    effect_HO_Sunflower_vs_Sunflower=mean(sc[b$treated])-mean(sc[b$control]),
    rho_inflammation_noncore=safe_spearman(sc, noncore_scores$inflammation_noncore_score))
}
random_df <- bind_rows(random_rows)
write.csv(random_df, file.path(outdir, '12wk_random_six_gene_null_distribution.csv'), row.names=FALSE)
observed_effects <- c(primary__HO_Safflower_vs_Safflower=mean(full_psps[group_indices$Safflower$treated])-mean(full_psps[group_indices$Safflower$control]),
                      primary__HO_Sunflower_vs_Sunflower=mean(full_psps[group_indices$Sunflower$treated])-mean(full_psps[group_indices$Sunflower$control]))
random_summary <- bind_rows(lapply(names(observed_effects), function(nm) {
  col <- if (grepl('Safflower', nm)) 'effect_HO_Safflower_vs_Safflower' else 'effect_HO_Sunflower_vs_Sunflower'
  null <- random_df[[col]]; obs <- observed_effects[[nm]]
  tibble(test='PSPS effect', contrast=nm, observed=obs, null_mean=mean(null), null_sd=sd(null), null_q025=quantile(null,.025), null_q975=quantile(null,.975),
         empirical_p_two_sided=(sum(abs(null) >= abs(obs))+1)/(length(null)+1), empirical_p_upper=(sum(null >= obs)+1)/(length(null)+1), n_perm=length(null))
}))
obs_rho <- safe_spearman(full_psps, noncore_scores$inflammation_noncore_score)
rho_null <- random_df$rho_inflammation_noncore
random_summary <- bind_rows(random_summary, tibble(test='PSPS vs inflammation_noncore rho', contrast='all_12wk_samples', observed=obs_rho, null_mean=mean(rho_null), null_sd=sd(rho_null), null_q025=quantile(rho_null,.025), null_q975=quantile(rho_null,.975), empirical_p_two_sided=(sum(abs(rho_null)>=abs(obs_rho))+1)/(length(rho_null)+1), empirical_p_upper=(sum(rho_null>=obs_rho)+1)/(length(rho_null)+1), n_perm=length(rho_null)))
write.csv(random_summary, file.path(outdir, '12wk_random_six_gene_null_summary.csv'), row.names=FALSE)

# Label permutation null for the fixed PSPS score in each primary contrast.
n_label_perm <- as.integer(Sys.getenv('PSPS_LABEL_PERM_N', '5000'))
label_summary <- bind_rows(lapply(names(observed_effects), function(nm) {
  oil <- if (grepl('Safflower', nm)) 'Safflower' else 'Sunflower'
  idx_t <- group_indices[[oil]]$treated; idx_c <- group_indices[[oil]]$control; idx <- c(idx_t, idx_c)
  nt <- length(idx_t); nc <- length(idx_c); obs <- observed_effects[[nm]]
  vals <- vapply(seq_len(n_label_perm), function(i) { z <- sample(c(rep(TRUE, nt), rep(FALSE, nc))); mean(full_psps[idx[z]]) - mean(full_psps[idx[!z]]) }, numeric(1))
  tibble(contrast=nm, observed=obs, null_mean=mean(vals), null_sd=sd(vals), null_q025=quantile(vals,.025), null_q975=quantile(vals,.975), empirical_p_two_sided=(sum(abs(vals)>=abs(obs))+1)/(length(vals)+1), empirical_p_upper=(sum(vals>=obs)+1)/(length(vals)+1), n_perm=length(vals))
}))
write.csv(label_summary, file.path(outdir, '12wk_psps_label_permutation_summary.csv'), row.names=FALSE)

# Bootstrap CIs: within-group resampling for effects, sample-level resampling for correlations.
n_boot <- as.integer(Sys.getenv('PSPS_BOOTSTRAP_B', '5000'))
boot_effects <- bind_rows(lapply(names(observed_effects), function(nm) {
  oil <- if (grepl('Safflower', nm)) 'Safflower' else 'Sunflower'
  idx_t <- group_indices[[oil]]$treated; idx_c <- group_indices[[oil]]$control
  vals <- vapply(seq_len(n_boot), function(b) mean(full_psps[sample(idx_t, length(idx_t), replace=TRUE)]) - mean(full_psps[sample(idx_c, length(idx_c), replace=TRUE)]), numeric(1))
  tibble(metric='PSPS_core effect', contrast=nm, estimate=observed_effects[[nm]], ci95_lower=quantile(vals,.025), ci95_upper=quantile(vals,.975), n_boot=n_boot, n_valid=sum(is.finite(vals)))
}))
rho_metrics <- list(
  PSPS_vs_inflammation_noncore=list(x=full_psps, y=noncore_scores$inflammation_noncore_score),
  PSPS_vs_pyroptosis_noncore=list(x=full_psps, y=noncore_scores$pyroptosis_noncore_score),
  PSPS_vs_inflammasome_noncore=list(x=full_psps, y=noncore_scores$inflammasome_noncore_score)
)
boot_rhos <- bind_rows(lapply(names(rho_metrics), function(nm) {
  x <- rho_metrics[[nm]]$x; y <- rho_metrics[[nm]]$y
  vals <- vapply(seq_len(n_boot), function(b) { ii <- sample(seq_along(x), length(x), replace=TRUE); safe_spearman(x[ii], y[ii]) }, numeric(1))
  tibble(metric=nm, contrast='all_12wk_samples', estimate=safe_spearman(x,y), ci95_lower=quantile(vals,.025,na.rm=TRUE), ci95_upper=quantile(vals,.975,na.rm=TRUE), n_boot=n_boot, n_valid=sum(is.finite(vals)))
}))
boot_summary <- bind_rows(boot_effects, boot_rhos)
write.csv(boot_summary, file.path(outdir, '12wk_bootstrap_confidence_intervals.csv'), row.names=FALSE)

# Forest plot combines the primary PSPS, LOO, non-core and ssGSEA effects.
forest_scores <- bind_rows(
  loo_contrasts %>% filter(grepl('PSPS_without_', score)) %>% transmute(score, logFC, CI95_lower, CI95_upper, contrast, source='limma'),
  noncore_contrasts %>% transmute(score, logFC, CI95_lower, CI95_upper, contrast, source='limma'),
  ssgsea_contrasts %>% filter(score %in% c('PSPS_core_ssGSEA','pyroptosis_noncore_ssGSEA','inflammasome_noncore_ssGSEA')) %>% transmute(score, logFC, CI95_lower, CI95_upper, contrast, source='limma'),
  tibble(score='PSPS_core_bootstrap', logFC=observed_effects[['primary__HO_Safflower_vs_Safflower']], CI95_lower=boot_effects$ci95_lower[boot_effects$contrast=='primary__HO_Safflower_vs_Safflower'], CI95_upper=boot_effects$ci95_upper[boot_effects$contrast=='primary__HO_Safflower_vs_Safflower'], contrast='primary__HO_Safflower_vs_Safflower', source='bootstrap'),
  tibble(score='PSPS_core_bootstrap', logFC=observed_effects[['primary__HO_Sunflower_vs_Sunflower']], CI95_lower=boot_effects$ci95_lower[boot_effects$contrast=='primary__HO_Sunflower_vs_Sunflower'], CI95_upper=boot_effects$ci95_upper[boot_effects$contrast=='primary__HO_Sunflower_vs_Sunflower'], contrast='primary__HO_Sunflower_vs_Sunflower', source='bootstrap')
)
forest_scores <- forest_scores %>% mutate(score_label=case_when(
  score=='PSPS_core_bootstrap' ~ 'PSPS core (bootstrap)',
  score=='pyroptosis_noncore_score' ~ 'Pyroptosis non-core mean-z',
  score=='inflammasome_noncore_score' ~ 'Inflammasome non-core mean-z',
  score=='inflammation_noncore_score' ~ 'Inflammation non-core mean-z',
  score=='PSPS_core_ssGSEA' ~ 'PSPS core ssGSEA',
  score=='pyroptosis_noncore_ssGSEA' ~ 'Pyroptosis non-core ssGSEA',
  score=='inflammasome_noncore_ssGSEA' ~ 'Inflammasome non-core ssGSEA',
  TRUE ~ gsub('PSPS_without_', 'PSPS without ', score)
))
forest_scores$score_label <- factor(forest_scores$score_label, levels=rev(unique(forest_scores$score_label)))
forest_scores$y <- as.numeric(forest_scores$score_label)
write.csv(forest_scores, file.path(outdir, '12wk_robustness_forest_plot_data.csv'), row.names=FALSE)
forest_plot <- ggplot(forest_scores, aes(x=logFC, y=y, colour=source)) +
  geom_vline(xintercept=0, linetype=2, colour='grey45') +
  geom_segment(aes(x=CI95_lower, xend=CI95_upper, y=y, yend=y), linewidth=.55) +
  geom_point(size=2.1) + scale_y_continuous(breaks=seq_along(levels(forest_scores$score_label)), labels=levels(forest_scores$score_label)) +
  facet_wrap(~contrast, scales='free_y', ncol=2) + scale_colour_manual(values=c(bootstrap='#D95F02', limma='#1B9E77')) +
  theme_bw(base_size=10) + theme(legend.position='bottom', panel.grid.minor=element_blank(), axis.title.y=element_blank()) +
  labs(title='PSPS robustness forest plot', x='Estimated score difference (high-oleic minus ordinary oil)', colour='CI source')
ggsave(file.path(figdir, '12wk_psps_robustness_forest_plot.pdf'), forest_plot, width=14, height=11)
ggsave(file.path(figdir, '12wk_psps_robustness_forest_plot.png'), forest_plot, width=14, height=11, dpi=220)

manifest <- c(
  'Analysis: 12-week PSPS robustness audit',
  paste0('Samples: ', ncol(expr), '; groups: ', paste(names(table(meta$group)), as.integer(table(meta$group)), collapse='; ')),
  paste0('Core PSPS genes: ', paste(core, collapse=', ')),
  paste0('Random six-gene permutations: ', n_perm), paste0('Label permutations: ', n_label_perm), paste0('Bootstrap replicates: ', n_boot),
  'Leave-one-gene-out: each core gene removed and PSPS recalculated from the remaining five genes.',
  'Random null: six genes sampled uniformly from expressed non-core genes with non-zero cross-sample variance.',
  'Bootstrap: within-group resampling for PSPS effects; sample-level resampling for Spearman correlations.',
  'Alternative scoring: GSVA 2.x ssGSEA implementation with PSPS_core and non-core pyroptosis/inflammasome sets.',
  'Forest plot: LOO, non-core mean-z, ssGSEA, and bootstrap PSPS effects with 95% CIs.',
  'Interpretation: robustness analyses support transcriptomic score stability and do not establish functional pyroptosis.'
)
writeLines(manifest, file.path(outdir, 'README.md'))
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo_robustness.txt'))
message('Completed PSPS robustness analysis: ', outdir)
