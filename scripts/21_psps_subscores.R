#!/usr/bin/env Rscript
# PSPS-Pyro and PSPS-Inf subscore analysis for the 12-week mouse cohort.
# PSPS six-gene core is not changed. Mean-z scores are primary; ssGSEA is a
# sensitivity analysis. LOO, within-group, group-adjusted lipid and conditional
# immune analyses are kept in this independent output module. Immune-module FDR
# families are new and do not modify the existing PSPS immune results.

suppressPackageStartupMessages({
  library(GSVA); library(limma); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork)
})
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
outdir <- file.path(root, 'results/12wk_primary/subscores')
figdir <- file.path(root, 'results/12wk_primary/lipid_external/figures')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

core <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
pyro <- c('Nlrp3','Casp1','Gsdmd')
inf <- c('Il1b','Il18')
planned <- c(
  primary__HO_Safflower_vs_Safflower = 'HO_Safflower_12wk-Safflower_12wk',
  primary__HO_Sunflower_vs_Sunflower = 'HO_Sunflower_12wk-Sunflower_12wk'
)

obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr) == nrow(meta), all(colnames(expr) == meta$matrix_column), ncol(expr) == 21)
logexpr <- log2(expr + 1)

score_z <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[genes,,drop=FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm=TRUE)
}
safe_cor <- function(x, y, min_n=3) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(rho=NA_real_, p=NA_real_, n=n))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method='spearman', exact=FALSE))
  c(rho=unname(ct$estimate), p=ct$p.value, n=n)
}
partial_spearman <- function(x, y, group) {
  ok <- is.finite(x) & is.finite(y) & !is.na(group)
  n <- sum(ok); g <- droplevels(factor(group[ok]))
  if (n < 5 || nlevels(g) < 2) return(c(rho=NA_real_, p=NA_real_, n=n, df=NA_real_))
  X <- model.matrix(~g)
  rx <- residuals(lm(rank(x[ok], ties.method='average') ~ X - 1))
  ry <- residuals(lm(rank(y[ok], ties.method='average') ~ X - 1))
  r <- suppressWarnings(cor(rx, ry))
  df <- n - qr(X)$rank - 1
  if (!is.finite(r) || abs(r) >= 1 || df <= 0) p <- NA_real_ else p <- 2*pt(abs(r)*sqrt(df/(1-r^2)), df=df, lower.tail=FALSE)
  c(rho=unname(r), p=p, n=n, df=df)
}
partial_spearman_covariate <- function(x, y, covariate) {
  ok <- is.finite(x) & is.finite(y) & is.finite(covariate)
  n <- sum(ok)
  if (n < 5) return(c(rho=NA_real_, p=NA_real_, n=n, df=NA_real_))
  xr <- rank(x[ok], ties.method='average')
  yr <- rank(y[ok], ties.method='average')
  cr <- rank(covariate[ok], ties.method='average')
  rx <- residuals(lm(xr ~ cr))
  ry <- residuals(lm(yr ~ cr))
  r <- suppressWarnings(cor(rx, ry))
  df <- n - 3
  if (!is.finite(r) || abs(r) >= 1 || df <= 0) p <- NA_real_ else p <- 2*pt(abs(r)*sqrt(df/(1-r^2)), df=df, lower.tail=FALSE)
  c(rho=unname(r), p=p, n=n, df=df)
}

# Confirm the frozen PSPS core is exactly the existing six-gene score.
old <- read.csv(file.path(root,'results/12wk_primary/scores/12wk_sample_scores.csv'), check.names=FALSE)
old <- old[match(colnames(expr), old$sample),]
psps <- score_z(logexpr, core)
if (anyNA(old$sample) || max(abs(psps - old$PSPS_core), na.rm=TRUE) > 1e-10) stop('Frozen PSPS mismatch.')

# Primary mean-z and ssGSEA subscores.
mean_scores <- data.frame(sample=colnames(expr),
  PSPS_Pyro_mean_z=score_z(logexpr, pyro),
  PSPS_Inf_mean_z=score_z(logexpr, inf),
  PSPS_core=psps, group=meta$group, check.names=FALSE)
ss_sets <- list(PSPS_Pyro=pyro, PSPS_Inf=inf)
ss_present <- lapply(ss_sets, function(g) intersect(g, rownames(logexpr)))
ss_param <- GSVA::ssgseaParam(logexpr, ss_present, alpha=0.25, normalize=TRUE)
ss_mat <- as.matrix(GSVA::gsva(ss_param, verbose=FALSE))
ss_scores <- data.frame(sample=colnames(expr),
  PSPS_Pyro_ssGSEA=as.numeric(ss_mat['PSPS_Pyro', colnames(expr)]),
  PSPS_Inf_ssGSEA=as.numeric(ss_mat['PSPS_Inf', colnames(expr)]),
  group=meta$group, check.names=FALSE)
scores <- left_join(mean_scores, ss_scores, by=c('sample','group'))
write.csv(scores, file.path(outdir,'12wk_psps_subscores_sample_scores.csv'), row.names=FALSE)
write.csv(data.frame(axis=c('PSPS-Pyro','PSPS-Inf'), mean_z_genes=c(paste(pyro,collapse=';'),paste(inf,collapse=';')), ssGSEA_present_n=c(length(ss_present$PSPS_Pyro),length(ss_present$PSPS_Inf))), file.path(outdir,'12wk_psps_subscore_definitions.csv'), row.names=FALSE)

# Shared limma contrast helper.
design <- model.matrix(~0+group, data=meta); colnames(design) <- sub('^group','',colnames(design))
contrast_table <- function(score_df, method) {
  sc <- setdiff(colnames(score_df), c('sample','group'))
  fit0 <- lmFit(t(as.matrix(score_df[,sc,drop=FALSE])), design)
  out <- bind_rows(lapply(names(planned), function(nm) {
    ff <- eBayes(contrasts.fit(fit0, makeContrasts(contrasts=planned[[nm]], levels=design)), robust=TRUE)
    tt <- topTable(ff, number=Inf, sort.by='none') %>% rownames_to_column('score')
    se <- ff$coefficients[,1]/ff$t[,1]; se[!is.finite(se)] <- NA_real_
    tt %>% transmute(method=method, contrast=nm, score, estimate=logFC, SE=unname(se[score]),
      CI95_lower=logFC-qt(.975,df=ff$df.total[1])*SE,
      CI95_upper=logFC+qt(.975,df=ff$df.total[1])*SE, p_value=P.Value)
  }))
  out$q_within_contrast <- ave(out$p_value, out$contrast, FUN=function(x) p.adjust(x,'BH'))
  out
}
contrast_all <- contrast_table(scores %>% select(sample, PSPS_Pyro_mean_z, PSPS_Inf_mean_z, PSPS_Pyro_ssGSEA, PSPS_Inf_ssGSEA, group), 'subscore')
write.csv(contrast_all, file.path(outdir,'12wk_psps_subscore_contrasts.csv'), row.names=FALSE)

# Required LOO: the primary mean-z score is recalculated after removing each gene.
loo <- data.frame(sample=colnames(expr), group=meta$group, check.names=FALSE)
for (g in pyro) loo[[paste0('PSPS_Pyro_without_',g)]] <- score_z(logexpr, setdiff(pyro,g))
for (g in inf) loo[[paste0('PSPS_Inf_without_',g)]] <- score_z(logexpr, setdiff(inf,g))
write.csv(loo, file.path(outdir,'12wk_psps_subscore_LOO_scores.csv'), row.names=FALSE)
loo_contrasts <- contrast_table(loo, 'LOO_mean_z')
write.csv(loo_contrasts, file.path(outdir,'12wk_psps_subscore_LOO_contrasts.csv'), row.names=FALSE)

# External lipid scores generated by the pre-specified external module.
ext <- read.csv(file.path(root,'results/12wk_primary/lipid_external/12wk_external_lipid_sample_scores.csv'), check.names=FALSE)
ext <- ext[match(colnames(expr), ext$sample),]
if (anyNA(ext$sample)) stop('External lipid scores do not match the 21-sample cohort.')
lipid_scores <- c('Hallmark_FA_Metabolism','Hallmark_Adipogenesis','KEGG_PPAR_Signaling','Reactome_FA_Metabolism')
subscore_names <- c('PSPS_Pyro_mean_z','PSPS_Inf_mean_z','PSPS_Pyro_ssGSEA','PSPS_Inf_ssGSEA')

# Pooled correlations: one external-lipid family of 16 tests
# (2 biological axes x 2 scoring methods x 4 prespecified lipid sets).
pooled <- bind_rows(lapply(subscore_names, function(sc) bind_rows(lapply(lipid_scores, function(lp) {
  z <- safe_cor(scores[[sc]], ext[[lp]])
  tibble(method=ifelse(grepl('mean_z',sc),'mean-z','ssGSEA'), subscore=sub('_(mean_z|ssGSEA)$','',sc), score=sc, lipid_score=lp, rho=z[['rho']], p_value=z[['p']], n=z[['n']])
}))))
pooled$q_external_lipid_family <- p.adjust(pooled$p_value,'BH')
write.csv(pooled, file.path(outdir,'12wk_psps_subscore_external_lipid_pooled.csv'), row.names=FALSE)

# Within-group Spearman correlations: all groups and all 80 prespecified pairs
# (5 groups x 2 biological axes x 2 scoring methods x 4 lipid sets). These
# estimates are retained for confounding diagnostics, not treated as primary
# evidence because each group has only 3-5 samples.
within <- bind_rows(lapply(levels(factor(meta$group)), function(gr) bind_rows(lapply(subscore_names, function(sc) bind_rows(lapply(lipid_scores, function(lp) {
  ix <- meta$group == gr; z <- safe_cor(scores[[sc]][ix], ext[[lp]][ix])
  tibble(group=gr, method=ifelse(grepl('mean_z',sc),'mean-z','ssGSEA'), subscore=sub('_(mean_z|ssGSEA)$','',sc), score=sc, lipid_score=lp, rho=z[['rho']], p_value=z[['p']], n=z[['n']])
}))))))
within$q_within_group_family <- p.adjust(within$p_value,'BH')
write.csv(within, file.path(outdir,'12wk_psps_subscore_external_lipid_within_group.csv'), row.names=FALSE)

# Group-adjusted partial Spearman: rank residuals after the five-level group
# factor; one 16-test external-lipid family, kept separate from the pooled
# and within-group families.
partial <- bind_rows(lapply(subscore_names, function(sc) bind_rows(lapply(lipid_scores, function(lp) {
  z <- partial_spearman(scores[[sc]], ext[[lp]], meta$group)
  tibble(method=ifelse(grepl('mean_z',sc),'mean-z','ssGSEA'), subscore=sub('_(mean_z|ssGSEA)$','',sc), score=sc, lipid_score=lp, partial_rho=z[['rho']], p_value=z[['p']], n=z[['n']], df=z[['df']])
}))))
partial$q_group_adjusted_family <- p.adjust(partial$p_value,'BH')
write.csv(partial, file.path(outdir,'12wk_psps_subscore_external_lipid_partial_group.csv'), row.names=FALSE)

# LOO external-lipid correlations and contrasts demonstrate that no single gene is required.
loo_cor <- bind_rows(lapply(setdiff(names(loo),c('sample','group')), function(sc) bind_rows(lapply(lipid_scores, function(lp) {
  z <- safe_cor(loo[[sc]], ext[[lp]])
  tibble(score=sc, lipid_score=lp, rho=z[['rho']], p_value=z[['p']], n=z[['n']])
}))))
loo_cor$q_BH <- p.adjust(loo_cor$p_value,'BH')
write.csv(loo_cor, file.path(outdir,'12wk_psps_subscore_LOO_external_lipid_correlations.csv'), row.names=FALSE)

# Immune modules are inherited from the existing core-excluded module scores.
imm <- read.csv(file.path(root,'results/12wk_primary/immune/12wk_immune_sample_scores.csv'), check.names=FALSE)
imm <- imm[match(colnames(expr), imm$sample),]
mods <- setdiff(colnames(imm), c('sample','PSPS_core','group'))
immune_cor <- bind_rows(lapply(subscore_names, function(sc) bind_rows(lapply(mods, function(md) {
  z <- safe_cor(scores[[sc]], imm[[md]])
  tibble(method=ifelse(grepl('mean_z',sc),'mean-z','ssGSEA'), subscore=sub('_(mean_z|ssGSEA)$','',sc), score=sc, module=md, rho=z[['rho']], p_value=z[['p']], n=z[['n']])
}))))
immune_cor$FDR_family <- paste0(immune_cor$method,'__',immune_cor$subscore,'__immune_modules')
immune_cor$q_within_family <- ave(immune_cor$p_value, immune_cor$FDR_family, FUN=function(x) p.adjust(x,'BH'))
write.csv(immune_cor, file.path(outdir,'12wk_psps_subscore_immune_correlations.csv'), row.names=FALSE)

# LOO versus immune modules, for audit completeness.
loo_immune <- bind_rows(lapply(setdiff(names(loo),c('sample','group')), function(sc) bind_rows(lapply(mods, function(md) {
  z <- safe_cor(loo[[sc]], imm[[md]])
  tibble(score=sc, module=md, rho=z[['rho']], p_value=z[['p']], n=z[['n']])
}))))
loo_immune$FDR_BH <- p.adjust(loo_immune$p_value,'BH')
write.csv(loo_immune, file.path(outdir,'12wk_psps_subscore_LOO_immune_correlations.csv'), row.names=FALSE)

# Conditional immune correlations: control the opposing biological axis in
# rank space, then report both directions for every module. The two directions
# are kept in one FDR family per scoring method (18 tests per method).
conditional <- bind_rows(lapply(list(mean_z=c('PSPS_Pyro_mean_z','PSPS_Inf_mean_z'),
                                      ssGSEA=c('PSPS_Pyro_ssGSEA','PSPS_Inf_ssGSEA')), function(pair) {
  method <- if (grepl('mean_z', pair[[1]])) 'mean-z' else 'ssGSEA'
  pyro_sc <- scores[[pair[[1]]]]; inf_sc <- scores[[pair[[2]]]]
  pyro_rows <- bind_rows(lapply(mods, function(md) {
    z0 <- safe_cor(pyro_sc, imm[[md]])
    z <- partial_spearman_covariate(pyro_sc, imm[[md]], inf_sc)
    tibble(method=method, target_subscore='PSPS_Pyro', conditioned_on='PSPS_Inf', module=md,
           zero_order_rho=z0[['rho']], partial_rho=z[['rho']], p_value=z[['p']], n=z[['n']], df=z[['df']])
  }))
  inf_rows <- bind_rows(lapply(mods, function(md) {
    z0 <- safe_cor(inf_sc, imm[[md]])
    z <- partial_spearman_covariate(inf_sc, imm[[md]], pyro_sc)
    tibble(method=method, target_subscore='PSPS_Inf', conditioned_on='PSPS_Pyro', module=md,
           zero_order_rho=z0[['rho']], partial_rho=z[['rho']], p_value=z[['p']], n=z[['n']], df=z[['df']])
  }))
  bind_rows(pyro_rows, inf_rows)
}))
conditional$FDR_family <- paste0(conditional$method,'__conditional_immune_modules')
conditional$q_within_family <- ave(conditional$p_value, conditional$FDR_family, FUN=function(x) p.adjust(x,'BH'))
conditional$partial_minus_zero <- conditional$partial_rho - conditional$zero_order_rho
write.csv(conditional, file.path(outdir,'12wk_psps_subscore_conditional_immune_correlations.csv'), row.names=FALSE)

# Full zero-order Spearman matrices are retained so the conditional analysis
# can be read alongside the Pyro-Inf and module correlations without selective
# reporting. One long-format matrix is written per scoring method.
zero_order_matrix <- bind_rows(lapply(list(mean_z=c('PSPS_Pyro_mean_z','PSPS_Inf_mean_z'),
                                            ssGSEA=c('PSPS_Pyro_ssGSEA','PSPS_Inf_ssGSEA')), function(pair) {
  method <- if (grepl('mean_z', pair[[1]])) 'mean-z' else 'ssGSEA'
  dat <- data.frame(PSPS_Pyro=scores[[pair[[1]]]], PSPS_Inf=scores[[pair[[2]]]],
                    imm[,mods,drop=FALSE], check.names=FALSE)
  mat <- cor(dat, method='spearman', use='pairwise.complete.obs')
  as.data.frame(as.table(mat), stringsAsFactors=FALSE) %>%
    transmute(method=method, variable_i=as.character(Var1), variable_j=as.character(Var2),
              spearman_rho=Freq)
}))
write.csv(zero_order_matrix, file.path(outdir,'12wk_psps_subscore_immune_zero_order_correlation_matrix.csv'), row.names=FALSE)
zero_order_wide <- bind_rows(lapply(split(zero_order_matrix, zero_order_matrix$method), function(d) {
  wide <- tidyr::pivot_wider(d, names_from=variable_j, values_from=spearman_rho)
  wide %>% select(variable_i, everything()) %>% mutate(method=unique(d$method), .before=1)
}))
write.csv(zero_order_wide, file.path(outdir,'12wk_psps_subscore_immune_zero_order_correlation_matrix_wide.csv'), row.names=FALSE)

# Combined panel: existing PSPS-vs-lipid evidence plus subscore-vs-lipid sensitivity.
label_map <- c(Hallmark_FA_Metabolism='Hallmark FA metabolism', Hallmark_Adipogenesis='Hallmark adipogenesis', KEGG_PPAR_Signaling='KEGG PPAR signalling', Reactome_FA_Metabolism='Reactome FA metabolism')
old_cor <- read.csv(file.path(root,'results/12wk_primary/lipid_external/12wk_psps_external_lipid_correlations.csv'), check.names=FALSE)
old_con <- read.csv(file.path(root,'results/12wk_primary/lipid_external/12wk_external_lipid_contrasts.csv'), check.names=FALSE)
cl <- c(primary__HO_Safflower_vs_Safflower='HO safflower vs safflower', primary__HO_Sunflower_vs_Sunflower='HO sunflower vs sunflower')
old_df <- bind_rows(
  old_cor %>% transmute(score=score, label=unname(label_map[score]), analysis='Frozen PSPS correlation', estimate=rho, lower=NA_real_, upper=NA_real_, q=q_global_external),
  old_con %>% transmute(score=score, label=unname(label_map[score]), analysis=unname(cl[contrast]), estimate=estimate, lower=CI95_lower, upper=CI95_upper, q=q_global_external)
) %>% mutate(label=factor(label, levels=rev(unname(label_map))), analysis=factor(analysis, levels=c('Frozen PSPS correlation','HO safflower vs safflower','HO sunflower vs sunflower')))
p_old <- ggplot(old_df,aes(y=label,x=estimate)) + geom_vline(xintercept=0,linetype=2,colour='grey55') + geom_segment(data=filter(old_df,is.finite(lower)),aes(x=lower,xend=upper,yend=label),linewidth=.7) + geom_point(size=2.2,colour='#0072B2') + geom_text(aes(x=ifelse(analysis=='Frozen PSPS correlation',estimate+.04,estimate),hjust=ifelse(analysis=='Frozen PSPS correlation',0,.5),label=paste0('q=',formatC(q,format='e',digits=2))),vjust=-.9,size=2.4) + facet_grid(.~analysis,scales='free_x') + theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(),axis.title.y=element_blank(),legend.position='none',strip.background=element_rect(fill='grey92')) + labs(x='rho or ssGSEA score difference',title='Frozen PSPS and external lipid programmes')
sub_df <- pooled %>% mutate(label=ifelse(method=='mean-z',paste0(subscore,' mean-z'),paste0(subscore,' ssGSEA')), lipid_label=unname(label_map[lipid_score]), label=factor(label, levels=rev(c('PSPS_Inf ssGSEA','PSPS_Pyro ssGSEA','PSPS_Inf mean-z','PSPS_Pyro mean-z'))), lipid_label=factor(lipid_label,levels=unname(label_map)))
p_sub <- ggplot(sub_df,aes(x=rho,y=label,colour=method)) + geom_vline(xintercept=0,linetype=2,colour='grey55') + geom_point(size=2.1) + geom_text(aes(label=paste0('q=',formatC(q_external_lipid_family,format='e',digits=2))),hjust=-.08,size=2.1,show.legend=FALSE) + facet_wrap(~lipid_label,nrow=1) + scale_colour_manual(values=c('mean-z'='#D55E00','ssGSEA'='#009E73')) + theme_bw(base_size=8.5) + theme(panel.grid.minor=element_blank(),axis.title.y=element_blank(),legend.position='bottom',strip.background=element_rect(fill='grey92')) + labs(x='Spearman rho',colour='Subscore method',title='PSPS-Pyro and PSPS-Inf versus external lipid programmes')
combined <- p_old / p_sub + plot_layout(heights=c(1.35,1))
ggsave(file.path(figdir,'12wk_external_lipid_subscore_combined_panel.pdf'),combined,width=14,height=10,useDingbats=FALSE)
ggsave(file.path(figdir,'12wk_external_lipid_subscore_combined_panel.png'),combined,width=14,height=10,dpi=220)

manifest <- c(
  'Analysis: frozen PSPS-Pyro and PSPS-Inf subscores in GSE292815 12-week mouse pancreas.',
  'PSPS-Pyro = mean gene-wise z-score of Nlrp3, Casp1, Gsdmd; PSPS-Inf = mean gene-wise z-score of Il1b, Il18.',
  'Primary subscores use mean-z; GSVA 2.x ssGSEA (alpha=0.25, normalize=TRUE) is sensitivity analysis.',
  'LOO: each Pyro gene and each Inf gene was removed and the mean-z subscore recalculated; GSVA 2.x does not support the one-gene Inf LOO gene set, so ssGSEA LOO was not treated as a stable gene-set score.',
  'External lipid analyses: pooled Spearman (one 16-test family), within-group Spearman (one 80-test diagnostic family; n=3-5 per group), and group-adjusted partial Spearman (one independent 16-test family).',
  'Primary-contrast q values are BH-adjusted within each contrast across the four subscore/method combinations; LOO contrast q values are adjusted separately within contrast.',
  'Immune correlations: modules are inherited from the existing six-core-excluded immune scores; each subscore/method has its own BH family across immune modules.',
  'Conditional immune correlations: rank-based partial Spearman correlations control the opposing subscore; both directions are reported for all modules, with one 18-test FDR family per scoring method. Full zero-order matrices are provided.',
  'Conditional results: 12wk_psps_subscore_conditional_immune_correlations.csv; zero-order matrices: 12wk_psps_subscore_immune_zero_order_correlation_matrix.csv and 12wk_psps_subscore_immune_zero_order_correlation_matrix_wide.csv.',
  'Existing PSPS and external-lipid outputs were not changed. Combined panel is supplementary to the existing external-lipid panel.'
)
writeLines(manifest,file.path(outdir,'README.md'))
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo_12wk_subscores.txt'))
message('Completed PSPS-Pyro and PSPS-Inf subscore analysis: ',outdir)
