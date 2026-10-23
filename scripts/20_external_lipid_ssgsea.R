#!/usr/bin/env Rscript
# External, pre-specified lipid-program validation for the frozen six-gene PSPS.
#
# This module deliberately does not change PSPS or the existing robustness audit.
# It scores four externally defined lipid programs with the same GSVA 2.x ssGSEA
# implementation, then reports all PSPS correlations and both planned dietary
# contrasts. All 12 external tests are retained and BH-adjusted as one declared
# external-validation family; within-family q-values are also provided for the
# correlation and each contrast family.

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
options(stringsAsFactors = FALSE)

root <- normalizePath('.', mustWork = TRUE)
outdir <- file.path(root, 'results/12wk_primary/lipid_external')
figdir <- file.path(outdir, 'figures')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

required_core <- c('Txnip', 'Nlrp3', 'Casp1', 'Gsdmd', 'Il1b', 'Il18')
planned <- c(
  primary__HO_Safflower_vs_Safflower = 'HO_Safflower_12wk-Safflower_12wk',
  primary__HO_Sunflower_vs_Sunflower = 'HO_Sunflower_12wk-Sunflower_12wk'
)

# ----------------------------- input and invariants --------------------------
obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr
meta <- obj$metadata
stopifnot(ncol(expr) == nrow(meta),
          all(colnames(expr) == meta$matrix_column),
          ncol(expr) == 21)
logexpr <- log2(expr + 1)

gs_obj <- readRDS(file.path(root, 'scripts/results/gene_sets_mouse/gene_sets.rds'))
core <- unique(gs_obj$core_axis)
if (!setequal(core, required_core) || length(core) != 6) {
  stop('PSPS core-axis invariant failed: expected exactly Txnip, Nlrp3, Casp1, Gsdmd, Il1b, Il18; found ',
       paste(core, collapse = ', '))
}
core_present <- intersect(core, rownames(logexpr))
if (length(core_present) != 6) stop('The frozen six-gene PSPS is not fully represented in the 12-week matrix.')

score_z <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}
psps <- score_z(logexpr, required_core)

# Reconcile the score with the already published 12-week PSPS output. This is a
# guard against accidental changes to the frozen six-gene implementation.
existing_file <- file.path(root, 'results/12wk_primary/scores/12wk_sample_scores.csv')
existing_psps <- read.csv(existing_file, check.names = FALSE)
existing_psps <- existing_psps[match(colnames(expr), existing_psps$sample), ]
if (anyNA(existing_psps$sample) || max(abs(psps - existing_psps$PSPS_core), na.rm = TRUE) > 1e-10) {
  stop('Recomputed PSPS does not exactly match results/12wk_primary/scores/12wk_sample_scores.csv.')
}

# ----------------------------- pre-specified gene sets ------------------------
# Hallmark, KEGG and Reactome are retrieved from the installed msigdbr release.
# Human MSigDB definitions are ortholog-mapped to mouse symbols by msigdbr.
fetch_set <- function(collection, subcollection, target_name) {
  args <- list(db_species = 'HS', species = 'Mus musculus', collection = collection)
  if (!is.null(subcollection) && !is.na(subcollection)) args$subcollection <- subcollection
  tbl <- do.call(msigdbr::msigdbr, args)
  hit <- tbl %>% filter(.data$gs_name == .env$target_name)
  if (!nrow(hit)) stop('Requested gene set was not found: ', target_name)
  genes <- unique(hit$gene_symbol[!is.na(hit$gene_symbol) & nzchar(hit$gene_symbol)])
  list(genes = genes, table = hit)
}

spec <- tibble(
  score = c('Hallmark_FA_Metabolism', 'Hallmark_Adipogenesis',
            'KEGG_PPAR_Signaling', 'Reactome_FA_Metabolism'),
  role = c('primary', 'sensitivity', 'primary', 'sensitivity'),
  collection = c('H', 'H', 'C2', 'C2'),
  subcollection = c(NA_character_, NA_character_, 'CP:KEGG_LEGACY', 'CP:REACTOME'),
  gs_name = c('HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_ADIPOGENESIS',
              'KEGG_PPAR_SIGNALING_PATHWAY', 'REACTOME_FATTY_ACID_METABOLISM'),
  concept = c('fatty-acid metabolism', 'adipocyte differentiation/infiltration',
              'PPAR lipid-homeostasis signalling', 'fatty-acid metabolism (Reactome)')
)

set_objs <- lapply(seq_len(nrow(spec)), function(i)
  fetch_set(spec$collection[i], spec$subcollection[i], spec$gs_name[i]))
names(set_objs) <- spec$score
sets_requested <- lapply(set_objs, function(x) x$genes)
sets_present <- lapply(sets_requested, function(g) intersect(g, rownames(logexpr)))

gene_provenance <- bind_rows(lapply(seq_len(nrow(spec)), function(i) {
  x <- set_objs[[spec$score[i]]]$table
  genes <- sets_requested[[spec$score[i]]]
  present <- sets_present[[spec$score[i]]]
  tibble(
    score = spec$score[i], role = spec$role[i], concept = spec$concept[i],
    gs_name = spec$gs_name[i], collection = spec$collection[i],
    subcollection = spec$subcollection[i],
    db_version = paste(unique(x$db_version), collapse = ';'),
    source_species = paste(unique(x$gs_source_species), collapse = ';'),
    exact_source = paste(unique(x$gs_exact_source), collapse = ';'),
    gs_url = paste(unique(x$gs_url), collapse = ';'),
    requested_n = length(genes), unique_n = length(unique(genes)),
    present_n = length(present), coverage = length(present) / length(genes),
    core_overlap = paste(intersect(genes, required_core), collapse = ';'),
    core_overlap_n = length(intersect(genes, required_core)),
    genes_present = paste(present, collapse = ';'),
    genes_missing = paste(setdiff(genes, rownames(logexpr)), collapse = ';')
  )
}))
write.csv(gene_provenance, file.path(outdir, 'external_lipid_gene_set_provenance.csv'), row.names = FALSE)

gene_sets_long <- bind_rows(lapply(names(sets_requested), function(nm)
  tibble(score = nm, requested = TRUE, gene_symbol = sets_requested[[nm]],
         present_in_matrix = sets_requested[[nm]] %in% rownames(logexpr))))
write.csv(gene_sets_long, file.path(outdir, 'external_lipid_gene_sets_long.csv'), row.names = FALSE)

# Pairwise overlap audit is reported on the requested gene sets, independently
# of expression coverage.
overlap <- bind_rows(lapply(seq_len(nrow(spec)), function(i) {
  bind_rows(lapply(seq_len(nrow(spec)), function(j) {
    a <- sets_requested[[spec$score[i]]]; b <- sets_requested[[spec$score[j]]]
    tibble(score1 = spec$score[i], score2 = spec$score[j],
           n_score1 = length(a), n_score2 = length(b),
           overlap_n = length(intersect(a, b)),
           union_n = length(union(a, b)),
           jaccard = length(intersect(a, b)) / length(union(a, b)))
  }))
}))
write.csv(overlap, file.path(outdir, 'external_lipid_gene_set_overlap.csv'), row.names = FALSE)

if (any(gene_provenance$present_n < 10)) stop('At least one external lipid gene set has fewer than 10 expressed genes.')

# ----------------------------- ssGSEA scores ----------------------------------
ssgsea_param <- GSVA::ssgseaParam(logexpr, sets_present, alpha = 0.25, normalize = TRUE)
ssgsea_mat <- as.matrix(GSVA::gsva(ssgsea_param, verbose = FALSE))
if (!all(names(sets_present) %in% rownames(ssgsea_mat))) stop('GSVA output names do not match requested lipid gene sets.')
ssgsea_mat <- ssgsea_mat[names(sets_present), colnames(expr), drop = FALSE]

sample_scores <- tibble(
  sample = colnames(expr), PSPS_core = psps,
  as.data.frame(t(ssgsea_mat), check.names = FALSE),
  group = meta$group, dataset = 'GSE292815_KC_12wk'
)
score_names <- spec$score
write.csv(sample_scores, file.path(outdir, '12wk_external_lipid_sample_scores.csv'), row.names = FALSE)

# ----------------------------- planned comparisons ----------------------------
design <- model.matrix(~0 + group, data = meta)
colnames(design) <- sub('^group', '', colnames(design))
score_mat <- t(as.matrix(sample_scores[, score_names, drop = FALSE]))
rownames(score_mat) <- score_names
fit0 <- limma::lmFit(score_mat, design)

contrast_results <- bind_rows(lapply(names(planned), function(nm) {
  ff <- limma::eBayes(limma::contrasts.fit(fit0, limma::makeContrasts(contrasts = planned[[nm]], levels = design)), robust = TRUE)
  tt <- limma::topTable(ff, number = Inf, sort.by = 'none') %>% rownames_to_column('score')
  se <- ff$coefficients[, 1] / ff$t[, 1]
  se[!is.finite(se)] <- NA_real_
  tt %>% transmute(
    test_type = 'dietary_contrast', contrast = nm, score,
    estimate = logFC, logFC, SE = unname(se[score]),
    CI95_lower = logFC - qt(0.975, df = ff$df.total[1]) * SE,
    CI95_upper = logFC + qt(0.975, df = ff$df.total[1]) * SE,
    p_value = P.Value,
    n_treated = sum(meta$group == ifelse(grepl('HO_Safflower', nm), 'HO_Safflower_12wk', 'HO_Sunflower_12wk')),
    n_control = sum(meta$group == ifelse(grepl('HO_Safflower', nm), 'Safflower_12wk', 'Sunflower_12wk'))
  )
}))
contrast_results$q_within_family <- ave(contrast_results$p_value, contrast_results$contrast, FUN = function(x) p.adjust(x, 'BH'))

# ----------------------------- PSPS correlations ------------------------------
cor_results <- bind_rows(lapply(score_names, function(nm) {
  ct <- suppressWarnings(cor.test(sample_scores$PSPS_core, sample_scores[[nm]], method = 'spearman', exact = FALSE))
  tibble(test_type = 'PSPS_correlation', contrast = 'all_12wk_samples', score = nm,
         estimate = unname(ct$estimate), rho = unname(ct$estimate),
         p_value = ct$p.value, n = nrow(sample_scores))
}))
cor_results$q_within_family <- p.adjust(cor_results$p_value, 'BH')

# All four scores, both primary contrasts and the PSPS correlation are retained.
all_tests <- bind_rows(
  cor_results %>% transmute(test_type, contrast, score, estimate, p_value),
  contrast_results %>% transmute(test_type, contrast, score, estimate, p_value)
)
all_tests$q_global_external <- p.adjust(all_tests$p_value, 'BH')
cor_results <- cor_results %>% left_join(all_tests %>% filter(test_type == 'PSPS_correlation') %>% select(test_type, contrast, score, q_global_external), by = c('test_type', 'contrast', 'score'))
contrast_results <- contrast_results %>% left_join(all_tests %>% filter(test_type == 'dietary_contrast') %>% select(test_type, contrast, score, q_global_external), by = c('test_type', 'contrast', 'score'))
write.csv(cor_results, file.path(outdir, '12wk_psps_external_lipid_correlations.csv'), row.names = FALSE)
write.csv(contrast_results, file.path(outdir, '12wk_external_lipid_contrasts.csv'), row.names = FALSE)
write.csv(all_tests, file.path(outdir, '12wk_external_lipid_all_tests_fdr.csv'), row.names = FALSE)

# Pairwise Spearman matrix among the frozen PSPS and all external lipid scores.
matrix_scores <- as.data.frame(sample_scores[, c('PSPS_core', score_names), drop = FALSE])
pairs <- expand.grid(score1 = colnames(matrix_scores), score2 = colnames(matrix_scores), stringsAsFactors = FALSE) %>%
  filter(score1 < score2)
pairwise <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  a <- matrix_scores[[pairs$score1[i]]]; b <- matrix_scores[[pairs$score2[i]]]
  ct <- suppressWarnings(cor.test(a, b, method = 'spearman', exact = FALSE))
  tibble(score1 = pairs$score1[i], score2 = pairs$score2[i], n = length(a),
         spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}))
pairwise$FDR_BH <- p.adjust(pairwise$p_value, 'BH')
write.csv(pairwise, file.path(outdir, '12wk_psps_external_lipid_pairwise_spearman.csv'), row.names = FALSE)

# ----------------------------- common panel -----------------------------------
label_map <- c(
  Hallmark_FA_Metabolism = 'Hallmark FA metabolism',
  Hallmark_Adipogenesis = 'Hallmark adipogenesis',
  KEGG_PPAR_Signaling = 'KEGG PPAR signalling',
  Reactome_FA_Metabolism = 'Reactome FA metabolism'
)
contrast_labels <- c(primary__HO_Safflower_vs_Safflower = 'HO safflower vs safflower',
                     primary__HO_Sunflower_vs_Sunflower = 'HO sunflower vs sunflower')
cor_plot <- cor_results %>% transmute(score, score_label = unname(label_map[score]),
  analysis = 'PSPS correlation', estimate, label_x = estimate + 0.04, text_hjust = 0, lower = NA_real_, upper = NA_real_,
  q = q_global_external, role = spec$role[match(score, spec$score)])
contrast_plot <- contrast_results %>% transmute(score, score_label = unname(label_map[score]),
  analysis = unname(contrast_labels[contrast]), estimate, label_x = estimate, text_hjust = 0.5, lower = CI95_lower, upper = CI95_upper,
  q = q_global_external, role = spec$role[match(score, spec$score)])
panel_df <- bind_rows(cor_plot, contrast_plot) %>%
  mutate(score_label = factor(score_label, levels = rev(unname(label_map))),
         analysis = factor(analysis, levels = c('PSPS correlation', 'HO safflower vs safflower', 'HO sunflower vs sunflower')))

# A single shared panel reports the PSPS correlation and both dietary contrasts
# for all four pre-specified scores.
panel <- ggplot(panel_df, aes(y = score_label, x = estimate, colour = role)) +
  geom_vline(xintercept = 0, linetype = 2, colour = 'grey55') +
  geom_segment(data = filter(panel_df, is.finite(lower)),
               aes(x = lower, xend = upper, yend = score_label), linewidth = .7) +
  geom_point(size = 2.5) +
  geom_text(aes(x = label_x, hjust = text_hjust, label = paste0('q=', formatC(q, format = 'e', digits = 2))),
            vjust = -0.9, size = 2.7, show.legend = FALSE) +
  facet_grid(. ~ analysis, scales = 'free_x') +
  scale_colour_manual(values = c(primary = '#0072B2', sensitivity = '#D55E00')) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
        axis.title.y = element_blank(), strip.background = element_rect(fill = 'grey92')) +
  labs(x = 'Spearman rho or ssGSEA score difference (95% CI)',
       colour = 'Pre-specified role',
       title = 'External lipid-program validation of frozen PSPS')
ggsave(file.path(figdir, '12wk_external_lipid_validation_panel.pdf'), panel, width = 13, height = 5.8, useDingbats = FALSE)
ggsave(file.path(figdir, '12wk_external_lipid_validation_panel.png'), panel, width = 13, height = 5.8, dpi = 220)

# ----------------------------- manifest and session ---------------------------
manifest <- c(
  'Analysis: external lipid-program validation of the unchanged six-gene PSPS.',
  paste0('Input: GSE292815 KC 12-week normalized matrix; samples=', ncol(expr)),
  paste0('Frozen PSPS genes: ', paste(required_core, collapse = ', ')),
  'External scores: GSVA 2.x ssGSEA (alpha=0.25, normalize=TRUE) on four pre-specified sets.',
  'Primary external sets: MSigDB Hallmark Fatty Acid Metabolism and KEGG PPAR signaling.',
  'Sensitivity sets: MSigDB Hallmark Adipogenesis and Reactome Fatty Acid Metabolism.',
  'Mouse scoring uses msigdbr human MSigDB definitions with msigdbr ortholog mapping to Mus musculus.',
  'All four scores are reported regardless of direction or significance; no score was selected after seeing results.',
  'Planned contrasts: HO safflower versus safflower and HO sunflower versus sunflower.',
  'FDR: q_within_family is BH within PSPS correlations (n=4) or each dietary contrast (n=4); q_global_external is the pre-specified BH family across all 12 external tests and is used in the figure.',
  'The existing PSPS score was recomputed from exactly six genes and required to match the prior output to 1e-10.',
  'This external module does not merge lipid genes into PSPS and does not alter existing robustness analyses.'
)
writeLines(manifest, file.path(outdir, 'README.md'))
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo_external_lipid.txt'))
message('Completed external lipid validation: ', outdir)
