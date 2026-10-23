#!/usr/bin/env Rscript

# Formal bulk deconvolution for TCGA-PAAD.
# xCell, EPIC, ESTIMATE and a local CIBERSORT-compatible LM22 SVR run.
# Official CIBERSORTx is a restricted web/Docker service and is tracked
# separately in deconvolution_status.txt.

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(EPIC)
  library(xCell)
  library(csCIBERSORTx)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
outdir <- file.path(root, 'results/tcga_paad_external/deconvolution')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

spearman <- function(x, y, min_n = 20L) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(n = n, rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = 'spearman', exact = FALSE))
  tibble(n = n, rho = unname(ct$estimate), p_value = ct$p.value)
}
partial_spearman <- function(x, y, z, min_n = 20L) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z); n <- sum(ok)
  df <- n - 3L
  p_method <- 'rank-residual partial Spearman; two-sided t approximation (df=n-3)'
  if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0 || sd(z[ok]) == 0) {
    return(tibble(n = n, rho = NA_real_, p_value = NA_real_, df = if (df > 0) df else NA_real_, p_value_method = p_method))
  }
  ex <- residuals(lm(rank(x[ok], ties.method = 'average') ~ rank(z[ok], ties.method = 'average')))
  ey <- residuals(lm(rank(y[ok], ties.method = 'average') ~ rank(z[ok], ties.method = 'average')))
  rho <- suppressWarnings(cor(ex, ey, method = 'pearson'))
  p_value <- if (!is.finite(rho) || abs(rho) >= 1 || df <= 0) NA_real_ else
    2 * pt(abs(rho) * sqrt(df / (1 - rho^2)), df = df, lower.tail = FALSE)
  tibble(n = n, rho = unname(rho), p_value = p_value, df = df, p_value_method = p_method)
}
extract_cols <- function(df, mapping) {
  out <- tibble(sample_id = as.character(df$sample_id))
  for (nm in names(mapping)) {
    hit <- intersect(mapping[[nm]], names(df))
    out[[nm]] <- if (length(hit)) suppressWarnings(as.numeric(df[[hit[1]]])) else NA_real_
  }
  out
}

# --------------------------- Input and PSPS ---------------------------
se <- readRDS(file.path(root, 'data/raw/TCGA_PAAD_RNAseq.rds'))
if (!inherits(se, 'SummarizedExperiment')) stop('TCGA input is not SummarizedExperiment')
expr0 <- assay(se, 'tpm_unstrand')
rd <- as.data.frame(rowData(se))
gene <- toupper(as.character(rd$gene_name))
gene[is.na(gene) | gene == ''] <- rownames(expr0)[is.na(gene) | gene == '']
gene <- sub('\\..*$', '', gene)
keep_gene <- !duplicated(gene)
expr <- expr0[keep_gene, , drop = FALSE]
rownames(expr) <- gene[keep_gene]
cd <- as.data.frame(colData(se))
tumor <- cd$sample_type == 'Primary Tumor' | cd$shortLetterCode == 'TP'
tumor[is.na(tumor)] <- FALSE
expr <- as.matrix(expr[, tumor, drop = FALSE])
cd <- cd[colnames(expr), , drop = FALSE]
expr[!is.finite(expr) | expr < 0] <- 0

scores <- read_csv(file.path(root, 'results/tcga_paad_external/sample_scores.csv'), show_col_types = FALSE)
scores$sample_id <- as.character(scores$sample_id)
if (!setequal(colnames(expr), scores$sample_id)) stop('Expression and PSPS sample IDs do not match')
scores <- scores[match(colnames(expr), scores$sample_id), , drop = FALSE]

# Upload-ready files for the official CIBERSORTx interface. The official
# service accepts a gene-symbol-by-sample mixture; TPM is retained as the
# primary input and a log2(TPM+1) copy is provided for audit/comparison.
write.table(cbind(GeneSymbol = rownames(expr), expr),
  file = file.path(outdir, 'CIBERSORTx_input_TPM.tsv'), sep = '\t',
  quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cbind(GeneSymbol = rownames(expr), log2(expr + 1)),
  file = file.path(outdir, 'CIBERSORTx_input_log2TPM.tsv'), sep = '\t',
  quote = FALSE, row.names = FALSE, col.names = TRUE)
writeLines(c(
  paste0('input_samples=', ncol(expr)), paste0('input_genes=', nrow(expr)),
  'primary_input=CIBERSORTx_input_TPM.tsv',
  'secondary_input=CIBERSORTx_input_log2TPM.tsv',
  'sample identifiers are exact TCGA RNA-seq barcodes',
  'official CIBERSORTx upload still requires a registered CIBERSORTx account/token.'
), file.path(outdir, 'CIBERSORTx_input_manifest.txt'))

# ------------------------------- xCell --------------------------------
message('Running xCell ...')
x <- xCell::xCellAnalysis(expr, rnaseq = TRUE, parallel.sz = 2, parallel.type = 'SOCK')
if (!is.matrix(x)) stop('xCell did not return a matrix')
x_tbl <- as.data.frame(t(x), check.names = FALSE) %>% tibble::rownames_to_column('sample_id')
write_csv(x_tbl, file.path(outdir, 'xcell_scores.csv'))

# -------------------------------- EPIC ---------------------------------
message('Running EPIC ...')
e <- EPIC::EPIC(expr, withOtherCells = TRUE)
e_tbl <- as.data.frame(e$cellFractions, check.names = FALSE) %>% tibble::rownames_to_column('sample_id')
if (!is.null(e$otherCells)) {
  other <- tibble(sample_id = rownames(e$otherCells), EPIC_otherCells = as.numeric(e$otherCells[, 1]))
  e_tbl <- left_join(e_tbl, other, by = 'sample_id')
}
write_csv(e_tbl, file.path(outdir, 'epic_fractions.csv'))
if (!is.null(e$fit.gof)) {
  gof <- as.data.frame(e$fit.gof, check.names = FALSE) %>% tibble::rownames_to_column('sample_id')
  write_csv(gof, file.path(outdir, 'epic_fit_gof.csv'))
}

# ----------------------- CIBERSORT-compatible LM22 ---------------------
message('Running local CIBERSORT-compatible LM22 SVR ...')
data(LM22, package = 'csCIBERSORTx')
cs <- csCIBERSORTx::CIBERSORT2(LM22_matrix = LM22, mixture = expr,
  opt_path = outdir, perm = 0, QN = FALSE)
cs_tbl <- as.data.frame(cs, check.names = FALSE) %>% tibble::rownames_to_column('sample_id')
write_csv(cs_tbl, file.path(outdir, 'cibersort_lm22_local.csv'))

# ------------------------------- ESTIMATE -------------------------------
message('Running ESTIMATE ...')
eref <- readRDS(file.path(root, 'data/reference/estimate_genes.rds'))
common <- eref$common_genes
sets <- eref$signature_geneset
idx <- match(common$GeneSymbol, rownames(expr)); ok <- !is.na(idx)
em <- expr[idx[ok], , drop = FALSE]; rownames(em) <- common$GeneSymbol[ok]
for (j in seq_len(ncol(em))) em[, j] <- rank(em[, j], ties.method = 'average')
em <- 10000 * em / nrow(em)
gs <- as.matrix(sets[, -1, drop = FALSE])
score <- matrix(NA_real_, nrow = 2, ncol = ncol(em), dimnames = list(c('StromalScore', 'ImmuneScore'), colnames(em)))
for (k in seq_len(2)) {
  gset <- as.character(gs[k, ]); gset <- gset[nzchar(gset)]
  ov <- intersect(gset, rownames(em))
  if (!length(ov)) next
  for (j in seq_len(ncol(em))) {
    gl <- order(em[, j], decreasing = TRUE); hit <- as.numeric(gl %in% match(ov, rownames(em)))
    miss <- 1 - hit; N <- length(gl); Nh <- sum(hit); Nm <- N - Nh
    w <- abs(em[gl, j])^0.25; sw <- sum(w[hit == 1])
    if (!sw || !Nm) next
    score[k, j] <- sum(cumsum(hit * w / sw - miss / Nm))
  }
}
est_score <- colSums(score, na.rm = TRUE)
est_purity <- cos(0.6049872018 + 0.0001467884 * est_score)
est_tbl <- tibble(sample_id = colnames(em), StromalScore = as.numeric(score[1, ]),
  ImmuneScore = as.numeric(score[2, ]), ESTIMATEScore = as.numeric(est_score),
  TumorPurity = as.numeric(est_purity))
est_tbl$TumorPurity[est_tbl$TumorPurity < 0] <- NA_real_
write_csv(est_tbl, file.path(outdir, 'estimate_scores.csv'))

# -------------------- Harmonised methods table ----------------------------
metric <- tibble(sample_id = scores$sample_id, PSPS_core = as.numeric(scores$PSPS_core),
  ABSOLUTE_Purity = suppressWarnings(as.numeric(cd$paper_ABSOLUTE.Purity)))
xc <- extract_cols(x_tbl, list(
  xCell_Macrophages = c('Macrophages', 'Macrophage'),
  xCell_Monocytes = c('Monocytes', 'Monocyte'),
  xCell_CD8_T = c('CD8+ T-cells', 'CD8_Tcells', 'CD8+T-cells'), xCell_B_cells = c('B-cells', 'B_cells', 'Bcells'),
  xCell_NK = c('NK cells', 'NK_cells', 'NKcells'), xCell_Fibroblasts = c('Fibroblasts', 'Fibroblast'),
  xCell_Endothelial = c('Endothelial cells', 'Endothelial_cells', 'Endothelial'),
  xCell_ImmuneScore = c('ImmuneScore'), xCell_StromaScore = c('StromaScore')))
xt_names <- grep('^(CD4|CD8).*T.cells$|^Tgd.cells$|^Tregs$', names(x_tbl), value = TRUE)
xc$xCell_T_cells_mean <- if (length(xt_names)) rowMeans(x_tbl[, xt_names, drop = FALSE], na.rm = TRUE) else NA_real_
metric <- left_join(metric, xc, by = 'sample_id')
metric <- left_join(metric, extract_cols(e_tbl, list(
  EPIC_Bcells = c('Bcells', 'B.cells'), EPIC_Macrophages = c('Macrophages'),
  EPIC_CD4_Tcells = c('CD4_Tcells', 'CD4 T cells'), EPIC_CD8_Tcells = c('CD8_Tcells', 'CD8 T cells'),
  EPIC_NKcells = c('NKcells', 'NK cells'), EPIC_CAFs = c('CAFs', 'Fibroblasts'),
  EPIC_Endothelial = c('Endothelial'), EPIC_otherCells = c('EPIC_otherCells', 'otherCells'))), by = 'sample_id')
metric <- left_join(metric, extract_cols(cs_tbl, list(
  CIBERSORT_M0 = c('Macrophages M0'), CIBERSORT_M1 = c('Macrophages M1'), CIBERSORT_M2 = c('Macrophages M2'),
  CIBERSORT_CD8_T = c('T cells CD8'), CIBERSORT_Treg = c('T cells regulatory (Tregs)'),
  CIBERSORT_B_naive = c('B cells naive'), CIBERSORT_NK = c('NK cells resting', 'NK cells activated'))), by = 'sample_id')
metric <- left_join(metric, est_tbl, by = 'sample_id')
write_csv(metric, file.path(outdir, 'deconvolution_harmonized_metrics.csv'))

# -------------------- PSPS associations and purity adjustment --------------
metric_names <- setdiff(names(metric), c('sample_id', 'PSPS_core', 'ABSOLUTE_Purity'))
assoc <- bind_rows(lapply(metric_names, function(nm) {
  raw <- spearman(metric$PSPS_core, metric[[nm]]) %>% mutate(metric = nm, analysis = 'unadjusted_spearman', covariate = 'none')
  adj <- partial_spearman(metric$PSPS_core, metric[[nm]], metric$ABSOLUTE_Purity) %>% mutate(metric = nm, analysis = 'partial_spearman_adjusted_for_ABSOLUTE_purity', covariate = 'paper_ABSOLUTE.Purity')
  bind_rows(raw, adj)
})) %>% select(metric, analysis, covariate, n, rho, p_value, df, p_value_method) %>% group_by(analysis) %>% mutate(fdr_BH = p.adjust(p_value, 'BH')) %>% ungroup()
write_csv(assoc, file.path(outdir, 'deconvolution_psps_associations.csv'))

concepts <- list(
  macrophage = c('xCell_Macrophages', 'EPIC_Macrophages', 'CIBERSORT_M0', 'CIBERSORT_M1', 'CIBERSORT_M2'),
  t_cell = c('xCell_T_cells_mean', 'xCell_CD8_T', 'EPIC_CD4_Tcells', 'EPIC_CD8_Tcells', 'CIBERSORT_CD8_T', 'CIBERSORT_Treg'),
  b_cell = c('xCell_B_cells', 'EPIC_Bcells', 'CIBERSORT_B_naive'),
  stromal = c('xCell_Fibroblasts', 'xCell_Endothelial', 'EPIC_CAFs', 'EPIC_Endothelial'),
  purity = c('ABSOLUTE_Purity', 'TumorPurity')
)
consistency <- bind_rows(lapply(names(concepts), function(concept) {
  cols <- intersect(concepts[[concept]], names(metric)); if (length(cols) < 2) return(tibble())
  bind_rows(lapply(combn(cols, 2, simplify = FALSE), function(pr) {
    s <- spearman(metric[[pr[1]]], metric[[pr[2]]])
    tibble(concept = concept, metric_1 = pr[1], metric_2 = pr[2], n = s$n, rho = s$rho, p_value = s$p_value)
  }))
})) %>% mutate(fdr_BH = p.adjust(p_value, 'BH'))
write_csv(consistency, file.path(outdir, 'deconvolution_method_consistency.csv'))

status <- c(
  'TCGA-PAAD formal deconvolution run completed.',
  paste0('Primary tumor samples: ', nrow(metric)),
  'xCell: xcell_scores.csv; enrichment scores, not direct cell fractions.',
  'EPIC: epic_fractions.csv; cell fractions and optional convergence diagnostics.',
  'ESTIMATE: estimate_scores.csv; immune/stromal/ESTIMATE scores and transformed TumorPurity.',
  'Local CIBERSORT-compatible LM22: cibersort_lm22_local.csv; permutations=0.',
  'Official CIBERSORTx web/Docker service was not called because its restricted interface/credentials were unavailable.',
  'The local LM22 output must not be labelled as official CIBERSORTx output.',
  'ABSOLUTE purity is matched by exact TCGA RNA-seq sample barcode.',
  'Partial-correlation p values: two-sided t approximation with df=n-3 for one covariate; rho is the Pearson correlation of rank residuals.',
  'Associations: deconvolution_psps_associations.csv; method consistency: deconvolution_method_consistency.csv.'
)
writeLines(status, file.path(outdir, 'deconvolution_status.txt'))
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
cat(paste(status, collapse = '\n'), '\n')
cat('\nSelected PSPS associations:\n')
print(assoc %>% filter(metric %in% c('xCell_Macrophages', 'EPIC_Macrophages', 'CIBERSORT_M0', 'CIBERSORT_M1', 'CIBERSORT_M2', 'ImmuneScore', 'StromalScore', 'TumorPurity', 'ABSOLUTE_Purity')))
