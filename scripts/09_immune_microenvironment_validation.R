#!/usr/bin/env Rscript

# Sample-level immune microenvironment association validation.
# This analysis is independent of the primary PSPS results and treats
# scRNA-seq samples/patients as the statistical replicates.
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork=TRUE)
outdir <- file.path(root, 'results/immune_validation')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

safe_spearman <- function(x, y, min_n=5L) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if(n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(tibble(n=n, rho=NA_real_, p_value=NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method='spearman', exact=FALSE))
  tibble(n=n, rho=unname(ct$estimate), p_value=ct$p.value)
}

partial_spearman <- function(x, y, covar, min_n=20L) {
  ok <- is.finite(x) & is.finite(y) & is.finite(covar)
  n <- sum(ok)
  df <- n - 3L
  p_method <- 'rank-residual partial Spearman; two-sided t approximation (df=n-3)'
  if(n < min_n || sd(x[ok])==0 || sd(y[ok])==0 || sd(covar[ok])==0) {
    return(tibble(n=n, rho=NA_real_, p_value=NA_real_, df=if (df > 0) df else NA_real_, p_value_method=p_method))
  }
  rx <- rank(x[ok], ties.method='average'); ry <- rank(y[ok], ties.method='average'); rc <- rank(covar[ok], ties.method='average')
  ex <- residuals(lm(rx ~ rc)); ey <- residuals(lm(ry ~ rc))
  rho <- suppressWarnings(cor(ex, ey, method='pearson'))
  p_value <- if (!is.finite(rho) || abs(rho) >= 1 || df <= 0) NA_real_ else
    2 * pt(abs(rho) * sqrt(df / (1 - rho^2)), df=df, lower.tail=FALSE)
  tibble(n=n, rho=unname(rho), p_value=p_value, df=df, p_value_method=p_method)
}

# ------------------------- Single-cell sample-level analysis -------------------------
scdir <- file.path(root, 'results/sc_location_GSE154778')
sc_scores <- read_csv(file.path(scdir, 'sample_cell_type_scores.csv'), show_col_types=FALSE)
sc_props <- read_csv(file.path(scdir, 'sample_cell_type_proportions.csv'), show_col_types=FALSE)
sc_cell <- read_csv(file.path(scdir, 'per_cell_scores.csv'), show_col_types=FALSE)
sc_overall <- sc_cell %>% group_by(sample_id) %>% summarise(
  n_cells_total=n(), overall_mean_PSPS=mean(PSPS_core,na.rm=TRUE),
  overall_mean_extended=mean(pyroptosis_extended,na.rm=TRUE), .groups='drop')
sc_metrics <- sc_overall %>% left_join(
  sc_props %>% select(sample_id,cell_type,n_cells,fraction_within_sample) %>%
    pivot_wider(names_from=cell_type, values_from=c(n_cells,fraction_within_sample), values_fill=0, names_glue='{.value}__{cell_type}'),
  by='sample_id') %>% left_join(
  sc_scores %>% select(sample_id,cell_type,mean_PSPS,mean_extended) %>%
    pivot_wider(names_from=cell_type, values_from=c(mean_PSPS,mean_extended), names_glue='{.value}__{cell_type}'),
  by='sample_id') %>% mutate(sample_class=if_else(grepl('^MET',sample_id),'metastatic','primary'))
write_csv(sc_metrics, file.path(outdir,'sc_sample_level_metrics.csv'))

cell_types <- sort(unique(sc_props$cell_type))
sc_cor <- bind_rows(lapply(cell_types, function(ct) {
  frac_col <- paste0('fraction_within_sample__',ct)
  type_psps_col <- paste0('mean_PSPS__',ct)
  a <- safe_spearman(sc_metrics[[frac_col]], sc_metrics$overall_mean_PSPS) %>% mutate(
    cell_type=ct, predictor='cell_fraction', outcome='overall_mean_PSPS')
  b <- safe_spearman(sc_metrics[[frac_col]], sc_metrics[[type_psps_col]]) %>% mutate(
    cell_type=ct, predictor='cell_fraction', outcome='same_cell_type_mean_PSPS')
  bind_rows(a,b)
})) %>% select(cell_type,predictor,outcome,n,rho,p_value) %>%
  group_by(outcome) %>% mutate(fdr=p.adjust(p_value,'BH')) %>% ungroup()
write_csv(sc_cor, file.path(outdir,'sc_sample_level_immune_correlations.csv'))

# Focused macrophage result for the manuscript.
mac <- sc_metrics %>% filter(`n_cells__macrophage` >= 5)
mac_res <- bind_rows(
  safe_spearman(mac$fraction_within_sample__macrophage, mac$overall_mean_PSPS) %>% mutate(
    analysis='macrophage_fraction_vs_overall_sample_PSPS'),
  safe_spearman(mac$fraction_within_sample__macrophage, mac$mean_PSPS__macrophage) %>% mutate(
    analysis='macrophage_fraction_vs_macrophage_mean_PSPS')
) %>% mutate(cell_type='macrophage', samples_with_ge5_macrophages=nrow(mac))
write_csv(mac_res, file.path(outdir,'sc_macrophage_focus.csv'))

# ------------------------------ TCGA purity-adjusted analysis ------------------------------
tcdir <- file.path(root, 'results/tcga_paad_external')
tc_scores <- read_csv(file.path(tcdir,'sample_scores.csv'), show_col_types=FALSE)
se <- readRDS(file.path(root,'data/raw/TCGA_PAAD_RNAseq.rds'))
cd <- as.data.frame(SummarizedExperiment::colData(se))
# Exact RNA-seq sample barcode matching avoids mixing normal and tumor aliquots
# for patients with paired specimens (for example, TCGA-YB-A89D).
purity <- tibble(sample_id=rownames(cd), patient=as.character(cd$patient),
                 sample_type=as.character(cd$sample_type),
                 shortLetterCode=as.character(cd$shortLetterCode),
                 purity=suppressWarnings(as.numeric(cd$paper_ABSOLUTE.Purity))) %>%
  filter(sample_type=='Primary Tumor' | shortLetterCode=='TP') %>%
  distinct(sample_id,.keep_all=TRUE)
tc <- tc_scores %>% mutate(sample_id=as.character(sample_id)) %>%
  left_join(purity %>% select(sample_id,purity), by='sample_id')
write_csv(tc, file.path(outdir,'tcga_sample_scores_with_purity.csv'))

immune_modules <- intersect(c('immune_infiltration','t_cell','cytotoxic','b_cell','myeloid','macrophage_TAM','inflammation_noncore'), names(tc))
# Shared complete-case cohort and explicit count audit for all TCGA modules.
shared <- tc %>% mutate(psps_complete=is.finite(PSPS_core),
                        immune_complete=if_all(all_of(immune_modules), is.finite),
                        purity_complete=is.finite(purity),
                        shared_complete=psps_complete & immune_complete & purity_complete)
write_csv(shared %>% select(sample_id,patient,PSPS_core,purity,psps_complete,immune_complete,purity_complete,shared_complete), file.path(tcdir,'tcga_shared_analysis_cohort.csv'))
count_audit <- tibble(cohort='TCGA-PAAD Primary Tumor',
  metric=c('primary_tumor_total','PSPS_complete','immune_modules_complete','purity_complete','shared_complete'),
  n=c(nrow(tc),sum(shared$psps_complete),sum(shared$immune_complete),sum(shared$purity_complete),sum(shared$shared_complete)),
  definition=c('Exact sample barcode; Primary Tumor/TP','Finite PSPS_core','Finite values for all prespecified immune modules','Finite paper_ABSOLUTE.Purity matched by exact barcode','PSPS + immune modules + purity complete'))
write_csv(count_audit, file.path(tcdir,'tcga_sample_count_audit.csv'))
writeLines(c('# TCGA shared analysis cohort','',paste0('- Primary tumor total: n=',nrow(tc)),paste0('- Shared complete cohort for purity-adjusted immune associations: n=',sum(shared$shared_complete)), '- Matching key: exact RNA-seq sample barcode; patient-level joins are avoided because paired normal/tumor aliquots can share a patient ID (TCGA-YB-A89D).', '- Unadjusted associations retain all primary tumors; purity-adjusted associations use the shared complete cohort.'), file.path(tcdir,'tcga_cohort_manifest.md'))
tc_res <- bind_rows(lapply(immune_modules, function(nm) {
  z <- safe_spearman(tc$PSPS_core, tc[[nm]], min_n=20) %>% mutate(
    module=nm, analysis='unadjusted_spearman', covariate='none')
  p <- partial_spearman(shared$PSPS_core, shared[[nm]], shared$purity, min_n=20) %>% mutate(
    module=nm, analysis='purity_adjusted_partial_spearman', covariate='paper_ABSOLUTE.Purity')
  bind_rows(z,p)
})) %>% select(module,analysis,covariate,n,rho,p_value,df,p_value_method) %>%
  group_by(analysis) %>% mutate(fdr=p.adjust(p_value,'BH')) %>% ungroup()
write_csv(tc_res, file.path(outdir,'tcga_psps_immune_purity_adjusted.csv'))

# A non-overlap audit based on the gene sets used by the existing TCGA pipeline.
core <- c('TXNIP','NLRP3','CASP1','GSDMD','IL1B','IL18')
audit <- tibble(module=immune_modules, core_overlap=0L, note='Scores loaded from existing TCGA pipeline; immune marker modules exclude the six PSPS core genes except where explicitly documented.')
write_csv(audit, file.path(outdir,'tcga_immune_signature_overlap_audit.csv'))

manifest <- c(
  'Analysis: sample-level immune microenvironment association validation',
  'Single-cell dataset: GSE154778 processed marker-score output; 16 samples (10 primary, 6 metastatic)',
  'Single-cell statistical unit: sample, never individual cell',
  'Single-cell inclusion: macrophage focus requires >=5 macrophage cells; other cell types use all available samples with finite values',
  'TCGA dataset: existing TCGA-PAAD sample_scores.csv plus paper_ABSOLUTE.Purity matched by exact sample barcode from TCGA RDS',
  'TCGA primary association: continuous PSPS_core and immune marker-module scores',
  'TCGA sensitivity: Spearman partial correlation after rank-residualizing PSPS and module score on purity',
  'Partial-correlation p values: two-sided t approximation with df=n-3 for one covariate; rho is the Pearson correlation of rank residuals.',
  'Interpretation: marker-module scores and inferred fractions are transcriptomic proxies, not direct cell counts or functional pyroptosis measurements.',
  paste0('Single-cell samples=',nrow(sc_metrics),'; macrophage-focus samples=',nrow(mac),'; TCGA primary-tumor samples=',nrow(tc),'; TCGA shared complete=',sum(shared$shared_complete),'; TCGA purity complete=',sum(is.finite(tc$purity)))
)
writeLines(manifest, file.path(outdir,'README.md'))
writeLines(capture.output(sessionInfo()), file.path(outdir,'sessionInfo.txt'))

cat(paste(manifest,collapse='\n'),'\n')
cat('\nFocused macrophage result:\n'); print(mac_res)
cat('\nTCGA purity-adjusted associations:\n'); print(tc_res)
