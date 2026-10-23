#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(tibble)})
root <- normalizePath('.', mustWork = TRUE)
dir <- file.path(root, 'results/tcga_paad_external/deconvolution')
metric <- read_csv(file.path(dir, 'deconvolution_harmonized_metrics.csv'), show_col_types = FALSE)
gof <- read_csv(file.path(dir, 'epic_fit_gof.csv'), show_col_types = FALSE)
good <- gof %>% filter(is.finite(convergeCode), convergeCode == 0) %>% pull(sample_id)
bad <- setdiff(gof$sample_id, good)
spearman <- function(x, y, z = NULL) {
  ok <- is.finite(x) & is.finite(y); if (!is.null(z)) ok <- ok & is.finite(z)
  n <- sum(ok); df <- if (is.null(z)) NA_integer_ else n - 3L
  p_method <- if (is.null(z)) 'Spearman cor.test (exact=FALSE)' else
    'rank-residual partial Spearman; two-sided t approximation (df=n-3)'
  if (n < 20) return(tibble(n=n, rho=NA_real_, p_value=NA_real_, df=df, p_value_method=p_method))
  if (is.null(z)) {
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method='spearman', exact=FALSE))
    rho <- unname(ct$estimate); p_value <- ct$p.value
  } else {
    ex <- residuals(lm(rank(x[ok], ties.method='average') ~ rank(z[ok], ties.method='average')))
    ey <- residuals(lm(rank(y[ok], ties.method='average') ~ rank(z[ok], ties.method='average')))
    rho <- suppressWarnings(cor(ex, ey, method='pearson'))
    p_value <- if (!is.finite(rho) || abs(rho) >= 1 || df <= 0) NA_real_ else
      2 * pt(abs(rho) * sqrt(df / (1 - rho^2)), df=df, lower.tail=FALSE)
  }
  tibble(n=n, rho=unname(rho), p_value=p_value, df=df, p_value_method=p_method)
}
sub <- metric %>% filter(sample_id %in% good)
epic_cols <- c('EPIC_Bcells','EPIC_Macrophages','EPIC_CD4_Tcells','EPIC_CD8_Tcells','EPIC_NKcells','EPIC_CAFs','EPIC_Endothelial','EPIC_otherCells')
res <- bind_rows(lapply(epic_cols, function(nm) {
  a <- spearman(sub$PSPS_core, sub[[nm]]) %>% mutate(metric=nm, analysis='EPIC_converged_unadjusted')
  b <- spearman(sub$PSPS_core, sub[[nm]], sub$ABSOLUTE_Purity) %>% mutate(metric=nm, analysis='EPIC_converged_partial_ABSOLUTE_purity')
  bind_rows(a,b)
})) %>% select(metric, analysis, n, rho, p_value, df, p_value_method) %>% group_by(analysis) %>% mutate(fdr_BH=p.adjust(p_value,'BH')) %>% ungroup()
res <- res %>% mutate(total_samples=nrow(metric), converged_samples=length(good), nonconverged_samples=length(bad), nonconverged_ids=paste(bad, collapse=';'))
write_csv(res, file.path(dir, 'epic_convergence_sensitivity.csv'))
writeLines(c(paste0('EPIC converged samples: ', length(good)), paste0('EPIC non-converged samples: ', length(bad)), 'Partial-correlation p values: two-sided t approximation with df=n-3 for one covariate; rho is the Pearson correlation of rank residuals.', paste0('Non-converged IDs: ', paste(bad, collapse=';'))), file.path(dir, 'epic_convergence_manifest.txt'))
print(res)
