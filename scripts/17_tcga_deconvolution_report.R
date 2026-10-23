#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(tibble); library(ggplot2)})
root <- normalizePath('.', mustWork = TRUE)
dir <- file.path(root, 'results/tcga_paad_external/deconvolution')
assoc <- read_csv(file.path(dir, 'deconvolution_psps_associations.csv'), show_col_types = FALSE)
cons <- read_csv(file.path(dir, 'deconvolution_method_consistency.csv'), show_col_types = FALSE)
metrics <- read_csv(file.path(dir, 'deconvolution_harmonized_metrics.csv'), show_col_types = FALSE)

focus <- assoc %>% filter(metric %in% c('xCell_Macrophages','EPIC_Macrophages','CIBERSORT_M0','CIBERSORT_M1','CIBERSORT_M2','xCell_ImmuneScore','ImmuneScore','StromalScore','TumorPurity'))
write_csv(focus, file.path(dir, 'deconvolution_focus_associations.csv'))

plot_dat <- focus %>% filter(analysis != 'unadjusted_spearman') %>%
  mutate(metric = factor(metric, levels = rev(unique(metric))))
p <- ggplot(plot_dat, aes(x = rho, y = metric, colour = p_value < 0.05)) +
  geom_vline(xintercept = 0, linetype = 2, colour = 'grey50') + geom_point(size = 2.4) +
  facet_wrap(~analysis, scales = 'free_y') + scale_colour_manual(values = c('TRUE'='#D95F02','FALSE'='grey60')) +
  theme_classic(base_size = 10) + theme(legend.position='none', panel.grid.minor=element_blank()) +
  labs(title='TCGA-PAAD: PSPS associations with formal deconvolution metrics', x='Spearman rho (purity-adjusted where shown)', y=NULL)
ggsave(file.path(dir, 'deconvolution_psps_association_forest.pdf'), p, width=9, height=5.8, useDingbats=FALSE)

lines <- c(
  '# TCGA-PAAD formal deconvolution report', '',
  paste0('- Samples: ', nrow(metrics), ' Primary Tumor RNA-seq samples.'),
  '- xCell was run as an RNA-seq enrichment method; its values are enrichment scores rather than calibrated cell percentages.',
  '- EPIC was run with its built-in TRef reference and other-cell component; convergence diagnostics are provided in epic_fit_gof.csv.',
  '- ESTIMATE was run with the public ESTIMATE gene sets; TumorPurity is its transformed estimate and is separate from paper_ABSOLUTE.Purity.',
  '- The local CIBERSORT-compatible LM22 SVR result is a sensitivity analysis. It is not official CIBERSORTx.',
  '- Official CIBERSORTx requires the registered web/Docker service and credentials; upload-ready input files and a manifest are provided.', '',
  '- Purity-adjusted rank-residual partial Spearman P values use a two-sided t approximation with df=n-3 for one covariate; BH correction is applied separately within each full metric family.',
  '## Selected purity-adjusted findings', '',
  paste0('- xCell macrophage enrichment: rho=', sprintf('%.3f', assoc$rho[assoc$metric=='xCell_Macrophages' & grepl('partial',assoc$analysis)]), ', FDR=', formatC(assoc$fdr_BH[assoc$metric=='xCell_Macrophages' & grepl('partial',assoc$analysis)], format='e', digits=2), '.'),
  paste0('- EPIC macrophage fraction: rho=', sprintf('%.3f', assoc$rho[assoc$metric=='EPIC_Macrophages' & grepl('partial',assoc$analysis)]), ', FDR=', formatC(assoc$fdr_BH[assoc$metric=='EPIC_Macrophages' & grepl('partial',assoc$analysis)], format='e', digits=2), '.'),
  paste0('- ESTIMATE ImmuneScore: rho=', sprintf('%.3f', assoc$rho[assoc$metric=='ImmuneScore' & grepl('partial',assoc$analysis)]), ', FDR=', formatC(assoc$fdr_BH[assoc$metric=='ImmuneScore' & grepl('partial',assoc$analysis)], format='e', digits=2), '.'),
  paste0('- ESTIMATE TumorPurity versus PSPS: rho=', sprintf('%.3f', assoc$rho[assoc$metric=='TumorPurity' & grepl('partial',assoc$analysis)]), ', FDR=', formatC(assoc$fdr_BH[assoc$metric=='TumorPurity' & grepl('partial',assoc$analysis)], format='e', digits=2), '.'), '',
  'These results support an association between PSPS and an immune-rich bulk transcriptomic context. They do not identify a causal cell source or prove functional pyroptosis. TCGA-PAAD remains a validation of the candidate axis clinical context in pancreatic cancer, not evidence that TCGA proves fatty pancreas.'
)
writeLines(lines, file.path(dir, 'DECONVOLUTION_REPORT_zh.md'))
