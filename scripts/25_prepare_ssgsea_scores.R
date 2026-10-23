#!/usr/bin/env Rscript
# Prepare the independent six-gene GSVA/ssGSEA scores for Figure 2e-iii.
# Scores are reused from the frozen 12-week robustness analysis.

root <- normalizePath('.', mustWork = TRUE)
source_file <- file.path(root, 'results/12wk_primary/robustness/12wk_ssgsea_scores.csv')
source <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c('sample', 'PSPS_core_ssGSEA') %in% names(source)))
stopifnot(nrow(source) == 21L, !anyDuplicated(source$sample))

# Match the compact IDs used in fig2f_six_gene_z_long.csv.
sample_id <- gsub('_M_', '_', source$sample, fixed = TRUE)
sample_id <- sub('^HO_Safflower_', 'HO-Safflower_', sample_id)
sample_id <- sub('^HO_Sunflower_', 'HO-Sunflower_', sample_id)
out <- data.frame(sample = sample_id, ssGSEA = source$PSPS_core_ssGSEA, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(out$sample), all(is.finite(out$ssGSEA)))

heatmap_file <- file.path(root, 'results/12wk_primary/scores/fig2f_six_gene_z_long.csv')
heatmap_samples <- unique(read.csv(heatmap_file, stringsAsFactors = FALSE)$sample)
stopifnot(setequal(out$sample, heatmap_samples))

out_file <- file.path(root, 'results/12wk_primary/scores/ssgsea_scores.csv')
write.csv(out, out_file, row.names = FALSE, quote = TRUE)
message(sprintf('Wrote %d ssGSEA scores aligned to the six-gene z table: %s', nrow(out), out_file))
