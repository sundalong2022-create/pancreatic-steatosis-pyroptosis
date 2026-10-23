#!/usr/bin/env Rscript
# Export the frozen 6-gene x 21-sample row-wise z values as a long CSV.
# The row z scores match the GSE292815 12-week PSPS analysis: log2(normalized
# expression + 1), standardized once per gene across all 21 samples using n-1 SD.

root <- normalizePath('.', mustWork = TRUE)
obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr
meta <- obj$metadata
stopifnot(ncol(expr) == 21L, nrow(meta) == 21L, all(colnames(expr) == meta$matrix_column))

genes <- c('Txnip', 'Nlrp3', 'Casp1', 'Gsdmd', 'Il1b', 'Il18')
stopifnot(all(genes %in% rownames(expr)))
logexpr <- log2(expr + 1)
z_mat <- t(scale(t(logexpr[genes, , drop = FALSE])))
z_mat[!is.finite(z_mat)] <- 0

group_map <- c(
  CD_12wk = 'CD',
  Safflower_12wk = 'Safflower',
  HO_Safflower_12wk = 'HO-Safflower',
  Sunflower_12wk = 'Sunflower',
  HO_Sunflower_12wk = 'HO-Sunflower'
)
sample_group <- unname(group_map[as.character(meta$group)])
stopifnot(!anyNA(sample_group), length(unique(sample_group)) == 5L)
sample_id <- sub('_M_', '_', meta$matrix_column, fixed = TRUE)
sample_id <- sub('^HO_Safflower_', 'HO-Safflower_', sample_id)
sample_id <- sub('^HO_Sunflower_', 'HO-Sunflower_', sample_id)
stopifnot(!anyDuplicated(sample_id))

# Preserve matrix sample order and the established six-gene order.
out <- expand.grid(
  gene = genes,
  sample_index = seq_along(sample_id),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
out$sample <- sample_id[out$sample_index]
out$group <- sample_group[out$sample_index]
out$z <- as.numeric(z_mat[cbind(match(out$gene, genes), out$sample_index)])
out <- out[, c('sample', 'group', 'gene', 'z')]

# Structural and normalization checks before writing.
stopifnot(identical(names(out), c('sample', 'group', 'gene', 'z')), nrow(out) == 126L)
stopifnot(all(table(out$gene) == 21L), all(table(out$sample) == 6L))
stopifnot(all(vapply(split(out$z, out$gene), function(x) abs(mean(x)) < 1e-12 && abs(sd(x) - 1) < 1e-12, logical(1))))
stopifnot(setequal(unique(out$group), unname(group_map)))

out_file <- file.path(root, 'results/12wk_primary/scores/fig2f_six_gene_z_long.csv')
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_file, row.names = FALSE, quote = TRUE, na = '')
message(sprintf('Wrote %d rows for %d samples and %d genes: %s', nrow(out), length(unique(out$sample)), length(unique(out$gene)), out_file))
