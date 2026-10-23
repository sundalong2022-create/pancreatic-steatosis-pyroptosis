#!/usr/bin/env Rscript
# AUCell cutoff audit for the fixed six-gene PRII in 12-week GSE292815 bulk samples.
# AUCell's conventional 5% cutoff is the prespecified reference; broader cutoffs
# are diagnostic sensitivity analyses and are not selected by correlation.

suppressPackageStartupMessages({
  library(AUCell)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
out_dir <- file.path(root, "results", "supplementary", "S3_prii_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(root, "results", "qc", "GSE292815_KC_12wk_matrix.rds"))
expr <- as.matrix(obj$expr)
core <- readRDS(file.path(root, "scripts", "results", "gene_sets_mouse", "gene_sets.rds"))$core_axis
sample_scores <- read.csv(file.path(root, "results", "12wk_primary", "scores", "12wk_sample_scores.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(is.numeric(expr), ncol(expr) == 21L, !anyDuplicated(rownames(expr)),
          all(core %in% rownames(expr)), all(colnames(expr) %in% sample_scores$sample))
sample_scores <- sample_scores[match(colnames(expr), sample_scores$sample), , drop = FALSE]
stopifnot(identical(as.character(sample_scores$sample), colnames(expr)),
          all(is.finite(sample_scores$PSPS_core)))

seed <- 20261011L
fractions <- c(0.05, 0.10, 0.20, 0.50, 0.80)
set.seed(seed) # AUCell breaks tied expression ranks randomly.
rankings_obj <- AUCell_buildRankings(expr, featureType = "genes", plotStats = FALSE,
                                     keepZeroesAsNA = FALSE, verbose = FALSE)
rankings <- SummarizedExperiment::assays(rankings_obj)[["ranking"]]
if (is.null(rankings) || !all(core %in% rownames(rankings))) stop("AUCell rankings are missing PRII genes.")

score_out <- data.frame(sample = colnames(expr), group = sample_scores$group,
                        PRII_mean_z = sample_scores$PSPS_core,
                        stringsAsFactors = FALSE, check.names = FALSE)
cutoff_rows <- vector("list", length(fractions))
gene_rows <- vector("list", length(fractions))

for (i in seq_along(fractions)) {
  fraction <- fractions[i]
  max_rank <- ceiling(fraction * nrow(expr))
  auc_obj <- AUCell_calcAUC(list(PRII = core), rankings_obj, aucMaxRank = max_rank,
                            normAUC = TRUE, verbose = FALSE)
  auc_mat <- AUCell::getAUC(auc_obj)
  auc <- as.numeric(auc_mat["PRII", colnames(expr)])
  score_col <- paste0("AUCell_top", sprintf("%02d", round(fraction * 100)), "pct")
  score_out[[score_col]] <- auc
  rank_core <- as.matrix(rankings[core, colnames(expr), drop = FALSE])
  in_window <- rank_core <= max_rank
  coverage_by_gene <- rowSums(in_window)
  cor_test <- suppressWarnings(cor.test(sample_scores$PSPS_core, auc,
                                        method = "spearman", exact = FALSE))
  cutoff_rows[[i]] <- data.frame(
    cutoff_fraction = fraction, cutoff_label = paste0(round(fraction * 100), "%"),
    aucMaxRank = max_rank, n_samples = ncol(expr), n_core_genes = length(core),
    n_genes_in_any_sample = sum(coverage_by_gene > 0L),
    genes_in_any_sample = paste(core[coverage_by_gene > 0L], collapse = ";"),
    n_genes_in_all_samples = sum(coverage_by_gene == ncol(expr)),
    genes_in_all_samples = paste(core[coverage_by_gene == ncol(expr)], collapse = ";"),
    min_genes_per_sample = min(colSums(in_window)),
    median_genes_per_sample = median(colSums(in_window)),
    max_genes_per_sample = max(colSums(in_window)),
    spearman_rho = unname(cor_test$estimate), p_value = cor_test$p.value,
    auc_min = min(auc), auc_max = max(auc),
    interpretation = if (fraction == 0.05 && sum(coverage_by_gene == ncol(expr)) < length(core))
      "Conventional 5% cutoff; incomplete six-gene coverage, not interpretable as full PRII agreement" else
      "Cutoff sensitivity only; do not select based on correlation",
    stringsAsFactors = FALSE
  )
  gene_rows[[i]] <- data.frame(
    cutoff_fraction = fraction, cutoff_label = paste0(round(fraction * 100), "%"),
    aucMaxRank = max_rank, gene = core, n_samples_in_AUC_window = as.integer(coverage_by_gene),
    total_samples = ncol(expr), fraction_samples_in_AUC_window = coverage_by_gene / ncol(expr),
    min_rank = apply(rank_core, 1, min), median_rank = apply(rank_core, 1, median),
    max_rank = apply(rank_core, 1, max), stringsAsFactors = FALSE
  )
}

cutoff_summary <- do.call(rbind, cutoff_rows)
gene_coverage <- do.call(rbind, gene_rows)
write.csv(score_out, file.path(out_dir, "S3_AUCell_cutoff_sensitivity_scores.csv"), row.names = FALSE)
write.csv(cutoff_summary, file.path(out_dir, "S3_AUCell_cutoff_sensitivity_summary.csv"), row.names = FALSE)
write.csv(gene_coverage, file.path(out_dir, "S3_AUCell_gene_coverage_by_cutoff.csv"), row.names = FALSE)

# Record the conventional AUCell 5% result in the implementation-agreement source
# while explicitly marking it non-comparable when the full signature is not covered.
agreement_path <- file.path(root, "results", "figure3_supplementary_tables",
                            "12wk_scoring_implementation_agreement.csv")
agreement <- read.csv(agreement_path, stringsAsFactors = FALSE, check.names = FALSE)
ix <- which(agreement$score_alt == "AUCell")
stopifnot(length(ix) == 1L)
primary <- cutoff_summary[cutoff_summary$cutoff_fraction == 0.05, , drop = FALSE]
agreement$n[ix] <- primary$n_samples
agreement$method[ix] <- "spearman"
agreement$rho[ix] <- primary$spearman_rho
agreement$p_value[ix] <- primary$p_value
agreement$scale[ix] <- "rank"
agreement$family_id[ix] <- "none"
agreement$notes[ix] <- paste0(
  "computed with conventional aucMaxRank=5% (", primary$aucMaxRank, " genes); ",
  primary$n_genes_in_all_samples, "/", length(core),
  " PRII genes entered the AUC window in every sample (only ",
  primary$genes_in_all_samples, "); rho is not comparable as a full-signature score. " ,
  "See S3_AUCell_cutoff_sensitivity_summary.csv. AUCell ", as.character(packageVersion("AUCell")), "."
)
write.csv(agreement, agreement_path, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_AUCell.txt"))
manifest <- c(
  "AUCell 12-week bulk PRII cutoff audit",
  paste0("Matrix: 12-week GSE292815 KC expression matrix; ", nrow(expr), " genes x ", ncol(expr), " samples."),
  paste0("PRII genes: ", paste(core, collapse = ", "), "."),
  "Input: the existing non-negative expression matrix used by the frozen score pipeline; AUCell ranks each sample's genes.",
  "Ties were handled by AUCell_buildRankings with a fixed seed; cutoff sensitivity used 5%, 10%, 20%, 50%, and 80% of ranked genes.",
  "The conventional AUCell reference cutoff is 5%. Broader cutoffs are diagnostic only and were not selected by correlation.",
  paste0("Seed: ", seed, ". AUCell version: ", as.character(packageVersion("AUCell")), "."),
  "The 5% score is not interpreted as full PRII agreement if all six genes do not enter the AUC window in each sample."
)
writeLines(manifest, file.path(out_dir, "AUCell_run_manifest.txt"))

cat("AUCell cutoff audit written to ", out_dir, "\n", sep = "")
print(cutoff_summary[, c("cutoff_label", "aucMaxRank", "n_genes_in_any_sample",
                         "n_genes_in_all_samples", "min_genes_per_sample",
                         "max_genes_per_sample", "spearman_rho", "p_value")], row.names = FALSE)
