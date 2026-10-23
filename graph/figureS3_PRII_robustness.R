#!/usr/bin/env Rscript
# Supplementary Figure S3: inferential stress tests for the fixed PRII score.
# Uses frozen 12-week robustness outputs only; it does not rerun or modify
# upstream analysis results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
robust_dir <- file.path(root, "results", "12wk_primary", "robustness")
out_dir <- file.path(root, "results", "supplementary", "S3_prii_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(name) {
  path <- file.path(robust_dir, name)
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

label_map <- c(
  primary__HO_Safflower_vs_Safflower = "HO-Safflower vs Safflower",
  primary__HO_Sunflower_vs_Sunflower = "HO-Sunflower vs Sunflower"
)
contrast_order <- unname(label_map)
group_a_n <- c("HO-Safflower vs Safflower" = 4L,
               "HO-Sunflower vs Sunflower" = 3L)
group_b_n <- c("HO-Safflower vs Safflower" = 5L,
               "HO-Sunflower vs Sunflower" = 4L)
contrast_n <- group_a_n + group_b_n
cols <- c("HO-Safflower vs Safflower" = "#B2182B",
          "HO-Sunflower vs Sunflower" = "#2166AC")

perm <- read_required("12wk_psps_label_permutation_summary.csv") %>%
  filter(contrast %in% names(label_map)) %>%
  mutate(contrast_label = unname(label_map[contrast]),
         contrast_label = factor(contrast_label, levels = rev(contrast_order)))
random <- read_required("12wk_random_six_gene_null_distribution.csv") %>%
  tidyr::pivot_longer(c(effect_HO_Safflower_vs_Safflower, effect_HO_Sunflower_vs_Sunflower),
                      names_to = "contrast", values_to = "effect") %>%
  mutate(contrast_label = recode(
    contrast,
    effect_HO_Safflower_vs_Safflower = "HO-Safflower vs Safflower",
    effect_HO_Sunflower_vs_Sunflower = "HO-Sunflower vs Sunflower"
  ), contrast_label = factor(contrast_label, levels = contrast_order))
random_summary_all <- read_required("12wk_random_six_gene_null_summary.csv")
random_summary <- random_summary_all %>%
  filter(test == "PSPS effect", contrast %in% names(label_map)) %>%
  mutate(contrast_label = unname(label_map[contrast]),
         contrast_label = factor(contrast_label, levels = contrast_order))
boot <- read_required("12wk_bootstrap_confidence_intervals.csv") %>%
  filter(metric == "PSPS_core effect", contrast %in% names(label_map)) %>%
  mutate(contrast_label = unname(label_map[contrast]),
         contrast_label = factor(contrast_label, levels = contrast_order))
scores <- read.csv(file.path(root, "results/12wk_primary/scores/12wk_primary_score_contrasts.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(score == "PSPS_core", contrast %in% names(label_map)) %>%
  select(contrast, estimate = logFC) %>%
  mutate(contrast_label = unname(label_map[contrast]),
         contrast_label = factor(contrast_label, levels = contrast_order))
loo <- read_required("12wk_psps_leave_one_gene_out_contrasts.csv") %>%
  filter(method == "leave_one_gene_out", contrast %in% names(label_map)) %>%
  mutate(
    contrast_label = unname(label_map[contrast]),
    loo_label = ifelse(score == "PSPS_core", "Full six-gene PRII",
                       paste0("without ", toupper(sub("^PSPS_without_", "", score)))),
    is_full = score == "PSPS_core",
    contrast_label = factor(contrast_label, levels = contrast_order),
    loo_label = factor(loo_label, levels = c("Full six-gene PRII", "without TXNIP", "without NLRP3",
                                              "without CASP1", "without GSDMD", "without IL1B", "without IL18"))
  )
agreement <- read.csv(file.path(root, "results/figure3_supplementary_tables/12wk_scoring_implementation_agreement.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(method_label = factor(sub("^PRII_", "", score_alt),
                               levels = c("ssGSEA", "UCell", "AddModuleScore", "AUCell")),
         status = case_when(score_alt == "AUCell" ~ "computed_not_comparable",
                            is.finite(rho) ~ "computed",
                            TRUE ~ "not computed"))
aucell_cutoff_summary <- read.csv(file.path(out_dir, "S3_AUCell_cutoff_sensitivity_summary.csv"),
                                  stringsAsFactors = FALSE, check.names = FALSE)
coverage_full <- read.csv(file.path(root, "results/12wk_primary/enrichment/12wk_gene_set_coverage.csv"),
                          stringsAsFactors = FALSE)
overlap_full <- read.csv(file.path(root, "results/12wk_primary/enrichment/12wk_gene_set_overlap.csv"),
                         stringsAsFactors = FALSE)
coverage_noncore <- read_required("12wk_robustness_gene_set_coverage.csv")

stopifnot(nrow(perm) == 2L, nrow(random_summary) == 2L, nrow(boot) == 2L,
          nrow(scores) == 2L, nrow(random) == 10000L, nrow(loo) == 14L,
          nrow(agreement) == 4L)
perm <- left_join(perm, scores %>% select(contrast_label, observed_score = estimate), by = "contrast_label")
random_summary <- left_join(random_summary, scores %>% select(contrast_label, observed_score = estimate), by = "contrast_label")

theme_s3 <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = .25, colour = "grey92"),
    axis.text = element_text(colour = "grey25", size = 7),
    axis.title = element_text(colour = "grey15", size = 7),
    panel.border = element_rect(linewidth = .35, colour = "grey35"),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    plot.tag = element_text(size = 10, face = "bold"),
    plot.tag.position = c(.01, .99),
    plot.margin = margin(4, 5, 4, 5, unit = "mm")
  )

# A. Permutation summary only: do not reconstruct a distribution from summary
# statistics. Bars show the frozen 2.5th-97.5th percentile null interval.
perm$annotation <- sprintf("Observed %.2f; P = %.3f", perm$observed,
                           perm$empirical_p_two_sided)
gA <- ggplot(perm, aes(y = contrast_label, colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .35, colour = "grey50") +
  geom_segment(aes(x = null_q025, xend = null_q975, yend = contrast_label),
               linewidth = 1.05, lineend = "round") +
  geom_point(aes(x = null_mean), shape = 21, size = 2.2, fill = "white", stroke = .55) +
  geom_point(aes(x = observed), shape = 23, size = 2.7, fill = "black", colour = "black") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(-1.2, 1.5), breaks = seq(-1, 1.5, .5), expand = c(0, 0)) +
  labs(title = "Group-label permutation (5,000 permutations)",
       x = "PRII contrast estimate (within-cohort z-score units)", y = NULL, tag = "a") +
  theme_s3

# B. Empirical null from 5,000 six-gene sets sampled from expressed non-core
# genes. Sampling was not stratified by expression or variance.
gB <- ggplot(random, aes(effect, fill = contrast_label)) +
  geom_histogram(bins = 36, colour = "white", linewidth = .2, alpha = .86) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .35, colour = "grey45") +
  geom_vline(data = random_summary, aes(xintercept = observed_score, colour = contrast_label),
             linewidth = .65, show.legend = FALSE) +
  facet_wrap(~contrast_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, guide = "none") +
  labs(title = "Random six-gene null",
       x = "Contrast estimate for sampled six-gene score", y = "Number of sets", tag = "b") +
  theme_s3 + theme(strip.text = element_text(size = 6.8, face = "bold"))

# C. Group-stratified bootstrap percentile CIs for the fixed PRII score.
gC <- ggplot(boot, aes(y = contrast_label, colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .35, colour = "grey50") +
  geom_segment(aes(x = ci95_lower, xend = ci95_upper, yend = contrast_label), linewidth = .75) +
  geom_point(aes(x = estimate), shape = 21, size = 2.6, fill = "white", stroke = .7) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.08, .08))) +
  labs(title = "Group-stratified bootstrap (5,000 resamples)",
       x = "PRII contrast estimate and percentile 95% CI", y = NULL, tag = "c") +
  theme_s3

loo$contrast_label <- factor(loo$contrast_label, levels = contrast_order)
gD <- ggplot(loo, aes(y = loo_label, colour = contrast_label)) +
  geom_vline(data = loo %>% filter(is_full), aes(xintercept = logFC, colour = contrast_label),
             linetype = "dotted", linewidth = .38, show.legend = FALSE) +
  geom_segment(aes(x = CI95_lower, xend = CI95_upper, yend = loo_label), linewidth = .52) +
  geom_point(aes(x = logFC, shape = is_full), size = 1.65, stroke = .55) +
  facet_wrap(~contrast_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.05, .06))) +
  labs(title = "Complete six-gene PRII leave-one-gene-out sensitivity",
       x = "Estimated contrast with 95% CI", y = NULL, tag = "d") +
  theme_s3 + theme(axis.text.y = element_text(size = 5.5),
                   strip.text = element_text(size = 6.2, face = "bold"),
                   panel.spacing = grid::unit(1.5, "mm"))

agreement$method_label <- factor(as.character(agreement$method_label),
                                 levels = c("ssGSEA", "UCell", "AddModuleScore", "AUCell"))
agreement_done <- agreement %>% filter(status == "computed")
agreement_missing <- agreement %>% filter(status != "computed")
gE <- ggplot(agreement_done, aes(x = rho, y = method_label)) +
  geom_segment(aes(x = 0, xend = rho, yend = method_label), colour = "grey68", linewidth = .55) +
  geom_point(shape = 21, size = 2.35, fill = "#542788", colour = "#542788") +
  geom_text(aes(label = sprintf("%.3f", rho)), nudge_x = .045, hjust = 0, size = 2.45, colour = "grey20") +
  geom_text(data = agreement_missing,
            aes(x = .47, y = method_label,
                label = ifelse(status == "computed_not_comparable",
                               "not comparable (1/6 genes)", "not computed")),
            inherit.aes = FALSE, hjust = 0, size = 2.4, colour = "grey35") +
  scale_x_continuous(limits = c(0, 1.13), breaks = c(0, .25, .5, .75, 1), expand = c(0, 0)) +
  labs(title = "Alternative full-score implementations (n = 21)",
       x = "Spearman ρ with mean-z PRII", y = NULL, tag = "e") +
  theme_s3

coverage_row <- function(label, requested, present, overlap, y) {
  data.frame(set_label = label, requested = requested, present = present,
             coverage_text = sprintf("%d/%d (%.1f%%)", present, requested, 100 * present/requested),
             overlap_text = as.character(overlap), y = y)
}
get_full_cov <- function(name) coverage_full[coverage_full$gene_set == name, , drop = FALSE]
get_noncore_cov <- function(name) coverage_noncore[coverage_noncore$set == name, , drop = FALSE]
get_core_overlap <- function(name) {
  z <- overlap_full[(overlap_full$set1 == name & overlap_full$set2 == "pyroptosis_core6") |
                    (overlap_full$set2 == name & overlap_full$set1 == "pyroptosis_core6"), , drop = FALSE]
  if (nrow(z)) z$n_overlap[1] else 0L
}
pyro_full <- get_full_cov("pyroptosis")
inflam_full <- get_full_cov("inflammasome")
pyro_noncore <- get_noncore_cov("pyroptosis_noncore")
inflam_noncore <- get_noncore_cov("inflammasome_noncore")
cov_tbl <- bind_rows(
  coverage_row("PRII core", 6, 6, "6 (self)", 5),
  coverage_row("Pyroptosis, full", pyro_full$n_genes, pyro_full$n_in_matrix,
               get_core_overlap("pyroptosis"), 4),
  coverage_row("Pyroptosis, non-core", pyro_noncore$requested, pyro_noncore$n_present,
               pyro_noncore$core_overlap, 3),
  coverage_row("Inflammasome, full", inflam_full$n_genes, inflam_full$n_in_matrix,
               get_core_overlap("inflammasome"), 2),
  coverage_row("Inflammasome, non-core", inflam_noncore$requested, inflam_noncore$n_present,
               inflam_noncore$core_overlap, 1)
)
gF <- ggplot(cov_tbl) +
  geom_hline(yintercept = 5.5, linewidth = .25, colour = "grey60") +
  geom_text(aes(x = .10, y = y, label = set_label), hjust = 0, size = 2.15, colour = "grey20") +
  geom_text(aes(x = 1.67, y = y, label = coverage_text), hjust = .5, size = 2.0, colour = "grey20") +
  geom_text(aes(x = 2.72, y = y, label = overlap_text), hjust = .5, size = 2.15, colour = "grey20") +
  annotate("text", x = .10, y = 5.82, label = "Signature", hjust = 0, fontface = "bold", size = 2.25) +
  annotate("text", x = 1.67, y = 5.82, label = "In matrix", fontface = "bold", size = 2.2) +
  annotate("text", x = 2.72, y = 5.82, label = "Core overlap", fontface = "bold", size = 2.2) +
  scale_x_continuous(limits = c(0, 3.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.55, 6.05), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Coverage and overlap with PRII core", x = NULL, y = NULL, tag = "f") +
  theme_void(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = .5, margin = margin(b = 3)),
        plot.tag = element_text(size = 10, face = "bold"),
        plot.tag.position = c(.01, .99), plot.margin = margin(4, 5, 4, 5, unit = "mm"))

fig <- (gA | gB) / gC / gD / (gE | gF) +
  plot_layout(heights = c(1.0, .72, 1.35, .9))
pdf_out <- file.path(script_dir, "Supplementary_Figure_S3_PRII_Robustness.pdf")
tiff_out <- file.path(script_dir, "Supplementary_Figure_S3_PRII_Robustness.tiff")
ggsave(pdf_out, fig, width = 180, height = 265, units = "mm", device = cairo_pdf)
tiff(tiff_out, width = 180, height = 265, units = "mm", res = 600, compression = "lzw", type = "cairo")
print(fig)
dev.off()

# Machine-readable values directly used in the figure.
figure_stats <- bind_rows(
  perm %>% transmute(panel = "A_group_label_permutation", contrast = as.character(contrast_label),
                     estimate = observed, null_mean, ci95_lower = null_q025, ci95_upper = null_q975,
                     empirical_p_two_sided, n_resamples = n_perm),
  random_summary %>% transmute(panel = "B_random_six_gene_null", contrast = as.character(contrast_label),
                               estimate = observed, null_mean, ci95_lower = null_q025, ci95_upper = null_q975,
                               empirical_p_two_sided, n_resamples = n_perm),
  boot %>% transmute(panel = "C_group_stratified_bootstrap", contrast = as.character(contrast_label),
                     estimate, null_mean = NA_real_, ci95_lower, ci95_upper,
                     empirical_p_two_sided = NA_real_, n_resamples = n_boot)
  ) %>%
  mutate(n_group_a = unname(group_a_n[contrast]),
         n_group_b = unname(group_b_n[contrast]),
         n = n_group_a + n_group_b)
write.csv(figure_stats %>% select(panel, contrast, estimate, null_mean, ci95_lower, ci95_upper,
                                  empirical_p_two_sided, n, n_group_a, n_group_b, n_resamples),
          file.path(out_dir, "S3_figure_statistics.csv"), row.names = FALSE)
write.csv(loo %>% transmute(contrast = as.character(contrast_label), score,
                            gene_removed = ifelse(is_full, "none", sub("^without ", "", as.character(loo_label))),
                            estimate = logFC, ci95_lower = CI95_lower, ci95_upper = CI95_upper,
                            p_value = P.Value, q_BH = padj_BH,
                            n_group_a = unname(group_a_n[as.character(contrast_label)]),
                            n_group_b = unname(group_b_n[as.character(contrast_label)]),
                            n = n_group_a + n_group_b),
          file.path(out_dir, "S3_full_PRII_leave_one_gene_out.csv"), row.names = FALSE)
write.csv(agreement %>% transmute(implementation = as.character(method_label),
                                  n = ifelse(is.finite(rho), n, NA_integer_),
                                  spearman_rho = rho, p_value, status,
                                  notes),
          file.path(out_dir, "S3_full_PRII_implementation_agreement.csv"), row.names = FALSE)
write.csv(cov_tbl, file.path(out_dir, "S3_signature_coverage_overlap.csv"), row.names = FALSE)
pairwise_full <- read.csv(file.path(root, "results/12wk_primary/enrichment/12wk_gene_set_overlap.csv"),
                          stringsAsFactors = FALSE) %>%
  mutate(scope = "full")
pair_core <- pairwise_full %>%
  filter((set1 == "pyroptosis" & set2 == "inflammasome") |
           (set1 == "inflammasome" & set2 == "pyroptosis"))
core_size <- as.integer(coverage_full$n_genes[coverage_full$gene_set == "pyroptosis_core6"][1])
if (nrow(pair_core) != 1L || is.na(core_size) || core_size != 6L ||
    get_core_overlap("pyroptosis") != core_size ||
    get_core_overlap("inflammasome") != core_size || pair_core$n_overlap[1] < core_size) {
  stop("Expected the full pyroptosis/inflammasome overlap to include all six PRII core genes.")
}
pairwise_noncore <- data.frame(
  set1 = "pyroptosis_noncore", set2 = "inflammasome_noncore",
  n_set1 = pair_core$n_set1[1] - core_size,
  n_set2 = pair_core$n_set2[1] - core_size,
  n_overlap = pair_core$n_overlap[1] - core_size,
  jaccard = (pair_core$n_overlap[1] - core_size) /
    ((pair_core$n_set1[1] - core_size) + (pair_core$n_set2[1] - core_size) -
       (pair_core$n_overlap[1] - core_size)),
  scope = "non-core"
)
pairwise_overlap <- bind_rows(pairwise_full, pairwise_noncore)
write.csv(pairwise_overlap, file.path(out_dir, "S3_pairwise_gene_set_overlap.csv"), row.names = FALSE)
write.csv(random_summary_all, file.path(out_dir, "S3_random_six_gene_null_summary_full.csv"), row.names = FALSE)
table_s3 <- bind_rows(
  figure_stats %>% transmute(result_family = panel, item = contrast, contrast,
                             estimate, null_mean, ci95_lower, ci95_upper,
                             p_value = empirical_p_two_sided, q_value = NA_real_, n,
                             n_group_a, n_group_b,
                             present = NA_integer_, requested = NA_integer_, core_overlap = NA_character_,
                             resamples = n_resamples, details = "Frozen resampling summary"),
  loo %>% transmute(result_family = "complete_PRII_leave_one_gene_out", item = as.character(loo_label),
                    contrast = as.character(contrast_label), estimate = logFC, null_mean = NA_real_,
                    ci95_lower = CI95_lower, ci95_upper = CI95_upper, p_value = P.Value,
                    q_value = padj_BH,
                    n = unname(contrast_n[as.character(contrast_label)]),
                    n_group_a = unname(group_a_n[as.character(contrast_label)]),
                    n_group_b = unname(group_b_n[as.character(contrast_label)]),
                    present = NA_integer_, requested = NA_integer_,
                    core_overlap = NA_character_, resamples = NA_integer_,
                    details = ifelse(is_full, "Full six-gene PRII", "One core gene removed")),
  agreement %>% transmute(result_family = "full_PRII_implementation_agreement",
                          item = as.character(method_label), contrast = "all samples",
                          estimate = rho, null_mean = NA_real_, ci95_lower = NA_real_, ci95_upper = NA_real_,
                          p_value, q_value = NA_real_, n = ifelse(is.finite(rho), n, NA_integer_),
                          n_group_a = NA_integer_, n_group_b = NA_integer_,
                          present = NA_integer_, requested = NA_integer_,
                          core_overlap = NA_character_, resamples = NA_integer_,
                          details = case_when(
                            status == "computed_not_comparable" ~ "Computed at conventional AUCell 5% cutoff but not interpretable as full six-gene agreement; see cutoff sensitivity",
                            status == "computed" ~ "computed",
                            TRUE ~ "not computed")),
  aucell_cutoff_summary %>% transmute(
    result_family = "AUCell_cutoff_sensitivity", item = paste0("AUCell_", cutoff_label),
    contrast = "all samples", estimate = spearman_rho, null_mean = NA_real_,
    ci95_lower = NA_real_, ci95_upper = NA_real_, p_value, q_value = NA_real_,
    n = n_samples, n_group_a = NA_integer_, n_group_b = NA_integer_,
    present = n_genes_in_all_samples, requested = n_core_genes,
    core_overlap = genes_in_all_samples, resamples = NA_integer_,
    details = paste0("aucMaxRank=", aucMaxRank, "; signature genes in any sample=",
                     n_genes_in_any_sample, "; per-sample genes=", min_genes_per_sample,
                     "–", max_genes_per_sample, "; ", interpretation)),
  cov_tbl %>% transmute(result_family = "signature_coverage_overlap", item = set_label,
                        contrast = NA_character_, estimate = present/requested, null_mean = NA_real_,
                        ci95_lower = NA_real_, ci95_upper = NA_real_, p_value = NA_real_, q_value = NA_real_,
                        n = NA_integer_, n_group_a = NA_integer_, n_group_b = NA_integer_,
                        present, requested, core_overlap = overlap_text,
                        resamples = NA_integer_, details = coverage_text),
  pairwise_overlap %>% transmute(result_family = "pairwise_gene_set_overlap", item = set1,
                                 contrast = set2, estimate = n_overlap, null_mean = NA_real_,
                                 ci95_lower = NA_real_, ci95_upper = NA_real_, p_value = NA_real_,
                                 q_value = NA_real_, n = NA_integer_, present = NA_integer_,
                                 n_group_a = NA_integer_, n_group_b = NA_integer_,
                                 requested = NA_integer_, core_overlap = NA_character_,
                                 resamples = NA_integer_,
                                 details = paste0(scope, " signature overlap; Jaccard=", sprintf("%.5f", jaccard))),
  random_summary_all %>% mutate(contrast_label = unname(label_map[contrast])) %>%
    transmute(result_family = "random_six_gene_null_full_summary", item = test,
                                   contrast = ifelse(is.na(contrast_label), "all 12-week samples", contrast_label),
                                   estimate = observed, null_mean, ci95_lower = null_q025,
                                   ci95_upper = null_q975, p_value = empirical_p_two_sided,
                                   q_value = NA_real_,
                                   n = ifelse(test == "PSPS effect", unname(contrast_n[contrast_label]), 21L),
                                   n_group_a = ifelse(test == "PSPS effect", unname(group_a_n[contrast_label]), NA_integer_),
                                   n_group_b = ifelse(test == "PSPS effect", unname(group_b_n[contrast_label]), NA_integer_),
                                   present = NA_integer_, requested = NA_integer_,
                                   core_overlap = NA_character_, resamples = n_perm,
                                   details = "Random six-gene null; see source limitations in legend")
)
write.csv(table_s3, file.path(out_dir, "Table_S3_PRII_robustness_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_S3.txt"))
cat("Wrote:\n", pdf_out, "\n", tiff_out, "\n", sep = "")
