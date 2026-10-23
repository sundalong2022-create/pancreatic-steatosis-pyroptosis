#!/usr/bin/env Rscript
# Figure S7: exploratory PRII-lipid-program analysis.
#
# This figure is descriptive. The prespecified lipid programs are shown because
# their associations with PRII are negative, discordant across contrasts, or
# sensitive to scoring/composition adjustment. Nothing here defines PRII as a
# lipid or fatty-pancreas index and no causal claim is made.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "graph/FigS7_PRII_lipid_program_exploratory.R"
root <- normalizePath(file.path(dirname(normalizePath(script_path)), ".."), mustWork = TRUE)
gd <- file.path(root, "graph")
ld <- file.path(root, "results/12wk_primary/lipid_external")
sd <- file.path(root, "results/12wk_primary/subscores")
cd <- file.path(root, "results/composition_lipid_integration")
out <- file.path(root, "results/supplementary/S7_exploratory_lipid_context")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read0 <- function(path) read_csv(path, show_col_types = FALSE, progress = FALSE)
cor0 <- read0(file.path(ld, "12wk_psps_external_lipid_correlations.csv"))
contrast0 <- read0(file.path(ld, "12wk_external_lipid_contrasts.csv"))
all0 <- read0(file.path(ld, "12wk_external_lipid_all_tests_fdr.csv"))
pooled <- read0(file.path(sd, "12wk_psps_subscore_external_lipid_pooled.csv"))
partial <- read0(file.path(sd, "12wk_psps_subscore_external_lipid_partial_group.csv"))
within <- read0(file.path(sd, "12wk_psps_subscore_external_lipid_within_group.csv"))
deconv <- read0(file.path(cd, "mouse_bulk_psps_lipid_partial_deconvolution.csv"))
weighted <- read0(file.path(cd, "sc_weighted_pseudobulk_partial_cell_proportions.csv"))
prov <- read0(file.path(ld, "external_lipid_gene_set_provenance.csv"))
coverage <- read0(file.path(cd, "gene_set_coverage_mouse_bulk.csv"))
overlap <- read0(file.path(ld, "external_lipid_gene_set_overlap.csv"))
pairwise <- read0(file.path(ld, "12wk_psps_external_lipid_pairwise_spearman.csv"))

stopifnot(nrow(cor0) == 4L, nrow(contrast0) == 8L, nrow(all0) == 12L,
          nrow(pooled) == 16L, nrow(partial) == 16L, nrow(within) == 80L,
          nrow(deconv) == 24L, nrow(weighted) == 16L)

ink <- "#243746"
teal <- "#2C7C7A"
coral <- "#CE6A58"
blue <- "#4B7C9A"
muted <- "#6F7F89"
gridc <- "#E5EBEE"
score_order <- c("Hallmark_FA_Metabolism", "Hallmark_Adipogenesis",
                 "KEGG_PPAR_Signaling", "Reactome_FA_Metabolism")
score_labels <- c(
  Hallmark_FA_Metabolism = "Hallmark FA metabolism",
  Hallmark_Adipogenesis = "Hallmark adipogenesis",
  KEGG_PPAR_Signaling = "KEGG PPAR signalling",
  Reactome_FA_Metabolism = "Reactome FA metabolism"
)
contrast_labels <- c(
  primary__HO_Safflower_vs_Safflower = "HO safflower vs safflower",
  primary__HO_Sunflower_vs_Sunflower = "HO sunflower vs sunflower"
)
fmtq <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < .001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
}
theme_s7 <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", colour = ink, size = 10),
      plot.subtitle = element_text(colour = muted, size = 7, margin = margin(b = 4)),
      axis.text = element_text(colour = ink, size = 7),
      axis.text.x = element_text(colour = muted, size = 7),
      axis.title = element_text(colour = muted, size = 7.3),
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(colour = gridc, linewidth = .3),
      plot.tag = element_text(face = "bold", colour = ink, size = 12),
      plot.tag.position = c(0, 1),
      plot.margin = margin(8, 8, 5, 8)
    )
}

# A: pooled fixed six-gene PRII-lipid correlations ---------------------------
a <- cor0 %>%
  mutate(
    score_label = factor(unname(score_labels[as.character(score)]),
                         levels = rev(unname(score_labels[score_order]))),
    direction = ifelse(rho < 0, "inverse", "positive"),
    stat = paste0("rho=", sprintf("%+.3f", rho), "\nq=", fmtq(q_global_external))
  )
write_csv(a %>% transmute(score = as.character(score), score_label = as.character(score_label),
                          rho, p_value, n, q_within_family, q_global_external),
          file.path(out, "S7A_pooled_PRII_lipid_correlations.csv"))
pA <- ggplot(a, aes(y = score_label, x = rho, colour = direction)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .4, colour = "#AAB6BC") +
  geom_segment(aes(x = 0, xend = rho, yend = score_label), linewidth = .7) +
  geom_point(size = 2.7, shape = 21, fill = "white", stroke = .7) +
  geom_text(aes(x = .05, label = stat), hjust = 0, size = 2.6, colour = ink) +
  scale_colour_manual(values = c(inverse = coral, positive = teal), guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1), expand = c(0, 0)) +
  labs(tag = "A", title = "Pooled PRII-lipid-program associations",
       subtitle = "Spearman correlations across 21 discovery samples; global BH family across 12 tests",
       x = "Spearman rho", y = NULL) + theme_s7() +
  theme(plot.margin = margin(8, 8, 5, 8))

# B: prespecified dietary contrasts (lipid-score differences only) -----------
b <- contrast0 %>%
  mutate(
    score_label = factor(unname(score_labels[as.character(score)]),
                         levels = rev(unname(score_labels[score_order]))),
    contrast_label = factor(
      paste0(unname(contrast_labels[as.character(contrast)]),
             ifelse(grepl("Safflower", contrast), " (n=4 vs 5)", " (n=3 vs 4)")),
      levels = rev(c("HO safflower vs safflower (n=4 vs 5)",
                     "HO sunflower vs sunflower (n=3 vs 4)"))
    ),
    direction = ifelse(estimate < 0, "negative", "positive"),
    stat = paste0("Delta=", sprintf("%+.3f", estimate), "\nq=", fmtq(q_global_external))
  )
write_csv(b %>% transmute(contrast = as.character(contrast_label), score = as.character(score),
                          score_label = as.character(score_label), estimate, CI95_lower, CI95_upper,
                          p_value, q_within_family, q_global_external, n_treated, n_control),
          file.path(out, "S7B_dietary_lipid_contrasts.csv"))
pB <- ggplot(b, aes(y = score_label, x = estimate, colour = direction)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .4, colour = "#AAB6BC") +
  geom_segment(aes(x = CI95_lower, xend = CI95_upper, yend = score_label), linewidth = .65) +
  geom_point(size = 2.4, shape = 21, fill = "white", stroke = .7) +
  geom_text(aes(x = .07, label = stat), hjust = 0, size = 2.3, colour = ink) +
  facet_wrap(~contrast_label, nrow = 1) +
  scale_colour_manual(values = c(negative = coral, positive = teal), guide = "none") +
  scale_x_continuous(limits = c(-.38, .55), breaks = c(-.3, -.15, 0, .15, .3, .45), expand = c(0, 0)) +
  labs(tag = "B", title = "Dietary contrasts in lipid-program scores",
       subtitle = "High-oleic minus comparator; 95% CI; global BH family across 12 tests",
       x = "ssGSEA score difference", y = NULL) + theme_s7(7.8) +
  theme(strip.background = element_rect(fill = "#F1F4F5", colour = NA),
        strip.text = element_text(face = "bold", colour = ink, size = 7),
        plot.margin = margin(8, 8, 5, 8))

# C: pooled versus group-adjusted subscore relationships ---------------------
pooled_c <- pooled %>% transmute(method, subscore, lipid_score, rho,
                                  p_value, n, q = q_external_lipid_family,
                                  analysis = "Pooled")
partial_c <- partial %>% transmute(method, subscore, lipid_score, rho = partial_rho,
                                   p_value, n, q = q_group_adjusted_family,
                                   analysis = "Group-adjusted")
c <- bind_rows(pooled_c, partial_c) %>%
  mutate(
    score_label = factor(unname(score_labels[as.character(lipid_score)]),
                         levels = rev(unname(score_labels[score_order]))),
    subscore = recode(subscore, PSPS_Pyro = "PRII-Pyro", PSPS_Inf = "PRII-Inf"),
    method = recode(method, `mean-z` = "mean-z", ssGSEA = "ssGSEA"),
    analysis = factor(analysis, levels = c("Pooled", "Group-adjusted"))
  )
write_csv(c %>% transmute(analysis = as.character(analysis), method, subscore,
                          lipid_score, score_label = as.character(score_label), rho, p_value, n, q),
          file.path(out, "S7C_subscore_pooled_group_adjusted.csv"))
pC <- ggplot(c, aes(y = score_label, x = rho, colour = analysis, shape = analysis)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .35, colour = "#AAB6BC") +
  geom_segment(aes(x = 0, xend = rho, yend = score_label), linewidth = .35, alpha = .45) +
  geom_point(size = 2.0, position = position_dodge(width = .34)) +
  facet_grid(method ~ subscore) +
  scale_colour_manual(values = c(Pooled = blue, `Group-adjusted` = coral)) +
  scale_shape_manual(values = c(Pooled = 16, `Group-adjusted` = 17)) +
  scale_x_continuous(limits = c(-1, .25), breaks = c(-1, -.75, -.5, -.25, 0), expand = c(0, 0)) +
  labs(tag = "C", title = "Subscore relationships after group adjustment",
       subtitle = "Separate 16-test BH families; most inverse component associations persist after group adjustment",
       x = "Spearman rho", y = NULL, colour = NULL, shape = NULL) + theme_s7(7.6) +
  theme(strip.background = element_rect(fill = "#F1F4F5", colour = NA),
        strip.text = element_text(face = "bold", colour = ink, size = 6.8),
        legend.position = "bottom", legend.text = element_text(size = 6.2))

# D: alternative scoring/deconvolution sensitivity ---------------------------
# These rows do not use the primary mean-z PRII. They use UCell or
# AddModuleScore on the six genes and cross-species xCell/EPIC proxies as
# covariates; the weighted pseudobulk uses 16 single-cell samples.
deconv_d <- deconv %>%
  filter(panel %in% c("xCell_core", "EPIC_core"), score_method %in% c("UCell", "AddModuleScore")) %>%
  mutate(
    source = recode(panel, xCell_core = "xCell covariate", EPIC_core = "EPIC covariate"),
    row_label = paste(source, score_method, sep = " | "),
    q_value = q_BH, n_source = n,
    lipid_score = recode(lipid_score, PPAR_Pathway = "KEGG_PPAR_Signaling")
  ) %>%
  select(score_method, lipid_score, rho, p_value, q_value, n_source, row_label)
weighted_d <- weighted %>%
  filter(panel == "major5_cell_proportions_plus_site") %>%
  mutate(
    row_label = paste("sc weighted pseudobulk", score_method, "+ proportions + site"),
    q_value = q_BH, n_source = n,
    lipid_score = recode(lipid_score, PPAR_Pathway = "KEGG_PPAR_Signaling")
  ) %>%
  select(score_method, lipid_score, rho, p_value, q_value, n_source, row_label)
d <- bind_rows(deconv_d, weighted_d) %>%
  mutate(
    row_label = factor(row_label, levels = rev(c(
      "xCell covariate | UCell", "xCell covariate | AddModuleScore",
      "EPIC covariate | UCell", "EPIC covariate | AddModuleScore",
      "sc weighted pseudobulk UCell + proportions + site",
      "sc weighted pseudobulk AddModuleScore + proportions + site"
    ))),
    score_label = factor(unname(score_labels[as.character(lipid_score)]),
                         levels = unname(score_labels[score_order])),
    tile_label = paste0(sprintf("%+.2f", rho), "\nq=", fmtq(q_value))
  )
write_csv(d %>% transmute(row_label = as.character(row_label), score_method, lipid_score,
                          score_label = as.character(score_label), rho, p_value, q_BH = q_value, n = n_source),
          file.path(out, "S7D_composition_adjusted_sensitivity.csv"))
pD <- ggplot(d, aes(x = score_label, y = row_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = tile_label), size = 1.85, lineheight = .86, colour = ink) +
  scale_fill_gradient2(low = coral, mid = "white", high = teal, midpoint = 0, limits = c(-1, 1)) +
  labs(tag = "D", title = "Alternative PRII scoring and composition sensitivity",
       subtitle = "UCell/AddModuleScore alternatives; xCell/EPIC proxies and scRNA pseudobulk (n=21/16)",
       x = NULL, y = NULL, fill = "rho") +
  theme_s7(7.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.0),
        axis.text.y = element_text(size = 5.9), panel.grid = element_blank(),
        legend.position = "right")

# Deliver the figure and exact source tables ---------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.0, 1.22)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))
pdf_path <- file.path(gd, "FigS7_PRII_lipid_program_exploratory.pdf")
tiff_path <- file.path(gd, "FigS7_PRII_lipid_program_exploratory.tiff")
png_path <- file.path(gd, "FigS7_PRII_lipid_program_exploratory_preview.png")
ggsave(pdf_path, fig, width = 12.2, height = 10.8, units = "in", device = cairo_pdf)
tiff(tiff_path, width = 12.2, height = 10.8, units = "in", res = 600,
     compression = "lzw", type = "cairo")
print(fig)
dev.off()
ggsave(png_path, fig, width = 12.2, height = 10.8, units = "in", dpi = 150)

# Preserve the exact analysis-level provenance alongside panel tables.
write_csv(all0, file.path(out, "S7_all_tests_global_FDR.csv"))
write_csv(prov, file.path(out, "S7_gene_set_provenance.csv"))
write_csv(coverage, file.path(out, "S7_gene_set_coverage_mouse_bulk.csv"))
write_csv(overlap, file.path(out, "S7_gene_set_overlap.csv"))
write_csv(pairwise, file.path(out, "S7_pairwise_program_correlations.csv"))
write_csv(deconv, file.path(out, "S7_full_bulk_deconvolution_sensitivity.csv"))
write_csv(weighted, file.path(out, "S7_full_sc_weighted_pseudobulk_sensitivity.csv"))
write_csv(within, file.path(out, "S7_full_within_group_subscore_results.csv"))

writeLines(c(
  "Figure S7 run manifest",
  paste0("generated_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R=", R.version.string),
  "input_cohort=GSE292815_KC_12wk; n=21",
  "lipid_scores=GSVA 2.x ssGSEA; alpha=0.25; normalize=TRUE; MSigDB=2026.1.Hs",
  "pooled_global_fdr=BH across 12 tests (4 pooled correlations + 8 dietary contrasts)",
  "subscore_pooled_fdr=BH across 16 tests; subscore_group_adjusted_fdr=BH across 16 tests",
  "within_group_fdr=BH across 80 diagnostic tests; group n=3-5; retained as CSV, not plotted",
  "composition_sensitivity=bulk xCell/EPIC covariates and GSE154778 weighted pseudobulk",
  "interpretation=exploratory metabolic-context analysis; no lipid mechanism or causality claim",
  "historical_source_label=PSPS/PSPS_core; manuscript_label=PRII"
), file.path(out, "S7_run_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "S7_sessionInfo.txt"))

message("Wrote ", pdf_path)
message("Wrote ", tiff_path)
message("Wrote ", png_path)
