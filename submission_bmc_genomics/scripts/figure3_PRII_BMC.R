# Figure 3 (final v2): PRII composition, sub-axes, and scoring robustness
#
# PRII is the frozen six-gene index formerly labelled PSPS_core in the
# analysis files.  This figure deliberately does not repeat the six-gene
# effect-size panel used in Figure 2.  Historical result filenames retain
# PSPS for traceability; all labels generated here use PRII terminology.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
})

set.seed(20261008)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
out_dir <- file.path(root, "submission_bmc_genomics", "figures")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

sub_dir <- file.path(root, "results", "12wk_primary", "subscores")
dat <- read_required(file.path(sub_dir, "12wk_psps_subscores_sample_scores.csv"))
sub_con <- read_required(file.path(sub_dir, "12wk_psps_subscore_contrasts.csv"))
loo_con <- read_required(file.path(sub_dir, "12wk_psps_subscore_LOO_contrasts.csv"))

clean_group <- function(x) gsub("_", "-", sub("_12wk$", "", as.character(x)))
grp_levels <- c("CD", "Safflower", "HO-Safflower", "Sunflower", "HO-Sunflower")
pal <- c("CD" = "#BDBDBD", "Safflower" = "#D6604D", "HO-Safflower" = "#B2182B",
         "Sunflower" = "#92C5DE", "HO-Sunflower" = "#2166AC")
dat$group_clean <- factor(clean_group(dat$group), levels = grp_levels)
dat$PRII <- dat$PSPS_core
dat$PRII_Pyro <- dat$PSPS_Pyro_mean_z
dat$PRII_Inf <- dat$PSPS_Inf_mean_z

theme_prii <- theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"),
    axis.text = element_text(colour = "grey25", size = 6.2),
    axis.title = element_text(colour = "grey15", size = 6.2),
    axis.ticks = element_line(linewidth = .3, colour = "grey40"),
    panel.border = element_rect(linewidth = .35, colour = "grey30"),
    strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = .2),
    strip.text = element_text(size = 6.3, colour = "grey15"),
    legend.position = "none",
    plot.tag = element_text(size = 9, face = "bold"),
    plot.tag.position = c(.01, .99),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

# -------------------------------------------------------------------------
# A. Frozen PRII composition and the two predefined sub-axes.
# -------------------------------------------------------------------------
gene_order <- c("Txnip", "Nlrp3", "Casp1", "Gsdmd", "Il1b", "Il18")
membership <- expand.grid(
  gene = gene_order,
  score = c("PRII", "PRII-Pyro", "PRII-Inf"),
  stringsAsFactors = FALSE
) %>%
  mutate(included = case_when(
    score == "PRII" ~ TRUE,
    score == "PRII-Pyro" & gene %in% c("Nlrp3", "Casp1", "Gsdmd") ~ TRUE,
    score == "PRII-Inf" & gene %in% c("Il1b", "Il18") ~ TRUE,
    TRUE ~ FALSE
  ),
  gene = factor(gene, levels = rev(gene_order)),
  score = factor(score, levels = c("PRII", "PRII-Pyro", "PRII-Inf")))

g3a <- ggplot(membership, aes(score, gene, fill = included)) +
  geom_tile(colour = "white", linewidth = .5, width = .95, height = .88) +
  geom_text(aes(label = ifelse(included, "●", "")), colour = "white", size = 3.1, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = "#542788", `FALSE` = "#E5E5E5"), guide = "none") +
  labs(x = NULL, y = NULL, tag = "a") +
  theme_prii +
  theme(axis.text.x = element_text(size = 6, face = "bold"),
        axis.text.y = element_text(face = "italic", size = 6.4),
        panel.grid = element_blank())

# -------------------------------------------------------------------------
# B. Sample-level distributions in the five dietary groups.
# -------------------------------------------------------------------------
dist_long <- bind_rows(
  data.frame(sample = dat$sample, group = dat$group_clean, axis = "PRII", score = dat$PRII),
  data.frame(sample = dat$sample, group = dat$group_clean, axis = "PRII-Pyro", score = dat$PRII_Pyro),
  data.frame(sample = dat$sample, group = dat$group_clean, axis = "PRII-Inf", score = dat$PRII_Inf)
)
dist_long$axis <- factor(dist_long$axis, levels = c("PRII", "PRII-Pyro", "PRII-Inf"))
g3b <- ggplot(dist_long, aes(group, score, colour = group)) +
  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .35, alpha = .16) +
  geom_point(position = position_jitter(width = .10, seed = 11), size = 1.05, alpha = .86) +
  facet_wrap(~axis, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal) +
  scale_x_discrete(labels = c("CD\n(n=5)", "Safflower\n(n=5)", "HO-Safflower\n(n=4)",
                              "Sunflower\n(n=4)", "HO-Sunflower\n(n=3)")) +
  labs(x = NULL, y = "Within-cohort score", tag = "b") +
  theme_prii +
  theme(axis.text.x = element_text(size = 5.0, angle = 45, hjust = 1),
        panel.spacing = unit(2, "mm"))

# -------------------------------------------------------------------------
# C. Sub-axis effects for the two prespecified primary contrasts.
# -------------------------------------------------------------------------
sc <- sub_con %>%
  filter(contrast %in% c("primary__HO_Safflower_vs_Safflower", "primary__HO_Sunflower_vs_Sunflower")) %>%
  mutate(
    axis = ifelse(grepl("Pyro", score), "PRII-Pyro", "PRII-Inf"),
    method_clean = ifelse(grepl("ssGSEA", score), "ssGSEA", "mean-z"),
    contrast_clean = ifelse(grepl("Safflower_vs_Safflower", contrast),
                            "HO-Safflower vs Safflower", "HO-Sunflower vs Sunflower")
  )
sc$axis <- factor(sc$axis, levels = c("PRII-Pyro", "PRII-Inf"))
sc$method_clean <- factor(sc$method_clean, levels = c("mean-z", "ssGSEA"))
sc$contrast_clean <- factor(sc$contrast_clean, levels = c("HO-Safflower vs Safflower", "HO-Sunflower vs Sunflower"))

g3c <- ggplot(sc, aes(estimate, axis, colour = method_clean, shape = method_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3, colour = "grey55") +
  geom_errorbar(aes(xmin = CI95_lower, xmax = CI95_upper),
                height = .16, linewidth = .42, position = position_dodge(width = .34)) +
  geom_point(size = 1.75, position = position_dodge(width = .34)) +
  facet_wrap(~contrast_clean, nrow = 1) +
  scale_colour_manual(values = c("mean-z" = "#542788", "ssGSEA" = "#7F7F7F")) +
  labs(x = "Estimated score difference", y = NULL, tag = "c") +
  theme_prii +
  theme(strip.text = element_text(size = 5.8), legend.position = "bottom",
        legend.title = element_blank(), legend.text = element_text(size = 5.2))

# -------------------------------------------------------------------------
# D. Leave-one-gene-out sensitivity for the HO-Safflower contrast.
# -------------------------------------------------------------------------
ld <- loo_con %>%
  filter(contrast == "primary__HO_Safflower_vs_Safflower", method == "LOO_mean_z") %>%
  mutate(
    axis = ifelse(grepl("Pyro", score), "PRII-Pyro", "PRII-Inf"),
    removed = sub(".*_without_", "", score)
  )
ld$axis <- factor(ld$axis, levels = c("PRII-Pyro", "PRII-Inf"))
ld$removed <- factor(ld$removed, levels = c("Nlrp3", "Casp1", "Gsdmd", "Il1b", "Il18"))
baseline <- sc %>% filter(contrast_clean == "HO-Safflower vs Safflower", method_clean == "mean-z")

g3d <- ggplot(ld, aes(estimate, removed, colour = axis)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .3, colour = "grey55") +
  geom_vline(data = baseline, aes(xintercept = estimate, colour = axis),
             linetype = "dotted", linewidth = .45, show.legend = FALSE) +
  geom_errorbar(aes(xmin = CI95_lower, xmax = CI95_upper), height = .16, linewidth = .42) +
  geom_point(size = 1.75) +
  facet_wrap(~axis, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("PRII-Pyro" = "#E66101", "PRII-Inf" = "#2C7BB6")) +
  labs(x = "Estimated difference after removing one gene", y = "Gene removed", tag = "d") +
  theme_prii +
  theme(axis.text.y = element_text(face = "italic", size = 6.2), panel.spacing = unit(2, "mm"))

# -------------------------------------------------------------------------
# E. Relationship between the two predefined biological sub-axes.
# -------------------------------------------------------------------------
rho_sub <- suppressWarnings(cor.test(dat$PRII_Pyro, dat$PRII_Inf, method = "spearman", exact = FALSE))
rho_label <- paste0("Spearman ρ = ", sprintf("%.2f", unname(rho_sub$estimate)),
                    "; n = ", nrow(dat), "\nP = ", format.pval(rho_sub$p.value, digits = 2, eps = 1e-3))
g3e <- ggplot(dat, aes(PRII_Pyro, PRII_Inf, colour = group_clean)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .4, colour = "grey45") +
  geom_point(size = 1.45, alpha = .9) +
  scale_colour_manual(values = pal) +
  annotate("text", x = -Inf, y = Inf, hjust = -.08, vjust = 1.25,
           label = rho_label, size = 1.85, lineheight = .9) +
  labs(x = "PRII-Pyro mean-z", y = "PRII-Inf mean-z", tag = "e") +
  theme_prii +
  theme(axis.title = element_text(size = 5.8))

# -------------------------------------------------------------------------
# F. Agreement between mean-z and ssGSEA implementations (descriptive).
# -------------------------------------------------------------------------
method_long <- bind_rows(
  data.frame(axis = "PRII-Pyro", mean_z = dat$PRII_Pyro, ssGSEA = dat$PSPS_Pyro_ssGSEA,
             group = dat$group_clean),
  data.frame(axis = "PRII-Inf", mean_z = dat$PRII_Inf, ssGSEA = dat$PSPS_Inf_ssGSEA,
             group = dat$group_clean)
)
method_stats <- method_long %>% group_by(axis) %>% summarise(
  rho = suppressWarnings(cor(mean_z, ssGSEA, method = "spearman", use = "complete.obs")),
  p = suppressWarnings(cor.test(mean_z, ssGSEA, method = "spearman", exact = FALSE)$p.value),
  .groups = "drop"
)
method_long <- left_join(method_long, method_stats, by = "axis")
method_long$axis <- factor(method_long$axis, levels = c("PRII-Pyro", "PRII-Inf"))
g3f <- ggplot(method_long, aes(mean_z, ssGSEA, colour = group)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .35, colour = "grey45") +
  geom_point(size = 1.35, alpha = .9) +
  facet_wrap(~axis, nrow = 1, scales = "free") +
  scale_colour_manual(values = pal) +
  geom_text(data = method_stats,
            aes(x = -Inf, y = Inf, label = paste0("Spearman ρ = ", sprintf("%.2f", rho))),
            inherit.aes = FALSE, hjust = -.08, vjust = 1.3, size = 1.8) +
  labs(x = "Mean-z subscore", y = "ssGSEA subscore", tag = "f") +
  theme_prii +
  theme(strip.text = element_text(size = 6), panel.spacing = unit(2, "mm"))

# -------------------------------------------------------------------------
# Export.  A 600 dpi LZW TIFF and vector PDF are supplied for submission.
# -------------------------------------------------------------------------
fig3 <- (g3a | g3b) / (g3c | g3d) / (g3e | g3f) +
  plot_layout(widths = c(1, 1.7), heights = c(1.02, 1.07, 1.08))

save_tiff <- function(p, f, w, h) {
  tiff(f, width = w, height = h, units = "mm", res = 600, compression = "lzw", type = "cairo")
  print(p)
  dev.off()
}

pdf_out <- file.path(out_dir, "Fig3_PRII_BMC_Genomics.pdf")
tiff_out <- file.path(out_dir, "Fig3_PRII_BMC_Genomics.tiff")
ggsave(pdf_out, fig3, width = 180, height = 178, units = "mm", device = cairo_pdf)
save_tiff(fig3, tiff_out, 180, 178)

# Machine-readable panel data and statistics used in the legend/methods.
write.csv(membership, file.path(out_dir, "Fig3_PRII_BMC_Genomics_membership.csv"), row.names = FALSE)
write.csv(sc, file.path(out_dir, "Fig3_PRII_BMC_Genomics_subscore_contrasts.csv"), row.names = FALSE)
write.csv(ld, file.path(out_dir, "Fig3_PRII_BMC_Genomics_LOO_subscore.csv"), row.names = FALSE)
write.csv(method_stats, file.path(out_dir, "Fig3_PRII_BMC_Genomics_method_agreement.csv"), row.names = FALSE)
write.csv(data.frame(metric = "PRII_Pyro_vs_PRII_Inf_Spearman", n = nrow(dat),
                     rho = unname(rho_sub$estimate), p_value = rho_sub$p.value),
          file.path(out_dir, "Fig3_PRII_BMC_Genomics_subscore_correlation.csv"), row.names = FALSE)

cat("Wrote:\n", pdf_out, "\n", tiff_out, "\n", sep = "")
