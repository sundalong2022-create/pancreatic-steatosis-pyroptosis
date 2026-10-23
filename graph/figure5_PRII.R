#!/usr/bin/env Rscript

# Figure 5. TCGA-PAAD clinical association and tumor-purity sensitivity.
# Run from the project root:
#   Rscript graph/figure5_PRII.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

root <- normalizePath(".", mustWork = TRUE)
outdir <- file.path(root, "graph")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

audit_file <- file.path(root, "results/tcga_paad_external/tcga_sample_count_audit.csv")
cox_file <- file.path(root, "results/tcga_paad_external/tcga_stepwise_cox_forest_data.csv")
cox_summary_file <- file.path(root, "results/tcga_paad_external/tcga_stepwise_cox_model_summary.csv")
purity_file <- file.path(root, "results/immune_validation/tcga_sample_scores_with_purity.csv")

required_files <- c(audit_file, cox_file, cox_summary_file, purity_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing input file(s): ", paste(missing_files, collapse = ", "))

audit <- read_csv(audit_file, show_col_types = FALSE)
cox <- read_csv(cox_file, show_col_types = FALSE)
cox_summary <- read_csv(cox_summary_file, show_col_types = FALSE)
purity_data <- read_csv(purity_file, show_col_types = FALSE) %>%
  transmute(PRII = PSPS_core, purity = as.numeric(purity)) %>%
  filter(is.finite(PRII), is.finite(purity))

get_n <- function(metric) {
  x <- audit$n[audit$metric == metric]
  if (length(x) != 1L) stop("Expected one audit entry for metric: ", metric)
  as.integer(x)
}

# Panel A uses analysis-specific denominators; these are not a serial attrition path.
model_n <- function(model) {
  x <- cox_summary$n_complete[cox_summary$model == model]
  if (length(x) != 1L) stop("Expected one model summary for: ", model)
  as.integer(x)
}

panel_a_data <- tibble(
  x = 1:5,
  n = c(
    get_n("primary_tumor_total"),
    get_n("shared_complete"),
    model_n("Model_1_univariable"),
    model_n("Model_3_plus_stage"),
    model_n("Model_4_plus_purity")
  ),
  label = c(
    "All primary\ntumors",
    "Purity-complete\nimmune cohort",
    "OS-complete\ncohort",
    "OS + age +\nstage",
    "OS + age + stage\n+ purity"
  )
)

ink <- "#243746"
teal <- "#327C7A"
coral <- "#D16B58"
muted <- "#70818D"
card_fill <- "#EEF3F4"

p_a <- ggplot(panel_a_data) +
  geom_rect(aes(xmin = x - 0.46, xmax = x + 0.46, ymin = 0.18, ymax = 0.95),
            fill = card_fill, colour = "#D5E0E3", linewidth = 0.45) +
  geom_text(aes(x = x, y = 0.70, label = n),
            colour = ink, fontface = "bold", size = 8.2) +
  geom_text(aes(x = x, y = 0.39, label = label),
            colour = muted, size = 3.5, lineheight = 0.95) +
  annotate("text", x = 3, y = 0.04,
           label = "Each n is specific to its analysis; counts are not sequential exclusions.",
           colour = muted, size = 3.1, fontface = "italic") +
  scale_x_continuous(limits = c(0.48, 5.52), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.03), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Analysis-specific TCGA-PAAD sample counts") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(colour = ink, face = "bold", size = 12,
                                  margin = margin(b = 8)),
        plot.margin = margin(5, 8, 3, 8))

cox_prii <- cox %>%
  filter(term == "PSPS_core") %>%
  mutate(
    model_order = match(model, c("Model_1_univariable", "Model_2_plus_age",
                                 "Model_3_plus_stage", "Model_4_plus_purity")),
    y = 5 - model_order,
    model_label = factor(
      c("Univariable", "+ age", "+ age + stage", "+ age + stage + purity")[model_order],
      levels = rev(c("Univariable", "+ age", "+ age + stage", "+ age + stage + purity"))
    ),
    nominal = p_value < 0.05,
    stat_label = sprintf("%.2f (%.2f–%.2f); P = %s",
                         HR, CI_low, CI_high,
                         ifelse(p_value < 0.001, "<0.001", sprintf("%.4f", p_value))),
    n_event_label = sprintf("%d / %d", n_complete, events)
  ) %>%
  arrange(model_order)

if (nrow(cox_prii) != 4L) stop("Expected four PRII rows in the nested Cox data.")

# Forest and model statistics share one coordinate system to keep all estimates aligned.
p_b <- ggplot(cox_prii, aes(y = y)) +
  geom_vline(xintercept = 1, linewidth = 0.55, linetype = "dashed", colour = "#AAB6BC") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", height = 0.16,
                linewidth = 0.8, colour = muted) +
  geom_point(aes(x = HR, fill = nominal), shape = 21, size = 3.3,
             colour = "white", stroke = 0.55) +
  geom_text(aes(x = 3.15, label = n_event_label), hjust = 0,
            colour = ink, size = 3.15) +
  geom_text(aes(x = 5.45, label = stat_label), hjust = 0,
            colour = ink, size = 3.05) +
  annotate("text", x = 3.15, y = 4.47, label = "n / events",
           hjust = 0, colour = muted, fontface = "bold", size = 3.0) +
  annotate("text", x = 5.45, y = 4.47, label = "HR (95% CI); P",
           hjust = 0, colour = muted, fontface = "bold", size = 3.0) +
  scale_fill_manual(values = c(`FALSE` = coral, `TRUE` = teal), guide = "none") +
  scale_x_log10(
    limits = c(0.42, 16),
    breaks = c(0.5, 1, 2),
    labels = c("0.5", "1", "2"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = cox_prii$y,
    labels = cox_prii$model_label,
    limits = c(0.45, 4.62),
    expand = c(0, 0)
  ) +
  labs(title = "Continuous PRII in nested overall-survival Cox models",
       x = "Hazard ratio per one-unit increase in PRII (log scale)", y = NULL) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(colour = ink, face = "bold", size = 12,
                              margin = margin(b = 7)),
    axis.text.y = element_text(colour = ink, size = 9.5),
    axis.text.x = element_text(colour = muted, size = 9),
    axis.title.x = element_text(colour = muted, size = 9, margin = margin(t = 8)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 5, 3, 5)
  )

purity_cor <- suppressWarnings(cor.test(purity_data$PRII, purity_data$purity,
                                        method = "spearman", exact = FALSE))
rho <- unname(purity_cor$estimate)
pval <- purity_cor$p.value
cor_label <- sprintf("Spearman ρ = %.2f\nP = %s; n = %d",
                     rho,
                     ifelse(pval < 0.001, "<0.001", sprintf("%.4f", pval)),
                     nrow(purity_data))

p_c <- ggplot(purity_data, aes(x = purity, y = PRII)) +
  geom_point(shape = 21, size = 2.05, alpha = 0.73,
             fill = teal, colour = "white", stroke = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = coral, fill = "#F1D5CF", linewidth = 0.9) +
  annotate("label", x = 0.12, y = Inf, label = cor_label,
           hjust = 0, vjust = 1.12, size = 3.4, colour = ink,
           fill = "white", linewidth = 0.25, label.padding = unit(0.16, "lines")) +
  scale_x_continuous(limits = c(0.08, 0.92), breaks = seq(0.2, 0.8, 0.2),
                     labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(title = "PRII and measured tumor purity",
       x = "ABSOLUTE tumor purity", y = "PRII") +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(colour = ink, face = "bold", size = 12,
                              margin = margin(b = 7)),
    axis.text = element_text(colour = muted, size = 9),
    axis.title = element_text(colour = muted, size = 9),
    plot.margin = margin(5, 10, 3, 5)
  )

figure <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.72, 1.55), widths = c(1, 0.82)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(colour = ink, face = "bold", size = 14),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(8, 10, 6, 10)
    )
  )

pdf_file <- file.path(outdir, "Fig5_PRII.pdf")
tiff_file <- file.path(outdir, "Fig5_PRII.tiff")
ggsave(pdf_file, figure, width = 13.2, height = 8.4, device = cairo_pdf,
       bg = "white", limitsize = FALSE)
tiff(tiff_file, width = 13.2, height = 8.4, units = "in", res = 600,
     compression = "lzw", type = "cairo", bg = "white")
print(figure)
dev.off()

write_csv(tibble(
  n = nrow(purity_data),
  spearman_rho = rho,
  p_value = pval,
  estimate_method = "Spearman correlation; complete PRII and ABSOLUTE purity"
), file.path(outdir, "Fig5_PRII_purity_correlation.csv"))

cat("Wrote:\n", pdf_file, "\n", tiff_file, "\n",
    "Purity correlation: n=", nrow(purity_data), ", rho=", round(rho, 3),
    ", P=", format.pval(pval, digits=3), "\n", sep = "")
