#!/usr/bin/env Rscript
# Figure S4: sample-aware single-cell context and cell-type localization of PRII.
# Original result columns retain the historical PSPS name; all new labels and exports use PRII.
# Cell-level tests are not used for inference.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "graph/FigS4_PRII_scRNA_context.R"
root <- normalizePath(file.path(dirname(normalizePath(script_path)), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results", "supplementary", "S4_single_cell_context")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sc_dir <- file.path(root, "results", "sc_location_GSE154778")
up_dir <- file.path(sc_dir, "upgrade")
per_cell <- read.csv(file.path(up_dir, "per_cell_ucell_scores.csv"), check.names = FALSE)
pb <- read.csv(file.path(up_dir, "sample_celltype_pseudobulk_psps.csv"), check.names = FALSE)
cell_summary <- read.csv(file.path(up_dir, "cell_type_ucell_summary.csv"), check.names = FALSE)
effects <- read.csv(file.path(up_dir, "pseudobulk_celltype_effects.csv"), check.names = FALSE)
detection <- read.csv(file.path(sc_dir, "cell_type_gene_detection.csv"), check.names = FALSE)
high_conf <- read.csv(file.path(up_dir, "cell_type_ucell_summary_margin_ge_0.1.csv"), check.names = FALSE)
epi_sensitivity <- read.csv(file.path(up_dir, "cell_type_ucell_summary_epithelial_compartment.csv"), check.names = FALSE)
paired_direction <- read.csv(file.path(up_dir, "pseudobulk_paired_direction.csv"), check.names = FALSE)

core <- c("TXNIP", "NLRP3", "CASP1", "GSDMD", "IL1B", "IL18")
cell_order <- c("macrophage", "epithelial_tumor", "acinar", "stellate", "ductal", "t_cell", "endothelial")
comparison_order <- c("epithelial_tumor", "ductal", "stellate", "t_cell", "acinar", "endothelial")
cell_labels <- c(
  macrophage = "Macrophage", epithelial_tumor = "Tumor epithelium",
  acinar = "Acinar", stellate = "Stellate", ductal = "Ductal",
  t_cell = "T cell", endothelial = "Endothelial", beta = "Beta", alpha = "Alpha"
)
class_colors <- c(primary = "#3B7C88", metastatic = "#D38A5D")
ink <- "#263746"
muted <- "#657782"
grid_col <- "#E6ECEF"

stopifnot(nrow(per_cell) == 14926L)
stopifnot(length(unique(per_cell$sample_id)) == 16L)
stopifnot(all(core %in% colnames(detection)))
stopifnot(all(c("pseudobulk_UCell", "pseudobulk_UCell_primary_only") %in% effects$model))
stopifnot(all(cell_order %in% cell_summary$cell_type))

rename_prii <- function(x) {
  names(x) <- sub("PSPS_UCell", "PRII_UCell", names(x), fixed = TRUE)
  names(x) <- sub("PSPS_core", "PRII", names(x), fixed = TRUE)
  x
}

# Panel A: show sample and cell denominators separately as 100% stacked bars.
sample_meta <- unique(per_cell[c("sample_id", "sample_class")])
class_order <- c("primary", "metastatic")
sample_tab <- table(factor(sample_meta$sample_class, levels = class_order))
cell_tab <- table(factor(per_cell$sample_class, levels = class_order))
composition <- rbind(
  data.frame(metric = "Samples (n = 16)", sample_class = class_order, n = as.integer(sample_tab)),
  data.frame(metric = "Cells (n = 14,926)", sample_class = class_order, n = as.integer(cell_tab))
)
composition$metric <- factor(composition$metric, levels = rev(c("Samples (n = 16)", "Cells (n = 14,926)")))
composition$sample_class <- factor(composition$sample_class, levels = class_order)
totals <- aggregate(n ~ metric, composition, sum)
composition <- merge(composition, totals, by = "metric", suffixes = c("", "_total"))
composition$share <- composition$n / composition$n_total
composition$label <- ifelse(composition$metric == "Samples (n = 16)",
                            sprintf("%d (%.0f%%)", composition$n, 100 * composition$share),
                            sprintf("%s (%.0f%%)", format(composition$n, big.mark = ","), 100 * composition$share))
write.csv(composition, file.path(out_dir, "S4A_cohort_composition.csv"), row.names = FALSE)

p_a <- ggplot(composition, aes(x = share, y = metric, fill = sample_class)) +
  geom_col(position = "stack", width = 0.58, colour = "white", linewidth = 0.45) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 2.8) +
  scale_fill_manual(values = class_colors, breaks = class_order,
                    labels = c("Primary", "Metastatic")) +
  scale_x_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(tag = "A", title = "Cohort composition", x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 10),
        axis.text = element_text(colour = muted, size = 7.2),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 7),
        plot.tag = element_text(face = "bold", colour = ink, size = 12),
        plot.tag.position = c(0, 1), plot.margin = margin(7, 7, 4, 7))

# Panel B: fraction of cells with detectable transcript; this is not pathway activation.
det_types <- c("macrophage", "epithelial_tumor", "ductal", "stellate", "acinar", "t_cell", "endothelial")
det <- detection[detection$cell_type %in% det_types, c("cell_type", core), drop = FALSE]
det_long <- do.call(rbind, lapply(core, function(g) {
  data.frame(cell_type = det$cell_type, gene = g, detection_fraction = det[[g]])
}))
det_long <- merge(det_long, cell_summary[c("cell_type", "n_cells", "n_samples")],
                  by = "cell_type", all.x = TRUE, sort = FALSE)
det_long$cell_type <- factor(det_long$cell_type, levels = rev(det_types))
det_long$gene <- factor(det_long$gene, levels = core)
det_long$detection_percent <- 100 * det_long$detection_fraction
det_long$label <- sprintf("%.0f", det_long$detection_percent)
write.csv(det_long, file.path(out_dir, "S4B_PRII_core_gene_detection.csv"), row.names = FALSE)

p_b <- ggplot(det_long, aes(x = gene, y = cell_type, fill = detection_percent)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.25, colour = "#263746") +
  scale_fill_gradient(low = "#F3F6F7", high = "#2E7585", limits = c(0, 100),
                      breaks = c(0, 25, 50, 75, 100), name = "Detected cells (%)") +
  scale_y_discrete(labels = cell_labels) +
  labs(tag = "B", title = "PRII-gene detection across cell types", x = NULL, y = NULL) +
  theme_minimal(base_size = 8.2) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 10),
        axis.text.x = element_text(colour = ink, face = "italic", size = 7, angle = 35, hjust = 1),
        axis.text.y = element_text(colour = ink, size = 7.2),
        axis.title = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 6.8),
        legend.text = element_text(size = 6.5), legend.key.width = unit(10, "mm"),
        plot.tag = element_text(face = "bold", colour = ink, size = 12),
        plot.tag.position = c(0, 1), plot.margin = margin(7, 7, 4, 7))

# Panel C: cell-level score distribution is descriptive; open circles are sample-level units.
per_cell <- per_cell[per_cell$cell_type %in% cell_order, , drop = FALSE]
pb_major <- pb[pb$cell_type %in% cell_order, , drop = FALSE]
cell_summary <- cell_summary[match(cell_order, cell_summary$cell_type), , drop = FALSE]
tick_labels <- c(macrophage = "Macrophage", epithelial_tumor = "Tumor\nepithelium",
                 acinar = "Acinar", stellate = "Stellate", ductal = "Ductal",
                 t_cell = "T cell", endothelial = "Endothelial")
per_cell$cell_type <- factor(per_cell$cell_type, levels = cell_order)
pb_major$cell_type <- factor(pb_major$cell_type, levels = cell_order)
pb_major <- rename_prii(pb_major)
per_cell <- rename_prii(per_cell)
cell_summary <- rename_prii(cell_summary)
write.csv(cell_summary, file.path(out_dir, "S4C_cell_type_PRII_UCell_summary.csv"), row.names = FALSE)
write.csv(pb_major, file.path(out_dir, "S4C_sample_celltype_PRII_UCell.csv"), row.names = FALSE)

p_c <- ggplot(per_cell, aes(x = cell_type, y = PRII_UCell, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.55, width = 0.88) +
  geom_point(data = pb_major, aes(x = cell_type, y = PRII_UCell), inherit.aes = FALSE,
             position = position_jitter(width = 0.10, height = 0, seed = 247),
             shape = 21, size = 1.55, fill = "white", colour = ink, stroke = 0.4) +
  scale_fill_manual(values = c(macrophage = "#D88757", epithelial_tumor = "#7EA9B1",
                               acinar = "#A7BAC0", stellate = "#9BAFB5", ductal = "#88A6AE",
                               t_cell = "#88A6AE", endothelial = "#A7BAC0"), guide = "none") +
  scale_x_discrete(labels = tick_labels) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "C", title = "PRII UCell scores by annotated cell type",
       subtitle = "Violin: cell-level distribution; open circles: sample × cell-type pseudobulk values",
       x = NULL, y = "PRII UCell score") +
  theme_classic(base_size = 8.8) +
  theme(plot.title = element_text(face = "bold", colour = ink, size = 10),
        plot.subtitle = element_text(colour = muted, size = 7.4, margin = margin(b = 5)),
        axis.text.x = element_text(colour = ink, size = 6.8, lineheight = 0.9),
        axis.text.y = element_text(colour = muted, size = 7.7),
        axis.title.y = element_text(colour = muted, size = 8),
        axis.line.x = element_line(colour = "#AAB6BC", linewidth = 0.35),
        plot.tag = element_text(face = "bold", colour = ink, size = 12),
        plot.tag.position = c(0, 1), plot.margin = margin(8, 9, 4, 7))

# Panel D: mixed-model contrasts; BH correction is within each six-contrast model family.
effects <- effects[effects$response == "PSPS_UCell" &
                     effects$model %in% c("pseudobulk_UCell", "pseudobulk_UCell_primary_only"), , drop = FALSE]
effects$comparison_label <- unname(cell_labels[effects$comparison])
effects$BH_q <- effects$BH_FDR
effects$effect_type <- "Macrophage minus comparator"
effects <- rename_prii(effects)
effects$source_model <- effects$model
effects$analysis_population <- ifelse(effects$source_model == "pseudobulk_UCell",
                                      "primary_and_metastatic", "primary_only")
effects$response <- "PRII_UCell"
effects$fdr_family <- paste0(effects$analysis_population, "__six_macrophage_reference_contrasts")
names(effects)[names(effects) == "BH_FDR"] <- "q_value"
write.csv(effects, file.path(out_dir, "S4D_pseudobulk_PRII_UCell_model_effects.csv"), row.names = FALSE)
write.csv(rename_prii(high_conf), file.path(out_dir, "S4E_high_confidence_label_sensitivity.csv"), row.names = FALSE)
write.csv(rename_prii(epi_sensitivity), file.path(out_dir, "S4E_epithelial_compartment_sensitivity.csv"), row.names = FALSE)
write.csv(rename_prii(paired_direction), file.path(out_dir, "S4E_paired_direction_summary.csv"), row.names = FALSE)

forest_plot <- function(model_id, title, subtitle, tag, significant_colour) {
  d <- effects[effects$model == model_id, , drop = FALSE]
  stopifnot(nrow(d) == length(comparison_order))
  d$y <- length(comparison_order) - match(d$comparison, comparison_order) + 1L
  d$significance <- ifelse(d$BH_q < 0.05, "BH q < 0.05", "BH q ≥ 0.05")
  qlab <- ifelse(d$BH_q < 0.001, "<0.001", sprintf("%.3f", d$BH_q))
  ggplot(d, aes(y = y)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#9EABB1") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.16,
                  linewidth = 0.65, colour = muted) +
    geom_point(aes(x = macrophage_minus_comparison, fill = significance),
               shape = 21, size = 2.3, colour = ink, stroke = 0.4) +
    geom_text(aes(x = 0.225, label = qlab), hjust = 0, colour = ink, size = 2.3) +
    annotate("text", x = 0.225, y = 6.65, label = "BH q", hjust = 0,
             colour = muted, fontface = "bold", size = 2.45) +
    scale_fill_manual(values = c("BH q < 0.05" = significant_colour,
                                 "BH q ≥ 0.05" = "white"), name = NULL, guide = "none") +
    scale_y_continuous(breaks = 1:6,
                       labels = unname(cell_labels[rev(comparison_order)]),
                       limits = c(0.45, 6.85), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-0.10, 0.34), breaks = c(-0.1, 0, 0.1, 0.2),
                       labels = c("−0.1", "0", "0.1", "0.2"), expand = c(0, 0)) +
    labs(tag = tag, title = title, subtitle = subtitle,
         x = "Macrophage − comparator difference", y = NULL) +
    theme_classic(base_size = 7.8) +
    theme(plot.title = element_text(face = "bold", colour = ink, size = 9.4),
          plot.subtitle = element_text(colour = muted, size = 6.8, margin = margin(b = 4)),
          axis.text.y = element_text(colour = ink, size = 6.9),
          axis.text.x = element_text(colour = muted, size = 6.7),
          axis.title.x = element_text(colour = muted, size = 7, margin = margin(t = 5)),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.3),
          legend.position = "none",
          plot.tag = element_text(face = "bold", colour = ink, size = 12),
          plot.tag.position = c(0, 1), plot.margin = margin(8, 7, 4, 7))
}

p_d1 <- forest_plot("pseudobulk_UCell", "All samples",
                    "65 sample × cell-type units; 16 samples; adjusted for sample class",
                    "D1", "#3B7C88")
p_d2 <- forest_plot("pseudobulk_UCell_primary_only", "Primary-only sensitivity",
                    "48 units; 10 primary samples; no BH q < 0.05",
                    "D2", "#D38A5D")

figure <- (p_a | p_b) / p_c / (p_d1 | p_d2) +
  plot_layout(heights = c(0.78, 1.14, 1.24), guides = "collect") +
  plot_annotation(
    title = "Single-cell context and cell-type localization of PRII-related transcript scores in GSE154778",
    subtitle = "Transcriptomic localization in PDAC; this cohort is not a dietary pancreas model",
    theme = theme(plot.title = element_text(face = "bold", colour = ink, size = 12, hjust = 0),
                  plot.subtitle = element_text(colour = muted, size = 8.7, hjust = 0),
                  plot.margin = margin(8, 10, 6, 10))
  ) & theme(legend.position = "bottom")

pdf_file <- file.path(root, "graph", "FigS4_PRII_scRNA_context.pdf")
tiff_file <- file.path(root, "graph", "FigS4_PRII_scRNA_context.tiff")
png_file <- file.path(root, "graph", "FigS4_PRII_scRNA_context_preview.png")
ggsave(pdf_file, figure, width = 7.5, height = 9.0, device = cairo_pdf, bg = "white", limitsize = FALSE)
tiff(tiff_file, width = 7.5, height = 9.0, units = "in", res = 600,
     compression = "lzw", type = "cairo", bg = "white")
print(figure)
dev.off()
ggsave(png_file, figure, width = 7.5, height = 9.0, dpi = 220, bg = "white", limitsize = FALSE)

writeLines(c(
  "Figure S4 PRII single-cell context run manifest",
  "dataset=GSE154778 processed pancreatic ductal adenocarcinoma scRNA-seq",
  "n_cells=14926; n_samples=16; primary_samples=10; metastatic_samples=6",
  "primary_cells=8000; metastatic_cells=6926",
  "cell_type_labels=existing marker-module assignments; ductal/epithelial overlap acknowledged",
  "cell_score=UCell on full gene-ranking background; six PRII core genes; maxRank=1500",
  "pseudobulk_unit=sample_id x cell_type; minimum 10 cells per unit",
  "mixed_model=score ~ cell_type + sample_class + (1|sample_id); primary-only sensitivity repeats without sample_class",
  "FDR=Benjamini-Hochberg within each six-comparison macrophage-reference model family",
  "gene_detection=per-cell transcript detection fraction; not pathway activation",
  "interpretation=transcriptomic context only; no dietary effect or functional pyroptosis inference",
  paste0("source_matrix=", file.path(root, "data/raw/GSE154778/GSE154778_dgeMtx.csv.gz")),
  paste0("created=", format(Sys.time(), tz = "UTC", usetz = TRUE))
), file.path(out_dir, "S4_run_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "S4_sessionInfo.txt"))

cat("Wrote Figure S4 files:\n", pdf_file, "\n", tiff_file, "\n", png_file,
    "\nSupplementary data: ", out_dir, "\n", sep = "")
