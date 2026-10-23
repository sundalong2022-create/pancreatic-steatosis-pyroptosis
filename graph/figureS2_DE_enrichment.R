#!/usr/bin/env Rscript
# Supplementary Figure S2: full 12-week differential expression and pathway
# enrichment. Inferential values come from frozen results; this script only
# reshapes and visualizes those results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(scales)
})
options(stringsAsFactors = FALSE)

root <- normalizePath(".", mustWork = TRUE)
de_dir <- file.path(root, "results/12wk_primary/de")
enr_dir <- file.path(root, "results/12wk_primary/enrichment")
fig_dir <- file.path(root, "graph")
table_dir <- file.path(root, "results/supplementary_S2")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

de_summary_file <- file.path(de_dir, "12wk_de_summary.csv")
fgsea_file <- file.path(enr_dir, "12wk_all_fgsea.csv")
stopifnot(file.exists(de_summary_file), file.exists(fgsea_file))
de_summary <- read.csv(de_summary_file, check.names = FALSE)
fgsea_all <- read.csv(fgsea_file, check.names = FALSE)

contrast_ids <- c(
  "primary__HO_Safflower_vs_Safflower",
  "primary__HO_Sunflower_vs_Sunflower",
  "secondary__HO_Safflower_vs_CD",
  "secondary__HO_Sunflower_vs_CD",
  "secondary__Safflower_vs_CD",
  "secondary__Sunflower_vs_CD",
  "secondary__Sunflower_vs_Safflower",
  "secondary__HO_Sunflower_vs_HO_Safflower",
  "secondary__HO_effect_Safflower_vs_Sunflower"
)
contrast_labels <- c(
  primary__HO_Safflower_vs_Safflower = "HO safflower vs safflower",
  primary__HO_Sunflower_vs_Sunflower = "HO sunflower vs sunflower",
  secondary__HO_Safflower_vs_CD = "HO safflower vs control",
  secondary__HO_Sunflower_vs_CD = "HO sunflower vs control",
  secondary__Safflower_vs_CD = "Safflower vs control",
  secondary__Sunflower_vs_CD = "Sunflower vs control",
  secondary__Sunflower_vs_Safflower = "Sunflower vs safflower",
  secondary__HO_Sunflower_vs_HO_Safflower = "HO sunflower vs HO safflower",
  secondary__HO_effect_Safflower_vs_Sunflower = "HO effect interaction: safflower vs sunflower"
)
stopifnot(all(contrast_ids %in% de_summary$contrast))

# Preserve all nine full transcriptome tables in one file without filtering
# or recalculating any inferential values.
de_files <- file.path(de_dir, paste0(contrast_ids, ".csv"))
stopifnot(all(file.exists(de_files)))
de_all <- bind_rows(lapply(seq_along(de_files), function(i) {
  x <- read.csv(de_files[i], check.names = FALSE)
  x$source_result_file <- basename(de_files[i])
  x
}))
write.csv(de_all,
          file.path(table_dir, "S2_full_differential_expression_all_contrasts.csv"),
          row.names = FALSE, na = "")
write.csv(de_summary,
          file.path(table_dir, "S2_differential_expression_summary.csv"),
          row.names = FALSE, na = "")

# Preserve the precomputed FGSEA results across all planned contrasts.
write.csv(fgsea_all,
          file.path(table_dir, "S2_full_fgsea_all_contrasts.csv"),
          row.names = FALSE, na = "")

# Combine existing GO-BP, KEGG, and custom-set ORA outputs and retain a
# manifest that records source files returning no terms.
ora_files <- list.files(
  enr_dir,
  pattern = "^(primary|secondary)__.*__ORA_(GO_BP|KEGG|custom)_(up|down)\\.csv$",
  full.names = TRUE
)
ora_manifest <- tibble(source_result_file = basename(ora_files)) %>%
  mutate(
    contrast = str_remove(source_result_file, "__ORA_.*$"),
    enrichment = str_match(source_result_file, "__ORA_(GO_BP|KEGG|custom)_")[, 2],
    direction = str_match(source_result_file, "_(up|down)\\.csv$")[, 2]
  ) %>%
  left_join(
    de_all %>%
      filter(significant %in% TRUE) %>%
      count(contrast, direction, name = "n_directional_significant_genes"),
    by = c("contrast", "direction")
  ) %>%
  mutate(n_directional_significant_genes = coalesce(n_directional_significant_genes, 0L))
ora_tables <- lapply(seq_along(ora_files), function(i) {
  if (file.info(ora_files[i])$size <= 3) return(NULL)
  x <- read.csv(ora_files[i], check.names = FALSE)
  if (!nrow(x)) return(NULL)
  x$contrast <- ora_manifest$contrast[i]
  x$enrichment <- ora_manifest$enrichment[i]
  x$direction <- ora_manifest$direction[i]
  x$source_result_file <- ora_manifest$source_result_file[i]
  x
})
ora_manifest$n_rows <- vapply(ora_files, function(f) {
  if (file.info(f)$size <= 3) return(0L)
  nrow(read.csv(f, check.names = FALSE))
}, integer(1))
ora_manifest$status <- case_when(
  ora_manifest$n_rows > 0 ~ "results present",
  ora_manifest$n_directional_significant_genes < 5 ~ "not_run_fewer_than_5_DEGs",
  TRUE ~ "no_terms_returned_after_source_cutoffs"
)
ora_manifest$status_note <- case_when(
  ora_manifest$status == "results present" ~ "Source CSV contains enrichment rows.",
  ora_manifest$status == "not_run_fewer_than_5_DEGs" ~ "The source script skipped ORA because fewer than five directional significant genes were available.",
  TRUE ~ "The source script attempted ORA but wrote no rows; enrichment errors were caught as NULL, so this status records an empty source result rather than proving that no term exists."
)
ora_all <- bind_rows(ora_tables)
write.csv(ora_all,
          file.path(table_dir, "S2_full_ORA_all_contrasts.csv"),
          row.names = FALSE, na = "")
write.csv(ora_manifest,
          file.path(table_dir, "S2_ORA_source_manifest.csv"),
          row.names = FALSE, na = "")

primary_id <- "primary__HO_Safflower_vs_Safflower"
primary_file <- file.path(de_dir, paste0(primary_id, ".csv"))
primary_de <- read.csv(primary_file, check.names = FALSE)
stopifnot(all(c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "significant", "direction") %in% names(primary_de)))
stopifnot(nrow(primary_de) == 17815L)

# Panel A: stored DE counts for all nine planned contrasts.
de_plot <- de_summary %>%
  mutate(contrast = factor(contrast, levels = rev(contrast_ids)),
         label = factor(unname(contrast_labels[as.character(contrast)]),
                        levels = unname(contrast_labels[rev(contrast_ids)]))) %>%
  select(contrast, label, n_up, n_down, n_significant) %>%
  pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
  mutate(direction = recode(direction, n_up = "Up", n_down = "Down"),
         direction = factor(direction, levels = c("Up", "Down")))
total_labels <- de_summary %>%
  mutate(label = factor(unname(contrast_labels[contrast]),
                        levels = unname(contrast_labels[rev(contrast_ids)])))
pA <- ggplot(de_plot, aes(x = n, y = label, fill = direction)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.15) +
  geom_text(data = total_labels,
            aes(x = n_significant + 45, y = label, label = comma(n_significant)),
            inherit.aes = FALSE, hjust = 0, size = 2.5, color = "#30343B") +
  scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC"), name = "DE direction") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(xlim = c(0, max(de_summary$n_significant) + 230), clip = "off") +
  labs(title = "Differential-expression yield across all planned contrasts",
       subtitle = "FDR < 0.05 and |log2 fold-change| ≥ 1; 17,815 genes tested per contrast; red = up, blue = down",
       x = "Number of differentially expressed genes", y = NULL) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, color = "#555555"),
        plot.margin = margin(4, 22, 4, 4))

# Panels B and C: the two matched high-oleic contrasts, plotted from the full
# stored DE tables. Colors use the precomputed significance flag and direction.
make_volcano <- function(contrast, title) {
  x <- read.csv(file.path(de_dir, paste0(contrast, ".csv")), check.names = FALSE)
  stopifnot(nrow(x) == 17815L)
  x <- x %>%
    mutate(plot_class = case_when(
      significant & direction == "up" ~ "Up",
      significant & direction == "down" ~ "Down",
      TRUE ~ "Not significant"
    ),
    neglog10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)))
  ggplot(x, aes(x = logFC, y = neglog10_fdr)) +
    geom_point(data = filter(x, plot_class == "Not significant"),
               color = "#B8BEC5", alpha = 0.38, size = 0.45) +
    geom_point(data = filter(x, plot_class != "Not significant"),
               aes(color = plot_class), alpha = 0.76, size = 0.58) +
    geom_vline(xintercept = c(-1, 1), color = "#666666", linewidth = 0.35, linetype = 2) +
    geom_hline(yintercept = -log10(0.05), color = "#666666", linewidth = 0.35, linetype = 2) +
    scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC")) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    labs(title = title, subtitle = "17,815 tested genes; red = up, blue = down",
         x = expression(log[2]~fold~change), y = expression(-log[10]~FDR)) +
    theme_bw(base_size = 8) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 8.5),
          plot.subtitle = element_text(size = 6.8, color = "#555555"),
          plot.margin = margin(4, 5, 4, 5))
}
pB <- make_volcano("primary__HO_Safflower_vs_Safflower",
                   "HO safflower vs safflower")
pC <- make_volcano("primary__HO_Sunflower_vs_Sunflower",
                   "HO sunflower vs sunflower")

# Panel D: stored fgsea normalized enrichment scores across contrasts.
stopifnot(all(c("contrast", "pathway", "NES", "padj") %in% names(fgsea_all)))
contrast_short <- c(
  primary__HO_Safflower_vs_Safflower = "HO-Saf vs Saf",
  primary__HO_Sunflower_vs_Sunflower = "HO-Sun vs Sun",
  secondary__HO_Safflower_vs_CD = "HO-Saf vs CD",
  secondary__HO_Sunflower_vs_CD = "HO-Sun vs CD",
  secondary__Safflower_vs_CD = "Saf vs CD",
  secondary__Sunflower_vs_CD = "Sun vs CD",
  secondary__Sunflower_vs_Safflower = "Sun vs Saf",
  secondary__HO_Sunflower_vs_HO_Safflower = "HO-Sun vs HO-Saf",
  secondary__HO_effect_Safflower_vs_Sunflower = "HO interaction"
)
pathway_short <- c(
  pyroptosis = "Pyroptosis",
  lipid_metabolism = "Lipid",
  inflammasome = "Inflammasome",
  pyroptosis_core6 = "PRII core",
  oxidative_stress_ROS = "Oxidative stress"
)
fgsea_all <- fgsea_all %>%
  mutate(contrast_label = unname(contrast_short[contrast]),
         contrast_label = factor(contrast_label, levels = unname(contrast_short[rev(contrast_ids)])),
         pathway_label = unname(pathway_short[pathway]),
         pathway_label = factor(pathway_label, levels = unique(unname(pathway_short[fgsea_all$pathway]))),
         sig_label = ifelse(!is.na(padj) & padj < 0.05, "*", ""),
         cell_label = sprintf("%.2f%s", NES, sig_label))
nes_lim <- max(abs(fgsea_all$NES[is.finite(fgsea_all$NES)]), na.rm = TRUE)
pD <- ggplot(fgsea_all, aes(x = pathway_label, y = contrast_label, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = cell_label), size = 2.5, color = "#202124") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-nes_lim, nes_lim), name = "NES") +
  labs(title = "Predefined gene-set GSEA across contrasts",
       subtitle = "* FDR-adjusted P < 0.05", x = NULL, y = NULL) +
  theme_bw(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.3), legend.position = "right",
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 7, color = "#555555"),
        plot.margin = margin(4, 4, 4, 4))

# Primary-contrast ORA panels use existing FDR, fold-enrichment, and gene-count
# fields. Only the five lowest stored adjusted P values per panel are displayed.
read_ora <- function(db, direction) {
  f <- file.path(enr_dir, paste0(primary_id, "__ORA_", db, "_", direction, ".csv"))
  x <- read.csv(f, check.names = FALSE)
  if (!nrow(x)) return(tibble())
  x %>%
    filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    transmute(term = Description, fold_enrichment = FoldEnrichment,
              adjusted_p = p.adjust, count = Count)
}
ora_plot <- function(data, title) {
  if (!nrow(data)) {
    return(ggplot() + theme_void(base_size = 8) +
             labs(title = title, subtitle = "No terms met the display threshold"))
  }
  data <- data %>%
    mutate(term = str_remove(term, " - Mus musculus \\(house mouse\\)$"),
           term = str_wrap(term, width = 29),
           term = factor(term, levels = rev(unique(term))))
  ggplot(data, aes(x = fold_enrichment, y = term)) +
    geom_point(aes(size = count, color = -log10(adjusted_p)), alpha = 0.9) +
    scale_color_gradient(low = "#B7C6D8", high = "#7F1D3A", name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(1.8, 4.5), name = "DE genes") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title, subtitle = "Bubble size = DE genes; color = −log10(FDR)",
         x = "Fold enrichment", y = NULL) +
    theme_bw(base_size = 7.2) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6.2, color = "#262626"),
          axis.text.x = element_text(size = 6.2),
          plot.title = element_text(face = "bold", size = 8),
          plot.subtitle = element_text(size = 5.8, color = "#555555"),
          legend.position = "none", plot.margin = margin(3, 3, 3, 3))
}
pE <- ora_plot(read_ora("GO_BP", "up"), "GO BP · up")
pF <- ora_plot(read_ora("GO_BP", "down"), "GO BP · down")
pG <- ora_plot(read_ora("KEGG", "up"), "KEGG · up")
pH <- ora_plot(read_ora("KEGG", "down"), "KEGG · down")

fig <- pA /
  (pB | pC) /
  pD /
  (pE | pF | pG | pH) +
  plot_layout(heights = c(1.15, 1.3, 1.05, 1.25)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 10),
                                plot.background = element_rect(fill = "white", color = NA)))

pdf_path <- file.path(fig_dir, "FigS2_DE_enrichment.pdf")
tiff_path <- file.path(fig_dir, "FigS2_DE_enrichment.tiff")
ggsave(pdf_path, fig, width = 280, height = 215, units = "mm",
       device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
grDevices::tiff(tiff_path, width = 280 / 25.4, height = 215 / 25.4,
                units = "in", res = 600, compression = "lzw", type = "cairo",
                bg = "white")
print(fig)
grDevices::dev.off()

manifest <- c(
  "Supplementary Figure S2 and full-result tables",
  "",
  "Figure source: graph/figureS2_DE_enrichment.R",
  "Figure outputs: graph/FigS2_DE_enrichment.pdf and graph/FigS2_DE_enrichment.tiff (landscape, 600 dpi, LZW)",
  "Analysis result source: results/12wk_primary/de and results/12wk_primary/enrichment",
  "",
  "Panel A: stored limma differential-expression summary for all nine planned 12-week contrasts.",
  "Panels B-C: all 17,815 tested genes for the two matched primary contrasts (HO safflower vs safflower; HO sunflower vs sunflower); point classes use stored significance flags and directions (FDR < 0.05 and |logFC| >= 1), and the y-axis uses precomputed BH-adjusted P values.",
  "Panel D: stored fgsea results for predefined mouse gene sets across all nine contrasts; cell labels are NES and an asterisk marks stored padj < 0.05. Contrast and pathway labels are shortened for display.",
  "Panels E-H: top five stored significant ORA terms for the primary HO safflower vs safflower contrast, split by GO BP or KEGG and up/down direction. The plot uses existing FoldEnrichment, Count, and p.adjust fields.",
  "",
  "Full results: results/supplementary_S2/S2_full_differential_expression_all_contrasts.csv, S2_full_fgsea_all_contrasts.csv, and S2_full_ORA_all_contrasts.csv. S2_ORA_source_manifest.csv records directional significant-gene counts, row counts, and status for every source ORA file: results present; not_run_fewer_than_5_DEGs; or no_terms_returned_after_source_cutoffs. The source enrichment calls catch errors as NULL, so the last status describes an empty returned file and does not prove that no eligible term exists.",
  "",
  "Differential-expression source method/cutoff: shared 12-week limma model with trend=TRUE and robust=TRUE; genes were filtered as rowSums(log2(expression+1) > 1) >= 2; significance was defined in the source analysis as BH-adjusted P < 0.05 and |logFC| >= 1.",
  "Enrichment source: fgsea on the stored moderated t-statistic ranks; ORA used the stored directional significant-gene lists with mapped expression genes as universe. The source script skips ORA when a direction contains fewer than five significant genes. GO-BP, cached KEGG, and predefined custom mouse gene sets follow scripts/06_12wk_integrated_analysis.R."
)
writeLines(manifest, file.path(table_dir, "S2_data_manifest.txt"))

message("Wrote FigS2_DE_enrichment.pdf/.tiff and full S2 result tables.")
