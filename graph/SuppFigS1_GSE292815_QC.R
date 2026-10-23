#!/usr/bin/env Rscript
# Supplementary Figure S1: sample inventory and expression-level QC for the
# three GSE292815 pancreatic matrices available in this project. The script
# reuses normalized GEO matrices; it does not rerun differential testing.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath('.', mustWork = TRUE)
fig_dir <- file.path(root, 'graph')
out_dir <- file.path(root, 'results/supplementary/S1_qc')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  '3-week KC' = file.path(root, 'results/qc/GSE292815_KC_3wk_matrix.rds'),
  '12-week KC' = file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'),
  '8-week KCO' = file.path(root, 'results/qc/GSE292815_KCO_8wk_matrix.rds')
)
stopifnot(all(file.exists(inputs)))

palette <- c('CD' = '#606C7A', 'Safflower' = '#D38A2E',
             'HO-Safflower' = '#B9573C', 'Sunflower' = '#D1B32B',
             'HO-Sunflower' = '#318A83', 'KCO-Safflower' = '#6D70A7')
group_levels <- names(palette)
group_short <- c('CD' = 'CD', 'Safflower' = 'Saf', 'HO-Safflower' = 'HO-Saf',
                 'Sunflower' = 'Sun', 'HO-Sunflower' = 'HO-Sun',
                 'KCO-Safflower' = 'KCO-Saf')
cohort_levels <- c('3-week KC', '12-week KC', '8-week KCO')

cohort_data <- lapply(names(inputs), function(cohort_name) {
  obj <- readRDS(inputs[[cohort_name]])
  expr <- as.matrix(obj$expr)
  meta <- obj$metadata
  stopifnot(ncol(expr) == nrow(meta), identical(colnames(expr), meta$matrix_column))
  stopifnot(all(is.finite(expr)), all(expr >= 0))
  meta$sample <- meta$matrix_column
  meta$cohort <- factor(cohort_name, levels = cohort_levels)
  meta$group_display <- gsub('_', '-', sub('_([0-9]+wk)$', '', meta$group))
  meta$group_display <- factor(meta$group_display, levels = group_levels)
  stopifnot(!anyNA(meta$group_display))
  meta$sample_index <- sub('^.*_', '', meta$matrix_column)
  meta$sample_short <- paste0(group_short[as.character(meta$group_display)], meta$sample_index)
  meta$sex <- ifelse(grepl(', male,', meta$sample_title, fixed = TRUE), 'Male',
                     ifelse(grepl(', female,', meta$sample_title, fixed = TRUE), 'Female', 'Not stated'))
  log_expr <- log2(expr + 1)
  keep <- apply(log_expr, 1, sd) > 0
  pca <- prcomp(t(log_expr[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- data.frame(
    sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    group = meta$group_display[match(rownames(pca$x), meta$sample)],
    sex = meta$sex[match(rownames(pca$x), meta$sample)],
    sample_short = meta$sample_short[match(rownames(pca$x), meta$sample)],
    cohort = factor(cohort_name, levels = cohort_levels), stringsAsFactors = FALSE
  )
  corr <- cor(log_expr, method = 'pearson', use = 'pairwise.complete.obs')
  distance <- 1 - corr
  ij <- which(upper.tri(distance), arr.ind = TRUE)
  pair_df <- data.frame(
    cohort = factor(cohort_name, levels = cohort_levels),
    sample_1 = rownames(distance)[ij[, 1]], sample_2 = colnames(distance)[ij[, 2]],
    correlation_distance = distance[ij], stringsAsFactors = FALSE
  )
  sum_df <- data.frame(
    sample = meta$sample, sample_short = meta$sample_short,
    cohort = factor(cohort_name, levels = cohort_levels), group = meta$group_display, sex = meta$sex,
    total_normalized_expression = colSums(expr), stringsAsFactors = FALSE
  )
  list(expr = expr, meta = meta, log_expr = log_expr, pca = pca_df,
       var_pct = var_pct, variable_genes = sum(keep), sums = sum_df,
       distance = distance, pairs = pair_df)
})
names(cohort_data) <- names(inputs)

# All-sample metadata and sample allocation counts, preserving source annotations.
meta_all <- do.call(rbind, lapply(cohort_data, `[[`, 'meta'))
counts <- aggregate(sample ~ cohort + group_display, data = meta_all,
                    FUN = length)
names(counts)[names(counts) == 'sample'] <- 'n'
counts$group_display <- factor(counts$group_display, levels = group_levels)
counts$cohort <- factor(counts$cohort, levels = cohort_levels)
write.csv(meta_all[, c('cohort', 'sample', 'gsm', 'sample_title', 'group', 'group_display', 'sex')],
          file.path(out_dir, 'S1_sample_metadata.csv'), row.names = FALSE)
write.csv(counts, file.path(out_dir, 'S1_group_counts.csv'), row.names = FALSE)

pca_all <- do.call(rbind, lapply(cohort_data, `[[`, 'pca'))
pca_all$cohort <- factor(pca_all$cohort, levels = cohort_levels)
write.csv(pca_all, file.path(out_dir, 'S1_PCA_coordinates.csv'), row.names = FALSE)
var_tbl <- do.call(rbind, lapply(names(cohort_data), function(nm) {
  z <- cohort_data[[nm]]
  data.frame(cohort = nm, samples = nrow(z$pca), variable_genes = z$variable_genes,
             PC1_percent = z$var_pct[1], PC2_percent = z$var_pct[2])
}))
write.csv(var_tbl, file.path(out_dir, 'S1_PCA_variance.csv'), row.names = FALSE)

sums_all <- do.call(rbind, lapply(cohort_data, `[[`, 'sums'))
sums_all$cohort <- factor(sums_all$cohort, levels = cohort_levels)
write.csv(sums_all, file.path(out_dir, 'S1_total_normalized_expression.csv'), row.names = FALSE)
pairs_all <- do.call(rbind, lapply(cohort_data, `[[`, 'pairs'))
pairs_all$cohort <- factor(pairs_all$cohort, levels = cohort_levels)
write.csv(pairs_all, file.path(out_dir, 'S1_pairwise_correlation_distances.csv'), row.names = FALSE)

base_theme <- theme_classic(base_size = 9, base_family = 'sans') +
  theme(plot.title = element_text(face = 'bold', size = 9.2, hjust = 0),
        plot.subtitle = element_text(size = 7.4, colour = 'grey30'),
        axis.title = element_text(size = 7.6),
        axis.text = element_text(size = 6.8, colour = 'grey20'),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 7.5),
        legend.title = element_blank(), legend.text = element_text(size = 7.2),
        legend.position = 'bottom', plot.margin = margin(4, 5, 4, 5))

p_a <- ggplot(counts, aes(group_display, n, fill = group_display)) +
  geom_col(width = .68, colour = 'white', linewidth = .25) +
  geom_text(aes(label = paste0('n=', n)), vjust = -.45, size = 2.7) +
  facet_wrap(~cohort, nrow = 1, scales = 'free_x', strip.position = 'top') +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_x_discrete(labels = group_short, drop = TRUE) +
  scale_y_continuous(limits = c(0, max(counts$n) + 1.5), breaks = seq(0, 10, 2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = 'Sample inventory by cohort and reported diet group', x = NULL, y = 'Samples') +
  base_theme + theme(legend.position = 'none', axis.text.x = element_text(angle = 25, hjust = 1))

pca_plot <- function(nm, tag) {
  d <- pca_all[pca_all$cohort == nm, , drop = FALSE]
  v <- var_tbl[var_tbl$cohort == nm, ]
  ggplot(d, aes(PC1, PC2, colour = group, shape = sex)) +
    geom_point(size = 2.25, alpha = .9) +
    scale_colour_manual(values = palette, drop = TRUE) +
    scale_shape_manual(values = c('Male' = 16, 'Female' = 17, 'Not stated' = 4), drop = TRUE) +
    labs(title = paste0(nm, ' (n=', nrow(d), ')'),
         subtitle = paste0('PC1 ', sprintf('%.1f', v$PC1_percent), '%; PC2 ', sprintf('%.1f', v$PC2_percent), '%'),
         x = paste0('PC1 (', sprintf('%.1f', v$PC1_percent), '%)'),
         y = paste0('PC2 (', sprintf('%.1f', v$PC2_percent), '%)')) +
    base_theme + theme(legend.position = if (nm == '8-week KCO') 'bottom' else 'none') +
    guides(colour = 'none')
}
p_b <- pca_plot('3-week KC', 'B')
p_c <- pca_plot('12-week KC', 'C')
p_d <- pca_plot('8-week KCO', 'D')

p_e <- ggplot(sums_all, aes(group, total_normalized_expression / 1e7,
                            colour = group)) +
  geom_boxplot(width = .6, outlier.shape = NA, linewidth = .4, alpha = .12) +
  geom_point(size = 1.7, alpha = .85,
             position = position_jitter(width = .07, height = 0)) +
  facet_wrap(~cohort, nrow = 1, scales = 'free') +
  scale_colour_manual(values = palette, drop = TRUE) +
  scale_x_discrete(labels = group_short, drop = TRUE) +
  labs(title = 'Total normalized expression per sample',
       subtitle = 'Descriptive sums of the supplied normalized GEO matrices; not raw sequencing depth',
       x = NULL, y = expression(sum~normalized~expression~(times~10^7))) +
  base_theme + theme(legend.position = 'none', axis.text.x = element_text(angle = 28, hjust = 1))

p_f <- ggplot(pairs_all, aes(cohort, correlation_distance, colour = cohort)) +
  geom_boxplot(width = .58, outlier.shape = NA, linewidth = .45, alpha = .18) +
  geom_point(size = 1.1, alpha = .32, position = position_jitter(width = .12, height = 0)) +
  scale_colour_manual(values = c('3-week KC' = '#6D70A7', '12-week KC' = '#318A83',
                                 '8-week KCO' = '#B9573C')) +
  labs(title = 'Pairwise sample correlation distance',
       subtitle = '1 − Pearson r across log2(normalized count + 1) genes; all pairs shown descriptively',
       x = NULL, y = '1 − Pearson r') +
  base_theme + theme(legend.position = 'none')

fig <- p_a / (p_b | p_c | p_d) / (p_e | p_f) +
  plot_layout(guides = 'collect', heights = c(.75, 1.0, 1.05)) +
  plot_annotation(tag_levels = 'A',
    title = 'GSE292815: sample inventory and expression-level quality overview',
    subtitle = 'Available 3-week KC, 12-week KC, and 8-week KCO normalized GEO expression matrices (n = 8, 21, and 10)',
    caption = paste0('Only normalized GEO matrices were available; raw FASTQ-level QC was not assessed. All samples are shown; none were removed and no objective outlier-exclusion rule was applied.\n',
      'Allocation/randomization procedures were not reported in the available study description or GEO metadata; this does not establish whether allocation was randomized.'),
    theme = theme(plot.title = element_text(face = 'bold', size = 12.2),
                  plot.subtitle = element_text(size = 8.5, colour = 'grey25'),
                  plot.caption = element_text(size = 6.7, colour = 'grey35', hjust = 0),
                  legend.position = 'bottom'))

pdf_file <- file.path(fig_dir, 'SuppFigS1_GSE292815_QC.pdf')
tiff_file <- file.path(fig_dir, 'SuppFigS1_GSE292815_QC.tiff')
ggsave(pdf_file, fig, width = 12, height = 10.5, device = cairo_pdf, bg = 'white')
raw_tiff <- tempfile(pattern = 'S1_raw_', fileext = '.tiff')
ggsave(raw_tiff, fig, width = 12, height = 10.5, dpi = 600, device = 'tiff', bg = 'white')
py_file <- tempfile(pattern = 'S1_finalize_', fileext = '.py')
writeLines(c(
  'from PIL import Image', 'import sys',
  'im = Image.open(sys.argv[1])',
  "kw = {'dpi': (600, 600), 'compression': 'tiff_lzw'}",
  "if 'icc_profile' in im.info: kw['icc_profile'] = im.info['icc_profile']",
  'im.save(sys.argv[2], **kw)'
), py_file)
python <- Sys.which('python3')
if (!nzchar(python)) { unlink(c(raw_tiff, py_file)); stop('python3 is required to write a compressed 600-dpi TIFF') }
status <- system2(python, args = c(shQuote(py_file), shQuote(raw_tiff), shQuote(tiff_file)))
unlink(c(raw_tiff, py_file))
if (status != 0L) stop('Pillow TIFF encoding failed')

writeLines(c(
  'Supplementary Figure S1: sample inventory and expression-level QC.',
  'Included normalized GEO matrices: 3-week KC (n=8), 12-week KC (n=21), and 8-week KCO (n=10).',
  paste0('PCA: log2(normalized count + 1), all nonzero-variance genes, centered and unscaled. ',
         paste(apply(var_tbl, 1, function(x) paste0(x['cohort'], ': PC1=', sprintf('%.2f', as.numeric(x['PC1_percent'])),
                     '%, PC2=', sprintf('%.2f', as.numeric(x['PC2_percent'])), '%')), collapse='; ')),
  'Pairwise sample distance: 1 - Pearson correlation across log2(normalized count + 1) genes; descriptive only.',
  'Panel E sums are total normalized expression, not raw sequencing depth.',
  'Only normalized GEO matrices were available; raw FASTQ-level QC was not assessed.',
  'All available samples were plotted; no objective outlier-exclusion rule was applied and no samples were removed.',
  'Allocation/randomization procedures were not reported in the available study description or GEO metadata.',
  'No differential expression or pathway enrichment statistics are recalculated by this script.'
), file.path(out_dir, 'S1_run_manifest.txt'))
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'S1_sessionInfo.txt'))
message('Wrote: ', pdf_file, ' and ', tiff_file)
