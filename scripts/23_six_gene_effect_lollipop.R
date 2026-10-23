#!/usr/bin/env Rscript
# Figure 2f: six-gene effect-size lollipop for the 12-week primary contrast.
# Inputs and normalization are frozen to scripts/05_psps_score.R.

suppressPackageStartupMessages(library(grid))
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
out_dir <- file.path(root, 'results/12wk_primary/figures')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr) == 21L, nrow(meta) == ncol(expr), all(colnames(expr) == meta$matrix_column))
genes <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
stopifnot(all(genes %in% rownames(expr)))
logexpr <- log2(expr + 1)
z <- t(scale(t(logexpr[genes,,drop=FALSE]))); z[!is.finite(z)] <- 0
treated <- which(meta$group == 'HO_Safflower_12wk'); control <- which(meta$group == 'Safflower_12wk')
stopifnot(length(treated) == 4L, length(control) == 5L)
delta_z <- rowMeans(z[,treated,drop=FALSE]) - rowMeans(z[,control,drop=FALSE])
log2fc <- rowMeans(logexpr[genes,treated,drop=FALSE]) - rowMeans(logexpr[genes,control,drop=FALSE])
# Advance RNG through the same pre-bootstrap path in scripts/14_psps_robustness.R.
n_boot <- 5000L; set.seed(292815)
universe <- rownames(logexpr)[apply(logexpr,1,sd) > 0]; universe <- setdiff(universe,genes)
for (i in seq_len(n_boot)) sample(universe,length(genes),replace=FALSE)
treated_sun <- which(meta$group == 'HO_Sunflower_12wk'); control_sun <- which(meta$group == 'Sunflower_12wk')
for (i in seq_len(n_boot)) { lab <- sample(c(rep(TRUE,length(treated)),rep(FALSE,length(control)))); idx <- c(treated,control); mean(colMeans(z[,idx[lab],drop=FALSE])) - mean(colMeans(z[,idx[!lab],drop=FALSE])) }
for (i in seq_len(n_boot)) { lab <- sample(c(rep(TRUE,length(treated_sun)),rep(FALSE,length(control_sun)))); idx <- c(treated_sun,control_sun); mean(colMeans(z[,idx[lab],drop=FALSE])) - mean(colMeans(z[,idx[!lab],drop=FALSE])) }
# Shared resampling indices for all six genes.
treated_boot <- matrix(NA_integer_,n_boot,length(treated)); control_boot <- matrix(NA_integer_,n_boot,length(control)); boot_delta_z <- matrix(NA_real_,n_boot,length(genes)); boot_log2fc <- boot_delta_z
for (b in seq_len(n_boot)) { it <- sample(treated,length(treated),TRUE); ic <- sample(control,length(control),TRUE); treated_boot[b,] <- it; control_boot[b,] <- ic; boot_delta_z[b,] <- rowMeans(z[,it,drop=FALSE]) - rowMeans(z[,ic,drop=FALSE]); boot_log2fc[b,] <- rowMeans(logexpr[genes,it,drop=FALSE]) - rowMeans(logexpr[genes,ic,drop=FALSE]) }
ci_z <- apply(boot_delta_z,2,quantile,probs=c(.025,.975),names=FALSE); ci_lfc <- apply(boot_log2fc,2,quantile,probs=c(.025,.975),names=FALSE)
delta_psps <- mean(delta_z); boot_delta_psps <- rowMeans(boot_delta_z); ci_psps <- quantile(boot_delta_psps,probs=c(.025,.975),names=FALSE)
psps_scores <- read.csv(file.path(root,'results/12wk_primary/scores/12wk_sample_scores.csv')); observed_psps <- with(psps_scores, mean(PSPS_core[group=='HO_Safflower_12wk']) - mean(PSPS_core[group=='Safflower_12wk']))
stopifnot(abs(delta_psps-observed_psps) < 1e-10)
fig3a <- read.csv(file.path(root,'results/12wk_primary/robustness/12wk_bootstrap_confidence_intervals.csv')); fig3a <- fig3a[fig3a$metric=='PSPS_core effect' & fig3a$contrast=='primary__HO_Safflower_vs_Safflower',,drop=FALSE]
if (nrow(fig3a)==1L) stopifnot(max(abs(ci_psps-c(fig3a$ci95_lower,fig3a$ci95_upper))) < 1e-10)
ord <- order(delta_z,decreasing=TRUE)
effect_table <- data.frame(gene=genes[ord],delta_z=unname(delta_z[ord]),ci95_lower=ci_z[1,ord],ci95_upper=ci_z[2,ord],log2FC=unname(log2fc[ord]),log2FC_ci95_lower=ci_lfc[1,ord],log2FC_ci95_upper=ci_lfc[2,ord],stringsAsFactors=FALSE)
effect_table$direction <- ifelse(effect_table$delta_z >= 0,'positive','negative'); effect_table$rank_delta_z <- seq_len(nrow(effect_table))
write.csv(effect_table,file.path(out_dir,'fig2f_six_gene_effects.csv'),row.names=FALSE)
write.csv(data.frame(contrast='HO_Safflower_12wk - Safflower_12wk',n_treated=length(treated),n_control=length(control),n_boot=n_boot,seed=292815,mean_delta_z=delta_psps,observed_psps_delta=observed_psps,additive_difference=delta_psps-observed_psps,psps_ci95_lower=ci_psps[1],psps_ci95_upper=ci_psps[2],fig3a_ci95_lower=if(nrow(fig3a)==1L) fig3a$ci95_lower else NA_real_,fig3a_ci95_upper=if(nrow(fig3a)==1L) fig3a$ci95_upper else NA_real_),file.path(out_dir,'fig2f_additivity_check.csv'),row.names=FALSE)
saveRDS(list(treated=treated_boot,control=control_boot,seed=292815,contrast='HO_Safflower_12wk - Safflower_12wk'),file.path(out_dir,'fig2f_bootstrap_indices.rds'))

draw_figure <- function(fontfamily = "sans") {
  grid.newpage()
  left <- .19; right <- .75; bottom <- .22; top <- .83
  x_lim <- c(-3, 3)
  x_map <- function(x) left + (x - x_lim[1]) / diff(x_lim) * (right - left)
  y_map <- function(i) top - (i - 1) / (nrow(effect_table) - 1) * (top - bottom)
  tick <- seq(-3, 3, 1)
  y <- vapply(seq_len(nrow(effect_table)), y_map, numeric(1))
  x0 <- x_map(0)
  colors <- ifelse(effect_table$direction == 'positive', '#B2182B', '#2166AC')
  clip <- function(x) pmax(x_lim[1], pmin(x_lim[2], x))
  gp_text <- function(size, face = 'plain', col = '#222222') gpar(fontfamily = fontfamily, fontsize = size, fontface = face, col = col)

  grid.text('2f', unit(.03, 'npc'), unit(.965, 'npc'), just = c('left', 'top'), gp = gp_text(8, 'bold'))
  grid.text('HO-Safflower vs Safflower', unit(left, 'npc'), unit(.965, 'npc'), just = c('left', 'top'), gp = gp_text(7.2))
  grid.text('log2FC', unit(.91, 'npc'), unit(top + .045, 'npc'), gp = gp_text(6.7, 'bold'))
  grid.lines(unit(c(left, right), 'npc'), unit(c(bottom - .025, bottom - .025), 'npc'), gp = gpar(col = '#222222', lwd = .5))
  grid.lines(unit(c(x0, x0), 'npc'), unit(c(bottom - .025, top + .005), 'npc'), gp = gpar(col = '#8C8C8C', lwd = .45))
  grid.segments(x0 = unit(x_map(delta_psps), 'npc'), x1 = unit(x_map(delta_psps), 'npc'), y0 = unit(bottom - .025, 'npc'), y1 = unit(top + .005, 'npc'), gp = gpar(col = '#4D4D4D', lwd = .55, lty = 3))
  grid.text(bquote(Delta*PSPS == .(sprintf('%.3f', delta_psps))), unit(x_map(delta_psps), 'npc'), unit(top + .022, 'npc'), gp = gp_text(5.7, col = '#4D4D4D'))

  for (i in seq_len(nrow(effect_table))) {
    grid.segments(x0 = unit(x0, 'npc'), x1 = unit(x_map(effect_table$delta_z[i]), 'npc'), y0 = unit(y[i], 'npc'), y1 = unit(y[i], 'npc'), gp = gpar(col = colors[i], lwd = 1.2))
    grid.segments(x0 = unit(x_map(clip(effect_table$ci95_lower[i])), 'npc'), x1 = unit(x_map(clip(effect_table$ci95_upper[i])), 'npc'), y0 = unit(y[i], 'npc'), y1 = unit(y[i], 'npc'), gp = gpar(col = colors[i], lwd = .8))
    cap <- .018
    for (xx in c(clip(effect_table$ci95_lower[i]), clip(effect_table$ci95_upper[i]))) {
      grid.segments(x0 = unit(x_map(xx), 'npc'), x1 = unit(x_map(xx), 'npc'), y0 = unit(y[i] - cap, 'npc'), y1 = unit(y[i] + cap, 'npc'), gp = gpar(col = colors[i], lwd = .8))
    }
    grid.points(unit(x_map(effect_table$delta_z[i]), 'npc'), unit(y[i], 'npc'), pch = 16, size = unit(1.65, 'mm'), gp = gpar(col = colors[i]))
    grid.text(effect_table$gene[i], unit(left - .025, 'npc'), unit(y[i], 'npc'), just = 'right', gp = gp_text(7.2, 'italic'))
    grid.text(sprintf('%+.2f', effect_table$log2FC[i]), unit(.91, 'npc'), unit(y[i], 'npc'), gp = gp_text(6.7))
  }
  tick_labels <- ifelse(tick == 0, '0', sprintf('%+.0f', tick))
  grid.text(tick_labels, unit(x_map(tick), 'npc'), unit(bottom - .055, 'npc'), gp = gp_text(6.2))
  grid.text(bquote(Delta*z~'(standardized effect size)'), unit((left + right) / 2, 'npc'), unit(.055, 'npc'), gp = gp_text(6.8))
  grid.text('dashed line: observed PSPS difference = (1/6) * sum(Delta-z); descriptive only', unit((left + right) / 2, 'npc'), unit(.115, 'npc'), gp = gp_text(5.2, col = '#555555'))
}
pdf(file.path(out_dir,'fig2f_six_gene_effect_lollipop.pdf'),width=85/25.4,height=58/25.4,useDingbats=FALSE); draw_figure('sans'); dev.off()
postscript(file.path(out_dir,'fig2f_six_gene_effect_lollipop.eps'),width=85/25.4,height=58/25.4,onefile=FALSE,paper='special',horizontal=FALSE); draw_figure('Helvetica'); dev.off()
ragg::agg_tiff(file.path(out_dir,'fig2f_six_gene_effect_lollipop.tiff'),width=85/25.4,height=58/25.4,units='in',res=600,compression='lzw',background='white'); draw_figure('sans'); dev.off()
message(sprintf('Wrote Figure 2f; mean Delta-z = %.6f; PSPS CI = [%.6f, %.6f]',delta_psps,ci_psps[1],ci_psps[2]))
