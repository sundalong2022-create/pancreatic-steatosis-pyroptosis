#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(fgsea); library(clusterProfiler); library(ggplot2); library(dplyr); library(tibble)})
root <- normalizePath('.', mustWork=TRUE); outdir <- file.path(root,'results/enrichment'); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
gs_obj <- readRDS(file.path(root,'scripts/results/gene_sets_mouse/gene_sets.rds')); pathways <- gs_obj$gene_sets
term2gene <- do.call(rbind, lapply(names(pathways), function(n) data.frame(term=n, gene=pathways[[n]], stringsAsFactors=FALSE)))
de_files <- list.files(file.path(root,'results/de'), pattern='^[^_].*\\.csv$', full.names=TRUE)
de_files <- de_files[!grepl('__significant\\.csv$',de_files) & basename(de_files) != 'differential_summary.csv']
run_one <- function(f) {
  tt <- read.csv(f,check.names=FALSE); tt <- tt[!is.na(tt$gene_symbol) & !duplicated(tt$gene_symbol),]
  stats <- tt$t; names(stats) <- tt$gene_symbol; stats <- stats[is.finite(stats)]; stats <- sort(stats,decreasing=TRUE)
  fg <- fgsea(pathways=pathways, stats=stats, minSize=10, maxSize=2000) %>% as.data.frame() %>% arrange(padj)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse=';', FUN.VALUE=character(1))
  stem <- sub('\\.csv$','',basename(f)); write.csv(fg,file.path(outdir,paste0(stem,'__fgsea.csv')),row.names=FALSE)
  sig <- tt %>% filter(!is.na(adj.P.Val),adj.P.Val<0.05,abs(logFC)>=1)
  for(direction in c('up','down')) { genes <- sig$gene_symbol[sig$direction==direction]; if(length(genes)<5) next; enr <- enricher(genes, TERM2GENE=term2gene, universe=tt$gene_symbol, minGSSize=10, maxGSSize=2000, pAdjustMethod='BH'); if(!is.null(enr)) write.csv(as.data.frame(enr),file.path(outdir,paste0(stem,'__ORA_',direction,'.csv')),row.names=FALSE) }
}
invisible(lapply(de_files,run_one)); writeLines(capture.output(sessionInfo()),file.path(root,'results/sessionInfo_enrichment.txt')); message('Completed custom gene-set GSEA and ORA.')
