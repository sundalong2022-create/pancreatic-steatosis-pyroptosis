#!/usr/bin/env Rscript
# Primary 12-week integrated analysis for GSE292815.
# Uses one shared five-group limma model and prespecified primary/secondary contrasts.

suppressPackageStartupMessages({
  library(limma); library(fgsea); library(clusterProfiler); library(org.Mm.eg.db)
  library(ggplot2); library(dplyr); library(tibble); library(tidyr)
})
set.seed(292815)

root <- normalizePath('.', mustWork=TRUE)
base <- file.path(root, 'results/12wk_primary')
de_dir <- file.path(base, 'de'); enr_dir <- file.path(base, 'enrichment'); fig_dir <- file.path(base, 'figures')
dir.create(de_dir, recursive=TRUE, showWarnings=FALSE); dir.create(enr_dir, recursive=TRUE, showWarnings=FALSE); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr) == 21, all(colnames(expr) == meta$matrix_column), all(!is.na(meta$group)))
logexpr <- log2(expr + 1)
keep_gene <- rowSums(logexpr > 1) >= 2
logexpr <- logexpr[keep_gene, , drop=FALSE]

design <- model.matrix(~0 + group, data=meta)
colnames(design) <- sub('^group', '', colnames(design))
contrasts <- c(
  primary__HO_Safflower_vs_Safflower = 'HO_Safflower_12wk-Safflower_12wk',
  primary__HO_Sunflower_vs_Sunflower = 'HO_Sunflower_12wk-Sunflower_12wk',
  secondary__HO_Safflower_vs_CD = 'HO_Safflower_12wk-CD_12wk',
  secondary__HO_Sunflower_vs_CD = 'HO_Sunflower_12wk-CD_12wk',
  secondary__Safflower_vs_CD = 'Safflower_12wk-CD_12wk',
  secondary__Sunflower_vs_CD = 'Sunflower_12wk-CD_12wk',
  secondary__Sunflower_vs_Safflower = 'Sunflower_12wk-Safflower_12wk',
  secondary__HO_Sunflower_vs_HO_Safflower = 'HO_Sunflower_12wk-HO_Safflower_12wk',
  secondary__HO_effect_Safflower_vs_Sunflower = '(HO_Safflower_12wk-Safflower_12wk)-(HO_Sunflower_12wk-Sunflower_12wk)'
)

gs_obj <- readRDS(file.path(root, 'scripts/results/gene_sets_mouse/gene_sets.rds'))
pathways <- gs_obj$gene_sets
pathways$pyroptosis_core6 <- intersect(gs_obj$core_axis, rownames(logexpr))
term2gene <- do.call(rbind, lapply(names(pathways), function(n) data.frame(term=n, gene=pathways[[n]], stringsAsFactors=FALSE)))
kegg_internal <- requireNamespace('KEGG.db', quietly=TRUE)
kegg_online <- requireNamespace('KEGGREST', quietly=TRUE)
use_internal_kegg <- kegg_internal
kegg_cache <- file.path(enr_dir, 'kegg_mmu_term2gene_name.rds')
kegg_map <- if (!kegg_internal && file.exists(kegg_cache)) tryCatch(readRDS(kegg_cache), error=function(e) NULL) else NULL
kegg_available <- kegg_internal || (!is.null(kegg_map) && nrow(kegg_map$term2gene) > 0)
kegg_status <- if (kegg_internal) 'KEGG status: available via KEGG.db' else if (kegg_available) 'KEGG status: available via cached KEGGREST mapping' else if (kegg_online) 'KEGG status: KEGGREST installed; cached mapping unavailable, KEGG CSVs are empty' else 'KEGG status: unavailable; install KEGGREST or KEGG.db, KEGG CSVs are empty'
writeLines(kegg_status, file.path(enr_dir, 'KEGG_status.txt'))
nm <- names(pathways)
overlap <- do.call(rbind, lapply(seq_len(length(nm)-1), function(i) do.call(rbind, lapply((i+1):length(nm), function(j) {
  a <- pathways[[i]]; b <- pathways[[j]]; data.frame(set1=nm[i], set2=nm[j], n_set1=length(a), n_set2=length(b), n_overlap=length(intersect(a,b)), jaccard=length(intersect(a,b))/length(union(a,b)))
}))))
write.csv(overlap, file.path(enr_dir, '12wk_gene_set_overlap.csv'), row.names=FALSE)
coverage <- data.frame(gene_set=names(pathways), n_genes=vapply(pathways,length,integer(1)), n_in_matrix=vapply(pathways,function(x)length(intersect(x,rownames(logexpr))),integer(1)))
coverage$coverage <- coverage$n_in_matrix/coverage$n_genes
write.csv(coverage, file.path(enr_dir, '12wk_gene_set_coverage.csv'), row.names=FALSE)

symbol_to_entrez <- function(symbols) {
  x <- suppressMessages(bitr(unique(symbols), fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Mm.eg.db)); unique(x$ENTREZID)
}
run_kegg_ora <- function(genes, universe) {
  if (!length(genes) || !kegg_available) return(NULL)
  ids <- symbol_to_entrez(genes); bg <- symbol_to_entrez(universe)
  if (!length(ids) || !length(bg)) return(NULL)
  if (kegg_internal) {
    return(tryCatch(enrichKEGG(gene=ids, organism='mmu', universe=bg, minGSSize=3, maxGSSize=2000, pAdjustMethod='BH', use_internal_data=TRUE), error=function(e) NULL))
  }
  tryCatch(enricher(gene=ids, TERM2GENE=kegg_map$term2gene, TERM2NAME=kegg_map$term2name, universe=bg, minGSSize=3, maxGSSize=2000, pAdjustMethod='BH'), error=function(e) NULL)
}

all_summary <- list(); all_fgsea <- list()
for (nm_contrast in names(contrasts)) {
  cf <- contrasts[[nm_contrast]]
  con <- makeContrasts(contrasts=cf, levels=design)
  fit <- eBayes(contrasts.fit(lmFit(logexpr, design), con), trend=TRUE, robust=TRUE)
  tt <- topTable(fit, number=Inf, sort.by='P') %>% rownames_to_column('gene_symbol')
  tt$direction <- ifelse(tt$logFC > 0, 'up', 'down'); tt$significant <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1; tt$contrast <- nm_contrast
  write.csv(tt, file.path(de_dir, paste0(nm_contrast, '.csv')), row.names=FALSE)
  write.csv(filter(tt, significant), file.path(de_dir, paste0(nm_contrast, '__significant.csv')), row.names=FALSE)
  write.csv(filter(tt, significant, direction=='up'), file.path(de_dir, paste0(nm_contrast, '__up.csv')), row.names=FALSE)
  write.csv(filter(tt, significant, direction=='down'), file.path(de_dir, paste0(nm_contrast, '__down.csv')), row.names=FALSE)
  p <- ggplot(tt, aes(AveExpr, logFC, color=significant)) + geom_point(alpha=.45, size=.7) + theme_bw() + scale_color_manual(values=c('FALSE'='grey70', 'TRUE'='#C0392B')) + labs(title=nm_contrast, x='Mean log2 expression', y='log2 fold-change')
  ggsave(file.path(fig_dir, paste0(nm_contrast, '__MA.pdf')), p, width=7, height=5)
  p <- ggplot(tt, aes(logFC, -log10(P.Value), color=significant)) + geom_point(alpha=.45, size=.7) + theme_bw() + scale_color_manual(values=c('FALSE'='grey70', 'TRUE'='#C0392B')) + labs(title=nm_contrast, x='log2 fold-change', y='-log10(P value)')
  ggsave(file.path(fig_dir, paste0(nm_contrast, '__volcano.pdf')), p, width=7, height=5)
  all_summary[[nm_contrast]] <- data.frame(contrast=nm_contrast, n_tested=nrow(tt), n_significant=sum(tt$significant), n_up=sum(tt$significant & tt$direction=='up'), n_down=sum(tt$significant & tt$direction=='down'), min_FDR=min(tt$adj.P.Val,na.rm=TRUE))
  stats <- tt %>% filter(is.finite(t), !duplicated(gene_symbol)); s <- stats$t; names(s)<-stats$gene_symbol; s<-sort(s,decreasing=TRUE)
  fg <- fgsea(pathways=pathways, stats=s, minSize=3, maxSize=2000) %>% as.data.frame() %>% arrange(padj)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse=';', FUN.VALUE=character(1)); fg$contrast <- nm_contrast
  write.csv(fg, file.path(enr_dir, paste0(nm_contrast, '__fgsea.csv')), row.names=FALSE); all_fgsea[[nm_contrast]] <- fg
  sig <- filter(tt, significant)
  for (direction in c('up','down')) {
    genes <- sig$gene_symbol[sig$direction==direction]
    custom <- if (length(genes) >= 5) enricher(genes, TERM2GENE=term2gene, universe=rownames(logexpr), minGSSize=3, maxGSSize=2000, pAdjustMethod='BH') else NULL
    if (is.null(custom)) write.csv(data.frame(), file.path(enr_dir, paste0(nm_contrast, '__ORA_custom_', direction, '.csv')), row.names=FALSE) else write.csv(as.data.frame(custom), file.path(enr_dir, paste0(nm_contrast, '__ORA_custom_', direction, '.csv')), row.names=FALSE)
    eg <- if (length(genes) >= 5) tryCatch(enrichGO(gene=symbol_to_entrez(genes), universe=symbol_to_entrez(rownames(logexpr)), OrgDb=org.Mm.eg.db, ont='BP', pAdjustMethod='BH', readable=TRUE), error=function(e) NULL) else NULL
    if (is.null(eg)) write.csv(data.frame(), file.path(enr_dir, paste0(nm_contrast, '__ORA_GO_BP_', direction, '.csv')), row.names=FALSE) else write.csv(as.data.frame(eg), file.path(enr_dir, paste0(nm_contrast, '__ORA_GO_BP_', direction, '.csv')), row.names=FALSE)
    ek <- if (length(genes) >= 5) run_kegg_ora(genes, rownames(logexpr)) else NULL
    if (is.null(ek)) write.csv(data.frame(), file.path(enr_dir, paste0(nm_contrast, '__ORA_KEGG_', direction, '.csv')), row.names=FALSE) else write.csv(as.data.frame(ek), file.path(enr_dir, paste0(nm_contrast, '__ORA_KEGG_', direction, '.csv')), row.names=FALSE)
  }
  core <- tt %>% filter(gene_symbol %in% gs_obj$core_axis) %>% select(gene_symbol, logFC, P.Value, adj.P.Val, significant, direction)
  write.csv(core, file.path(de_dir, paste0(nm_contrast, '__core_axis.csv')), row.names=FALSE)
}
summary <- bind_rows(all_summary); write.csv(summary, file.path(de_dir, '12wk_de_summary.csv'), row.names=FALSE)
fg_all <- bind_rows(all_fgsea); write.csv(fg_all, file.path(enr_dir, '12wk_all_fgsea.csv'), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(base, 'sessionInfo_12wk_integrated.txt'))
message('Completed 12-week shared-model DE, custom GSEA/ORA, GO-BP/KEGG ORA, overlap audit, and core-axis tables.')
