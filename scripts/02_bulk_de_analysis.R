#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(limma); library(ggplot2); library(pheatmap); library(dplyr); library(tibble); library(stringr)})
root <- normalizePath('.', mustWork=TRUE)
qc_dir <- file.path(root,'results/qc'); de_dir <- file.path(root,'results/de')
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE); dir.create(de_dir, recursive=TRUE, showWarnings=FALSE)
map <- read.delim(file.path(root,'data/metadata/matrix_sample_map.tsv'), check.names=FALSE, stringsAsFactors=FALSE)
map$file <- basename(map$file)
read_matrix <- function(path) {
  x <- read.csv(gzfile(path), check.names=FALSE, stringsAsFactors=FALSE)
  sample_cols <- setdiff(names(x), c('Gene_symbol','ENTREZ_ID','Gene_name'))
  mat <- as.matrix(x[,sample_cols,drop=FALSE]); storage.mode(mat) <- 'numeric'
  genes <- trimws(x$Gene_symbol); keep <- !is.na(genes) & genes != ''
  mat <- mat[keep,,drop=FALSE]; genes <- genes[keep]
  sums <- rowsum(mat, genes, reorder=FALSE); counts <- table(genes)[rownames(sums)]
  mat <- sums / as.numeric(counts); rownames(mat) <- rownames(sums); mat
}
save_qc <- function(mat, meta, label) {
  logmat <- log2(mat+1); lib <- colSums(mat)
  write.csv(data.frame(sample=names(lib), group=meta$group, total_normalized_expression=as.numeric(lib)), file.path(qc_dir,paste0(label,'_library_size.csv')), row.names=FALSE)
  g <- ggplot(data.frame(sample=names(lib),total=as.numeric(lib),group=meta$group),aes(sample,total,fill=group))+geom_col()+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position='bottom')
  ggsave(file.path(qc_dir,paste0(label,'_library_size.pdf')),g,width=10,height=5)
  rv <- apply(logmat,1,var); top <- order(rv,decreasing=TRUE)[seq_len(min(2000,nrow(logmat)))]
  d <- as.matrix(dist(t(logmat[top,,drop=FALSE]))); rownames(d)<-colnames(d)<-colnames(logmat)
  pheatmap(d,filename=file.path(qc_dir,paste0(label,'_sample_distance.pdf')),width=8,height=7)
  pc <- prcomp(t(logmat[top,,drop=FALSE]),scale.=TRUE); pca <- data.frame(sample=rownames(pc$x),PC1=pc$x[,1],PC2=pc$x[,2],group=meta$group)
  write.csv(pca,file.path(qc_dir,paste0(label,'_PCA.csv')),row.names=FALSE)
  gp <- ggplot(pca,aes(PC1,PC2,color=group,label=sample))+geom_point(size=3)+theme_bw()+theme(legend.position='bottom')
  ggsave(file.path(qc_dir,paste0(label,'_PCA.pdf')),gp,width=8,height=6)
}
run_contrast <- function(mat,meta,label,contrast_name,contrast_formula) {
  keep <- !is.na(meta$group) & meta$group != ''; mat <- mat[,keep,drop=FALSE]; meta <- meta[keep,,drop=FALSE]
  if(length(unique(meta$group))<2 || any(table(meta$group)<2)) return(NULL)
  logmat <- log2(mat+1); logmat <- logmat[rowSums(logmat>1)>=2,,drop=FALSE]
  design <- model.matrix(~0+group,data=meta); colnames(design) <- sub('^group','',colnames(design))
  fit <- eBayes(contrasts.fit(lmFit(logmat,design), makeContrasts(contrasts=contrast_formula,levels=design)), trend=TRUE, robust=TRUE)
  tt <- topTable(fit,number=Inf,sort.by='P') %>% rownames_to_column('gene_symbol')
  tt$direction <- ifelse(tt$logFC>0,'up','down'); tt$significant <- !is.na(tt$adj.P.Val) & tt$adj.P.Val<0.05 & abs(tt$logFC)>=1
  tt$dataset <- label; tt$contrast <- contrast_name
  write.csv(tt,file.path(de_dir,paste0(label,'__',contrast_name,'.csv')),row.names=FALSE)
  write.csv(filter(tt,significant),file.path(de_dir,paste0(label,'__',contrast_name,'__significant.csv')),row.names=FALSE)
  p <- ggplot(tt,aes(AveExpr,logFC,color=significant))+geom_point(alpha=.45,size=.7)+theme_bw(); ggsave(file.path(de_dir,paste0(label,'__',contrast_name,'__MA.pdf')),p,width=7,height=5)
  p <- ggplot(tt,aes(logFC,-log10(P.Value),color=significant))+geom_point(alpha=.45,size=.7)+theme_bw(); ggsave(file.path(de_dir,paste0(label,'__',contrast_name,'__volcano.pdf')),p,width=7,height=5)
}
for(f in unique(map$file)) {
  path <- file.path(root,'data/raw/GSE292815',f); label <- str_replace(f,'_NormalizedCounts.csv.gz','')
  mat <- read_matrix(path); mm <- filter(map,file==f); mm <- mm[match(colnames(mat),mm$matrix_column),,drop=FALSE]; stopifnot(!anyNA(mm$group))
  saveRDS(list(expr=mat,metadata=mm),file.path(qc_dir,paste0(label,'_matrix.rds'))); save_qc(mat,mm,label)
  if(grepl('KC_12wk',f,fixed=TRUE)) {
    # Primary diet effects relative to CD.
    run_contrast(mat,mm,label,'Safflower12wk_vs_CD12wk','Safflower_12wk-CD_12wk')
    run_contrast(mat,mm,label,'Sunflower12wk_vs_CD12wk','Sunflower_12wk-CD_12wk')
    run_contrast(mat,mm,label,'HO_Safflower12wk_vs_CD12wk','HO_Safflower_12wk-CD_12wk')
    run_contrast(mat,mm,label,'HO_Sunflower12wk_vs_CD12wk','HO_Sunflower_12wk-CD_12wk')
    # Matched high-oleic effects within each oil background.
    run_contrast(mat,mm,label,'HO_Safflower12wk_vs_Safflower12wk','HO_Safflower_12wk-Safflower_12wk')
    run_contrast(mat,mm,label,'HO_Sunflower12wk_vs_Sunflower12wk','HO_Sunflower_12wk-Sunflower_12wk')
    # Oil-type contrasts among ordinary and high-oleic diets.
    run_contrast(mat,mm,label,'Sunflower12wk_vs_Safflower12wk','Sunflower_12wk-Safflower_12wk')
    run_contrast(mat,mm,label,'HO_Sunflower12wk_vs_HO_Safflower12wk','HO_Sunflower_12wk-HO_Safflower_12wk')
    # Difference-in-differences: does the high-oleic effect depend on oil type?
    run_contrast(mat,mm,label,'HO_effect_Safflower_vs_Sunflower','(HO_Safflower_12wk-Safflower_12wk)-(HO_Sunflower_12wk-Sunflower_12wk)')
  } else if(grepl('KC_3wk',f,fixed=TRUE)) run_contrast(mat,mm,label,'HO_Safflower3wk_vs_Safflower3wk','HO_Safflower_3wk-Safflower_3wk')
}
summary_files <- list.files(de_dir, pattern='__significant[.]csv$', full.names=TRUE)
summary <- do.call(rbind, lapply(summary_files, function(z) {
  x <- read.csv(z, check.names=FALSE)
  data.frame(file=basename(z), n_significant=nrow(x), n_up=sum(x$direction == 'up'), n_down=sum(x$direction == 'down'))
}))
write.csv(summary, file.path(de_dir, 'differential_summary.csv'), row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'results/sessionInfo_bulk_de.txt'))
message('Completed QC and limma differential analysis for GSE292815.')
