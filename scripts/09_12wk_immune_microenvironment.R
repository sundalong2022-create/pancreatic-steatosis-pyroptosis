#!/usr/bin/env Rscript
# 12-week mouse bulk immune microenvironment association analysis.
# Immune markers are prespecified and all six PSPS_core genes are excluded
# before scoring to avoid direct score overlap.
suppressPackageStartupMessages({library(limma); library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(pheatmap)})
set.seed(292815)
root <- normalizePath('.', mustWork=TRUE)
base <- file.path(root,'results/12wk_primary'); outdir <- file.path(base,'immune'); figdir <- file.path(outdir,'figures')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(root,'results/qc/GSE292815_KC_12wk_matrix.rds')); expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr)==nrow(meta), all(colnames(expr)==meta$matrix_column))
logexpr <- log2(expr+1)
ps <- read.csv(file.path(base,'scores/12wk_sample_scores.csv'),check.names=FALSE,stringsAsFactors=FALSE)
ps <- ps[match(colnames(expr),ps$sample),,drop=FALSE]; stopifnot(all(ps$sample==colnames(expr)))
# Compact cell-lineage/state modules. Names are mouse gene symbols.
marker_sets <- list(
  leukocyte=c('Ptprc','H2-Aa','H2-Ab1','H2-Eb1','Cd74','Lsp1'),
  t_cell=c('Cd3d','Cd3e','Cd3g','Cd4','Cd8a','Il7r','Ltb'),
  cytotoxic=c('Nkg7','Gzmb','Gzmk','Prf1','Ctsw','Gzma'),
  b_cell=c('Ms4a1','Cd79a','Cd79b','Cd37','Cd74','H2-Aa','H2-Ab1','Cd19','Cd83'),
  myeloid=c('Lst1','Tyrobp','Aif1','Ctss','Ctsd','Lgals3','S100a10'),
  macrophage_TAM=c('Adgre1','Cd68','Csf1r','C1qa','C1qb','C1qc','Ms4a7','Mertk','Mafb','Cd163','Mrc1','Marco','Lgals3'),
  neutrophil=c('S100a8','S100a9','Lcn2','Retnlg','Cxcr2','Csf3r','Ngp','Camp'),
  NK=c('Ncr1','Klrb1c','Nkg7','Gzmb','Prf1','Xcl1'),
  inflammation_noncore=c('Tnf','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl2','Ccl3','Ccl4','Ccl5','Nfkb1','Rela','Tlr4','Tgfb1','Ptgs2','S100a8','S100a9','Nfkbia','Tnfaip3')
)
core_psps <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
marker_sets_filtered <- lapply(marker_sets,function(x)setdiff(unique(x),core_psps))
score_module <- function(mat,genes){genes<-intersect(genes,rownames(mat)); if(!length(genes)) return(rep(NA_real_,ncol(mat))); z<-t(scale(t(mat[genes,,drop=FALSE]))); z[!is.finite(z)]<-0; colMeans(z,na.rm=TRUE)}
coverage <- bind_rows(lapply(names(marker_sets),function(nm){req<-unique(marker_sets[[nm]]); filt<-marker_sets_filtered[[nm]]; present<-intersect(filt,rownames(expr)); tibble(module=nm,n_requested=length(req),n_core_psps_removed=sum(req%in%core_psps),removed_core_genes=paste(intersect(req,core_psps),collapse=';'),n_after_core_exclusion=length(filt),n_present_in_matrix=length(present),coverage=length(present)/length(filt),genes_present=paste(present,collapse=';'),genes_missing=paste(setdiff(filt,present),collapse=';'))}))
write.csv(coverage,file.path(outdir,'12wk_immune_marker_coverage.csv'),row.names=FALSE)
write.csv(bind_rows(lapply(names(marker_sets_filtered),function(nm)tibble(module=nm,gene=marker_sets_filtered[[nm]]))),file.path(outdir,'12wk_immune_marker_sets_filtered.csv'),row.names=FALSE)
mod <- names(marker_sets_filtered); sm <- sapply(marker_sets_filtered,function(g)score_module(logexpr,g)); colnames(sm)<-mod
scores <- data.frame(sample=colnames(expr),sm,PSPS_core=ps$PSPS_core,group=meta$group,check.names=FALSE)
write.csv(scores,file.path(outdir,'12wk_immune_sample_scores.csv'),row.names=FALSE)
# Spearman correlations against continuous PSPS_core; BH across modules.
cor_tbl <- bind_rows(lapply(mod,function(nm){ok<-is.finite(scores$PSPS_core)&is.finite(scores[[nm]]); ct<-suppressWarnings(cor.test(scores$PSPS_core[ok],scores[[nm]][ok],method='spearman',exact=FALSE)); tibble(module=nm,n=sum(ok),spearman_rho=unname(ct$estimate),p_value=ct$p.value)})); cor_tbl$padj_BH<-p.adjust(cor_tbl$p_value,'BH'); write.csv(cor_tbl,file.path(outdir,'12wk_immune_psps_spearman.csv'),row.names=FALSE)
# Shared five-group limma design and planned matched high-oleic contrasts.
design <- model.matrix(~0+group,data=meta); colnames(design)<-sub('^group','',colnames(design)); planned<-c(primary__HO_Safflower_vs_Safflower='HO_Safflower_12wk-Safflower_12wk',primary__HO_Sunflower_vs_Sunflower='HO_Sunflower_12wk-Sunflower_12wk')
fit0<-lmFit(t(as.matrix(scores[,mod,drop=FALSE])),design); contrast_tbl<-bind_rows(lapply(names(planned),function(nm){ff<-eBayes(contrasts.fit(fit0,makeContrasts(contrasts=planned[[nm]],levels=design)),robust=TRUE); topTable(ff,number=Inf,sort.by='none')%>%rownames_to_column('module')%>%mutate(contrast=nm)})); contrast_tbl$significant_BH<-contrast_tbl$adj.P.Val<.05; write.csv(contrast_tbl,file.path(outdir,'12wk_immune_group_contrasts.csv'),row.names=FALSE)
mc <- cor(scores[,mod,drop=FALSE],method='spearman',use='pairwise.complete.obs'); write.csv(mc,file.path(outdir,'12wk_immune_module_spearman_matrix.csv')); pdf(file.path(figdir,'12wk_immune_module_spearman_heatmap.pdf'),width=8,height=7); print(pheatmap(mc,display_numbers=TRUE,number_format='%.2f')); dev.off()
long<-scores%>%pivot_longer(cols=all_of(mod),names_to='module',values_to='score'); p<-ggplot(long,aes(group,score,fill=group))+geom_boxplot(outlier.shape=NA,alpha=.75)+geom_jitter(width=.12,size=1.4)+facet_wrap(~module,scales='free_y',ncol=3)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='none')+labs(title='12-week immune marker-module scores',x=NULL,y='Within-matrix gene-wise z-score mean'); ggsave(file.path(figdir,'12wk_immune_module_boxplots.pdf'),p,width=12,height=9)
p2<-ggplot(scores,aes(PSPS_core,leukocyte))+geom_point(size=2,alpha=.8,colour='#2C7FB8')+geom_smooth(method='lm',se=TRUE,colour='#D95F0E')+theme_classic(base_size=11)+labs(title='12-week PSPS_core and leukocyte module',x='PSPS_core',y='Leukocyte module score'); ggsave(file.path(figdir,'12wk_PSPS_vs_leukocyte.pdf'),p2,width=5.5,height=4.5)
manifest<-c('Analysis: 12-week GSE292815 mouse pancreatic bulk immune microenvironment association','Input: results/qc/GSE292815_KC_12wk_matrix.rds; log2(normalized expression + 1)','PSPS input: results/12wk_primary/scores/12wk_sample_scores.csv',paste0('Samples: ',ncol(expr),'; groups: ',paste(unique(meta$group),collapse=', ')),paste0('Core PSPS genes excluded before scoring: ',paste(core_psps,collapse=', ')),'Score: arithmetic mean of gene-wise z-scores for present marker genes','Correlation: Spearman with continuous PSPS_core; BH correction across modules','Group tests: same five-group limma design; planned matched high-oleic contrasts; BH-adjusted limma P values','Interpretation: marker scores are transcriptomic proxies, not direct cell counts or causal immune activation evidence.')
writeLines(manifest,file.path(outdir,'README.md')); writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo_12wk_immune.txt')); message('Completed immune analysis: ',outdir)
