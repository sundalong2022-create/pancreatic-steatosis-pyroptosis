#!/usr/bin/env Rscript
# 12-week PSPS and pathway-score analysis.
suppressPackageStartupMessages({library(limma); library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(pheatmap)})
root <- normalizePath('.', mustWork=TRUE)
base <- file.path(root, 'results/12wk_primary'); score_dir <- file.path(base, 'scores'); fig_dir <- file.path(base, 'figures')
dir.create(score_dir, recursive=TRUE, showWarnings=FALSE); dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds')); expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr)==21, all(colnames(expr)==meta$matrix_column))
logexpr <- log2(expr + 1)
gs_obj <- readRDS(file.path(root, 'scripts/results/gene_sets_mouse/gene_sets.rds')); pathways <- gs_obj$gene_sets; core_axis <- gs_obj$core_axis
score_sets <- c('pyroptosis','inflammasome','lipid_metabolism','oxidative_stress_ROS')
score_names <- c('PSPS_core','pyroptosis_score','inflammasome_score','lipid_metabolism_score','oxidative_stress_ROS_score')
score_matrix <- function(mat, genes) { genes <- intersect(unique(genes), rownames(mat)); if(!length(genes)) return(rep(NA_real_,ncol(mat))); z <- t(scale(t(mat[genes,,drop=FALSE]))); z[!is.finite(z)]<-0; colMeans(z,na.rm=TRUE) }
scores <- data.frame(sample=colnames(expr), stringsAsFactors=FALSE); scores$PSPS_core <- score_matrix(logexpr,core_axis)
for(nm in score_sets) scores[[paste0(nm,'_score')]] <- score_matrix(logexpr,pathways[[nm]])
scores <- left_join(scores,meta %>% select(matrix_column,group),by=c('sample'='matrix_column')); scores$dataset <- 'GSE292815_KC_12wk'
write.csv(scores,file.path(score_dir,'12wk_sample_scores.csv'),row.names=FALSE)
score_group_summary <- scores %>% pivot_longer(all_of(score_names),names_to='score',values_to='value') %>% group_by(group,score) %>% summarise(n=n(),mean=mean(value),sd=sd(value),median=median(value),.groups='drop')
write.csv(score_group_summary,file.path(score_dir,'12wk_score_group_summary.csv'),row.names=FALSE)
design <- model.matrix(~0+group,data=meta); colnames(design)<-sub('^group','',colnames(design))
planned <- c(primary__HO_Safflower_vs_Safflower='HO_Safflower_12wk-Safflower_12wk',primary__HO_Sunflower_vs_Sunflower='HO_Sunflower_12wk-Sunflower_12wk')
score_mat <- t(as.matrix(scores[,score_names])); rownames(score_mat)<-score_names; score_fit<-lmFit(score_mat,design)
score_stats <- bind_rows(lapply(names(planned),function(nm){
  fit <- eBayes(contrasts.fit(score_fit,makeContrasts(contrasts=planned[[nm]],levels=design)),robust=TRUE)
  tt <- topTable(fit,number=Inf,sort.by='none') %>% rownames_to_column('score')
  ## limma's moderated t statistic gives the standard error of each score contrast.
  ## Report a 95% CI for the difference in within-matrix z-score means.
  se <- fit$coefficients[,1] / fit$t[,1]
  df <- fit$df.total
  tt$SE <- unname(se[tt$score])
  tt$CI95_lower <- tt$logFC - qt(0.975, df=df) * tt$SE
  tt$CI95_upper <- tt$logFC + qt(0.975, df=df) * tt$SE
  tt %>% mutate(contrast=nm)
}))
score_stats$padj_BH <- ave(score_stats$P.Value,score_stats$contrast,FUN=function(x)p.adjust(x,'BH')); write.csv(score_stats,file.path(score_dir,'12wk_primary_score_contrasts.csv'),row.names=FALSE)
kw <- bind_rows(lapply(score_names,function(sc){kt<-kruskal.test(scores[[sc]]~scores$group); data.frame(score=sc,n=nrow(scores),n_groups=nlevels(factor(scores$group)),kruskal_chisq=unname(kt$statistic),kruskal_p=kt$p.value)}))
write.csv(kw,file.path(score_dir,'12wk_score_omnibus_tests.csv'),row.names=FALSE)
long <- scores %>% pivot_longer(all_of(score_names),names_to='score',values_to='value')
g <- ggplot(long,aes(group,value,fill=group))+geom_boxplot(outlier.shape=NA,alpha=.75)+geom_jitter(width=.13,size=1.5)+facet_wrap(~score,scales='free_y',ncol=3)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='none')+labs(title='12-week pathway and PSPS scores',x=NULL,y='Within-matrix gene-wise z-score mean')
ggsave(file.path(fig_dir,'12wk_primary_score_boxplots.pdf'),g,width=12,height=7)
pairs <- expand.grid(score1=score_names,score2=score_names,stringsAsFactors=FALSE) %>% filter(score1<score2)
cor_out <- bind_rows(lapply(seq_len(nrow(pairs)),function(i){a<-scores[[pairs$score1[i]]]; b<-scores[[pairs$score2[i]]]; ct<-suppressWarnings(cor.test(a,b,method='spearman',exact=FALSE)); data.frame(score1=pairs$score1[i],score2=pairs$score2[i],n=length(a),spearman_rho=unname(ct$estimate),p_value=ct$p.value)}))
cor_out$padj_BH<-p.adjust(cor_out$p_value,method='BH'); write.csv(cor_out,file.path(score_dir,'12wk_score_correlations_spearman.csv'),row.names=FALSE)
## Stratified check of the PSPS--ROS relationship by oil background. CD is
## excluded; each background combines its ordinary- and high-oleic-acid groups.
oil_background <- ifelse(grepl('Safflower', scores$group), 'Safflower_background',
                         ifelse(grepl('Sunflower', scores$group), 'Sunflower_background', NA_character_))
strat <- bind_rows(lapply(na.omit(unique(oil_background)), function(bg) {
  ix <- which(oil_background == bg); ct <- suppressWarnings(cor.test(scores$PSPS_core[ix], scores$oxidative_stress_ROS_score[ix], method='spearman', exact=FALSE))
  data.frame(oil_background=bg, n=length(ix), spearman_rho=unname(ct$estimate), p_value=ct$p.value)
}))
strat$padj_BH <- p.adjust(strat$p_value, method='BH')
write.csv(strat, file.path(score_dir,'12wk_psps_ros_by_oil_background.csv'), row.names=FALSE)
cor_mat<-cor(scores[,score_names],method='spearman',use='pairwise.complete.obs'); write.csv(cor_mat,file.path(score_dir,'12wk_score_correlation_matrix.csv')); pheatmap(cor_mat,display_numbers=TRUE,number_format='%.2f',filename=file.path(fig_dir,'12wk_score_correlation_heatmap.pdf'),width=7,height=6)
core_present<-intersect(core_axis,rownames(expr)); write.csv(data.frame(requested_core_gene=core_axis,present=core_axis%in%rownames(expr)),file.path(score_dir,'12wk_core_gene_presence.csv'),row.names=FALSE)
core_log<-logexpr[core_present,,drop=FALSE]; core_z<-t(scale(t(core_log))); core_z[!is.finite(core_z)]<-0
ann<-meta %>% select(matrix_column,group) %>% as.data.frame(); rownames(ann)<-ann$matrix_column; ann$matrix_column<-NULL
pheatmap(core_z[,rownames(ann),drop=FALSE],annotation_col=ann,cluster_cols=TRUE,cluster_rows=FALSE,filename=file.path(fig_dir,'12wk_core_gene_heatmap.pdf'),width=12,height=4.5,main='12-week core pyroptosis genes; row z-score')
core_long<-as.data.frame(core_log) %>% rownames_to_column('gene') %>% pivot_longer(-gene,names_to='sample',values_to='log2_expression') %>% left_join(meta %>% select(matrix_column,group),by=c('sample'='matrix_column'))
write.csv(core_long,file.path(score_dir,'12wk_core_gene_expression_long.csv'),row.names=FALSE); write.csv(core_long %>% group_by(gene,group) %>% summarise(n=n(),mean=mean(log2_expression),sd=sd(log2_expression),median=median(log2_expression),.groups='drop'),file.path(score_dir,'12wk_core_gene_group_summary.csv'),row.names=FALSE)
g<-ggplot(core_long,aes(group,log2_expression,fill=group))+geom_boxplot(outlier.shape=NA,alpha=.75)+geom_jitter(width=.13,size=1.2)+facet_wrap(~gene,scales='free_y',ncol=3)+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='none')+labs(title='12-week core pyroptosis genes',x=NULL,y='log2(normalized expression + 1)')
ggsave(file.path(fig_dir,'12wk_core_gene_grouped.pdf'),g,width=11,height=7)
core_primary<-data.frame(); for(i in seq_along(planned)){ff<-eBayes(contrasts.fit(lmFit(logexpr[core_present,,drop=FALSE],design),makeContrasts(contrasts=planned[[i]],levels=design)),robust=TRUE); core_primary<-bind_rows(core_primary,topTable(ff,number=Inf,sort.by='none') %>% rownames_to_column('gene_symbol') %>% mutate(contrast=names(planned)[i]))}; core_primary$significant_FDR<-core_primary$adj.P.Val<.05; write.csv(core_primary,file.path(score_dir,'12wk_core_gene_primary_contrasts.csv'),row.names=FALSE)
core_cor<-bind_rows(lapply(core_present,function(gene){ct<-suppressWarnings(cor.test(core_log[gene,],scores$PSPS_core,method='spearman',exact=FALSE)); data.frame(gene=gene,spearman_rho=unname(ct$estimate),p_value=ct$p.value)})); core_cor$padj_BH<-p.adjust(core_cor$p_value,'BH'); write.csv(core_cor,file.path(score_dir,'12wk_core_gene_vs_PSPS_correlations.csv'),row.names=FALSE)
for(nm in names(planned)){fgfile<-file.path(base,'enrichment',paste0(nm,'__fgsea.csv')); if(!file.exists(fgfile)) next; fg<-read.csv(fgfile,check.names=FALSE); genes<-unique(unlist(strsplit(fg$leadingEdge[fg$padj<.05],';',fixed=TRUE))); genes<-intersect(genes,rownames(expr)); if(!length(genes)) next; if(length(genes)>60) genes<-genes[order(apply(logexpr[genes,,drop=FALSE],1,var),decreasing=TRUE)[1:60]]; zz<-t(scale(t(logexpr[genes,,drop=FALSE]))); zz[!is.finite(zz)]<-0; pheatmap(zz[,rownames(ann),drop=FALSE],annotation_col=ann,cluster_cols=TRUE,cluster_rows=TRUE,show_rownames=length(genes)<=50,filename=file.path(fig_dir,paste0(nm,'__leading_edge_heatmap.pdf')),width=12,height=max(6,length(genes)*.12)); write.csv(data.frame(gene=genes),file.path(score_dir,paste0(nm,'__leading_edge_genes.csv')),row.names=FALSE)}
writeLines(capture.output(sessionInfo()),file.path(base,'sessionInfo_psps_12wk.txt')); message('Completed 12-week PSPS, score contrasts, correlations, core-gene plots, and leading-edge heatmaps.')
