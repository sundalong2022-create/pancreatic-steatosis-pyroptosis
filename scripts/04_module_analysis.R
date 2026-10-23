#!/usr/bin/env Rscript
# 12-week gene-set modules plus correlation-network module fallback.
suppressPackageStartupMessages({library(limma); library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(pheatmap); library(clusterProfiler); library(org.Mm.eg.db)})
root <- normalizePath('.', mustWork=TRUE); base <- file.path(root,'results/12wk_primary'); mod_dir <- file.path(base,'modules'); fig_dir <- file.path(base,'figures')
dir.create(mod_dir,recursive=TRUE,showWarnings=FALSE); dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(root,'results/qc/GSE292815_KC_12wk_matrix.rds')); expr<-obj$expr; meta<-obj$metadata; logexpr<-log2(expr+1)
scores <- read.csv(file.path(base,'scores/12wk_sample_scores.csv'),check.names=FALSE); rownames(scores)<-scores$sample
score_cols <- c('PSPS_core','pyroptosis_score','inflammasome_score','lipid_metabolism_score','oxidative_stress_ROS_score')
writeLines(c('method=correlation_distance_hierarchical_fallback','reason=WGCNA and dynamicTreeCut are unavailable','samples=21','scope=12-week KC matrix only'),file.path(mod_dir,'module_method.txt'))
kegg_internal <- requireNamespace('KEGG.db', quietly=TRUE); kegg_online <- requireNamespace('KEGGREST', quietly=TRUE); kegg_cache <- file.path(base,'enrichment','kegg_mmu_term2gene_name.rds'); kegg_map <- if(!kegg_internal && file.exists(kegg_cache)) tryCatch(readRDS(kegg_cache),error=function(e)NULL) else NULL; kegg_available <- kegg_internal || (!is.null(kegg_map) && nrow(kegg_map$term2gene)>0); use_internal_kegg <- kegg_internal
kegg_status <- if(kegg_internal) 'KEGG status: available via KEGG.db' else if(kegg_available) 'KEGG status: available via cached KEGGREST mapping' else if(kegg_online) 'KEGG status: KEGGREST installed; cached mapping unavailable, KEGG CSVs are empty' else 'KEGG status: unavailable; install KEGGREST or KEGG.db, KEGG CSVs are empty'
writeLines(kegg_status, file.path(mod_dir,'KEGG_status.txt'))
run_kegg_ora <- function(ids) {
  ids <- unique(as.character(ids)); if(!length(ids) || !kegg_available) return(NULL)
  if(kegg_internal) return(tryCatch(enrichKEGG(ids,organism='mmu',minGSSize=3,maxGSSize=2000,pAdjustMethod='BH',use_internal_data=TRUE),error=function(e)NULL))
  tryCatch(enricher(gene=ids,TERM2GENE=kegg_map$term2gene,TERM2NAME=kegg_map$term2name,minGSSize=3,maxGSSize=2000,pAdjustMethod='BH'),error=function(e)NULL)
}
module_scores <- scores %>% dplyr::select(sample,group,all_of(score_cols)); write.csv(module_scores,file.path(mod_dir,'gene_set_module_scores.csv'),row.names=FALSE)
cor_score <- cor(module_scores[,score_cols],method='spearman',use='pairwise.complete.obs'); write.csv(cor_score,file.path(mod_dir,'gene_set_module_correlations.csv')); pheatmap(cor_score,display_numbers=TRUE,number_format='%.2f',filename=file.path(fig_dir,'12wk_gene_set_module_correlations.pdf'),width=7,height=6)
keep <- rowSums(logexpr>1) >= 6; vv <- apply(logexpr[keep,,drop=FALSE],1,var); genes <- names(sort(vv,decreasing=TRUE))[seq_len(min(3000,sum(keep)))]
gene_cor <- cor(t(logexpr[genes,,drop=FALSE]),method='spearman',use='pairwise.complete.obs'); gene_dist<-as.dist(1-abs(gene_cor)); hc<-hclust(gene_dist,method='average')
ks <- 4:8; sil <- data.frame(k=ks,mean_silhouette=NA_real_)
if(requireNamespace('cluster',quietly=TRUE)) for(i in seq_along(ks)) sil$mean_silhouette[i]<-mean(cluster::silhouette(cutree(hc,k=ks[i]),gene_dist)[,'sil_width'])
if(all(is.na(sil$mean_silhouette))) selected_k<-6 else selected_k<-sil$k[which.max(sil$mean_silhouette)]
write.csv(sil,file.path(mod_dir,'module_cluster_k_selection.csv'),row.names=FALSE); writeLines(paste0('selected_k=',selected_k),file.path(mod_dir,'module_selected_k.txt'))
membership <- data.frame(gene=names(cutree(hc,k=selected_k)),module=paste0('M',cutree(hc,k=selected_k)),stringsAsFactors=FALSE); write.csv(membership,file.path(mod_dir,'coexpression_module_gene_membership.csv'),row.names=FALSE)
module_genes <- split(membership$gene,membership$module)
eig <- sapply(module_genes,function(gs){z<-prcomp(t(logexpr[gs,,drop=FALSE]),center=TRUE,scale.=FALSE)$x[,1]; if(cor(z,colMeans(logexpr[gs,,drop=FALSE]))<0) z<- -z; z})
eig <- as.data.frame(eig); rownames(eig)<-colnames(logexpr); colnames(eig)<-names(module_genes)
eig_out <- data.frame(sample=rownames(eig),group=meta$group,eig,check.names=FALSE); write.csv(eig_out,file.path(mod_dir,'coexpression_module_eigengenes.csv'),row.names=FALSE)
trait_long <- bind_rows(lapply(colnames(eig),function(m) bind_rows(lapply(score_cols,function(sc){ct<-suppressWarnings(cor.test(eig[[m]],scores[[sc]],method='spearman',exact=FALSE)); data.frame(module=m,trait=sc,n=21,rho=unname(ct$estimate),p_value=ct$p.value)}))))
trait_long$padj_BH<-p.adjust(trait_long$p_value,'BH'); write.csv(trait_long,file.path(mod_dir,'module_trait_spearman_correlations.csv'),row.names=FALSE)
design<-model.matrix(~0+group,data=meta); colnames(design)<-sub('^group','',colnames(design)); planned<-c(primary__HO_Safflower_vs_Safflower='HO_Safflower_12wk-Safflower_12wk',primary__HO_Sunflower_vs_Sunflower='HO_Sunflower_12wk-Sunflower_12wk')
eig_mat<-t(as.matrix(eig)); rownames(eig_mat)<-colnames(eig); eig_contrast<-bind_rows(lapply(names(planned),function(nm){ff<-eBayes(contrasts.fit(limma::lmFit(eig_mat,design),makeContrasts(contrasts=planned[[nm]],levels=design)),robust=TRUE); topTable(ff,number=Inf,sort.by='none') %>% rownames_to_column('module') %>% mutate(contrast=nm)})); write.csv(eig_contrast,file.path(mod_dir,'module_primary_contrasts.csv'),row.names=FALSE)
stability <- bind_rows(lapply(names(module_genes),function(m){gs<-module_genes[[m]]; full<-eig[[m]]; cors<-vapply(seq_len(ncol(logexpr)),function(i){z<-logexpr[gs,-i,drop=FALSE]; ee<-prcomp(t(z),center=TRUE,scale.=FALSE)$x[,1]; if(cor(ee,colMeans(z))<0) ee<- -ee; cor(full[-i],ee)},numeric(1)); data.frame(module=m,n_genes=length(gs),median_leave_one_out_rho=median(cors),min_leave_one_out_rho=min(cors))}))
write.csv(stability,file.path(mod_dir,'module_leave_one_out_stability.csv'),row.names=FALSE)
top_modules <- unique(c(trait_long$module[order(-abs(trait_long$rho))[seq_len(min(3,nrow(trait_long)))]],trait_long$module[trait_long$trait=='PSPS_core' & abs(trait_long$rho)>=.5]))
for(m in top_modules){gs<-module_genes[[m]]; ids<-suppressMessages(bitr(gs,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Mm.eg.db))$ENTREZID; go<-tryCatch(enrichGO(ids,OrgDb=org.Mm.eg.db,ont='BP',pAdjustMethod='BH',readable=TRUE),error=function(e)NULL); if(is.null(go)) write.csv(data.frame(),file.path(mod_dir,paste0(m,'__GO_BP.csv')),row.names=FALSE) else write.csv(as.data.frame(go),file.path(mod_dir,paste0(m,'__GO_BP.csv')),row.names=FALSE); ek<-run_kegg_ora(ids); if(is.null(ek)) write.csv(data.frame(),file.path(mod_dir,paste0(m,'__KEGG.csv')),row.names=FALSE) else write.csv(as.data.frame(ek),file.path(mod_dir,paste0(m,'__KEGG.csv')),row.names=FALSE)}
writeLines(capture.output(sessionInfo()),file.path(base,'sessionInfo_modules_12wk.txt')); message('Completed 12-week gene-set modules and exploratory correlation-network modules.')
