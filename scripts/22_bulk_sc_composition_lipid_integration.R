#!/usr/bin/env Rscript
# Integrate mouse bulk lipid/PSPS scoring with xCell/EPIC adjustment and
# GSE154778 cell-level / pseudobulk composition decomposition.
suppressPackageStartupMessages({
  library(UCell); library(AUCell); library(data.table); library(dplyr)
  library(tidyr); library(tibble); library(ggplot2); library(pheatmap)
  library(lme4); library(msigdbr); library(xCell); library(EPIC)
  library(babelgene)
})
set.seed(292815)
options(stringsAsFactors=FALSE)
root <- normalizePath('.', mustWork=TRUE)
out <- file.path(root,'results','composition_lipid_integration')
fig <- file.path(out,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE); dir.create(fig,recursive=TRUE,showWarnings=FALSE)

core_mouse <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
core_human <- toupper(core_mouse)
spec <- tibble(
  score=c('Hallmark_FA_Metabolism','Hallmark_Adipogenesis','Reactome_FA_Metabolism','PPAR_Pathway'),
  collection=c('H','H','C2','C2'),
  subcollection=c(NA_character_,NA_character_,'CP:REACTOME','CP:KEGG_LEGACY'),
  gs_name=c('HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_ADIPOGENESIS','REACTOME_FATTY_ACID_METABOLISM','KEGG_PPAR_SIGNALING_PATHWAY'),
  concept=c('fatty-acid metabolism','adipogenesis','fatty-acid metabolism','PPAR signalling')
)
fetch_set <- function(i, species) {
  args <- list(db_species='HS',species=species,collection=spec$collection[i])
  if (!is.na(spec$subcollection[i])) args$subcollection <- spec$subcollection[i]
  tb <- do.call(msigdbr::msigdbr,args)
  hit <- tb %>% filter(gs_name==spec$gs_name[i])
  if (!nrow(hit)) stop('Missing gene set: ',spec$gs_name[i])
  genes <- unique(as.character(hit$gene_symbol))
  list(genes=genes, table=hit)
}
set_mouse_obj <- lapply(seq_len(nrow(spec)),fetch_set,species='Mus musculus')
set_human_obj <- lapply(seq_len(nrow(spec)),fetch_set,species='Homo sapiens')
names(set_mouse_obj) <- spec$score; names(set_human_obj) <- spec$score
sets_mouse <- lapply(set_mouse_obj,'[[','genes'); sets_human <- lapply(set_human_obj,'[[','genes')
names(sets_mouse) <- names(sets_human) <- spec$score
sets_mouse$PSPS_core <- core_mouse; sets_human$PSPS_core <- core_human
set_order <- c('PSPS_core',spec$score)

gene_prov <- bind_rows(lapply(seq_len(nrow(spec)),function(i){
  m <- set_mouse_obj[[i]]; h <- set_human_obj[[i]]
  tibble(score=spec$score[i],concept=spec$concept[i],gs_name=spec$gs_name[i],
    collection=spec$collection[i],subcollection=spec$subcollection[i],
    db_version=paste(unique(m$table$db_version),collapse=';'),
    source_species=paste(unique(m$table$gs_source_species),collapse=';'),
    target_mouse_n=length(m$genes),target_human_n=length(h$genes),
    mapping='MSigDB human definitions queried via msigdbr; mouse gene symbols returned by msigdbr species=Mus musculus')
}))
gene_prov <- bind_rows(
  tibble(score='PSPS_core',concept='PSPS six-gene core',gs_name='manual',
    collection=NA_character_,subcollection=NA_character_,db_version=NA_character_,
    source_species='manual mouse symbols',target_mouse_n=length(core_mouse),target_human_n=length(core_human),mapping='mouse core genes used directly'),
  gene_prov)
write.csv(gene_prov,file.path(out,'gene_set_provenance.csv'),row.names=FALSE)

# Seurat AddModuleScore-compatible implementation.
add_module_score <- function(logmat, genes, nbin=24L, ctrl=100L, seed=292815) {
  genes <- intersect(unique(genes),rownames(logmat))
  if (!length(genes)) return(rep(NA_real_,ncol(logmat)))
  avg <- rowMeans(logmat,na.rm=TRUE)
  ord <- order(avg,na.last=NA)
  bins <- integer(length(avg)); bins[ord] <- pmin(nbin,ceiling(seq_along(ord)/length(ord)*nbin))
  names(bins) <- names(avg)
  target <- colMeans(logmat[genes,,drop=FALSE],na.rm=TRUE)
  ctrl_genes <- character()
  set.seed(seed)
  for (g in genes) {
    pool <- names(bins)[bins==bins[g] & !(names(bins) %in% genes)]
    if (!length(pool)) pool <- setdiff(names(bins),genes)
    ctrl_genes <- c(ctrl_genes,sample(pool,size=min(ctrl,length(pool)),replace=length(pool)<ctrl))
  }
  ctrl_genes <- unique(ctrl_genes)
  ctrl_score <- if(length(ctrl_genes)) colMeans(logmat[ctrl_genes,,drop=FALSE],na.rm=TRUE) else 0
  as.numeric(target-ctrl_score)
}
fisher_ci <- function(r,n) {
  if (!is.finite(r) || n<=3 || abs(r)>=1) return(c(lower=NA_real_,upper=NA_real_))
  z <- atanh(r); se <- 1/sqrt(n-3)
  c(lower=tanh(z-1.96*se),upper=tanh(z+1.96*se))
}
spearman_row <- function(x,y,method_label='',scope_label='') {
  ok <- is.finite(x)&is.finite(y); n <- sum(ok)
  if(n<4 || sd(x[ok])==0 || sd(y[ok])==0) return(tibble(method=method_label,scope=scope_label,n=n,rho=NA_real_,p_value=NA_real_,CI95_low=NA_real_,CI95_high=NA_real_))
  ct <- suppressWarnings(cor.test(x[ok],y[ok],method='spearman',exact=FALSE))
  ci <- fisher_ci(unname(ct$estimate),n)
  tibble(method=method_label,scope=scope_label,n=n,rho=unname(ct$estimate),p_value=ct$p.value,CI95_low=ci[['lower']],CI95_high=ci[['upper']])
}
partial_multi <- function(x,y,covars,panel,score_method,psps_name,lipid_name) {
  covars <- as.data.frame(covars); ok <- is.finite(x)&is.finite(y)&apply(covars,1,function(z)all(is.finite(z)))
  n <- sum(ok); k <- ncol(covars)
  if(n <= k+3) return(tibble(panel=panel,score_method=score_method,psps_score=psps_name,lipid_score=lipid_name,n=n,n_covariates=k,df=NA_real_,rho=NA_real_,p_value=NA_real_))
  xr <- rank(x[ok],ties.method='average'); yr <- rank(y[ok],ties.method='average'); z <- as.matrix(covars[ok,,drop=FALSE])
  X <- cbind(1,apply(z,2,rank,ties.method='average'))
  rx <- lm.fit(X,xr)$residuals; ry <- lm.fit(X,yr)$residuals
  r <- suppressWarnings(cor(rx,ry)); df <- n-ncol(X)-1
  p <- if(is.finite(r)&&abs(r)<1&&df>0) 2*pt(abs(r)*sqrt(df/(1-r^2)),df=df,lower.tail=FALSE) else NA_real_
  tibble(panel=panel,score_method=score_method,psps_score=psps_name,lipid_score=lipid_name,n=n,n_covariates=k,df=df,rho=r,p_value=p)
}
score_matrix <- function(mat, sets, method_prefix) {
  feats <- lapply(sets,function(g)intersect(g,rownames(mat)))
  u <- ScoreSignatures_UCell(matrix=mat,features=feats,maxRank=min(1500L,nrow(mat)),chunk.size=200L,ncores=1L,ties.method='average')
  u <- as.data.frame(u); u <- u[,paste0(names(feats),'_UCell'),drop=FALSE]
  names(u) <- paste0(method_prefix,'_',names(feats))
  u
}
coverage_table <- function(sets,mat,context) {
  bind_rows(lapply(names(sets),function(nm) {
    req <- unique(sets[[nm]]); present <- intersect(req,rownames(mat))
    tibble(context=context,score=nm,requested_n=length(req),present_n=length(present),coverage=length(present)/length(req),genes_present=paste(present,collapse=';'),genes_missing=paste(setdiff(req,present),collapse=';'))
  }))
}

# --------------------------- mouse bulk -------------------------------------
obj <- readRDS(file.path(root,'results','qc','GSE292815_KC_12wk_matrix.rds'))
expr <- as.matrix(obj$expr); meta <- as.data.frame(obj$metadata)
logexpr <- log2(expr+1)
bulk_cov <- coverage_table(sets_mouse,expr,'mouse_bulk')
bulk_add <- as.data.frame(lapply(sets_mouse,function(g)add_module_score(logexpr,g)),check.names=FALSE)
names(bulk_add) <- paste0('AddModuleScore_',names(sets_mouse))
bulk_uc <- score_matrix(expr,sets_mouse,'UCell')
bulk_scores <- tibble(sample=colnames(expr),group=meta$group) %>% bind_cols(bulk_uc[match(colnames(expr),rownames(bulk_uc)),,drop=FALSE]) %>% bind_cols(bulk_add)
# Include frozen PSPS score for audit.
old_ps <- read.csv(file.path(root,'results','12wk_primary','scores','12wk_sample_scores.csv'),check.names=FALSE)
bulk_scores$frozen_PSPS_core <- old_ps$PSPS_core[match(bulk_scores$sample,old_ps$sample)]
write.csv(bulk_scores,file.path(out,'mouse_bulk_scores_UCell_AddModuleScore.csv'),row.names=FALSE)
write.csv(bulk_cov,file.path(out,'gene_set_coverage_mouse_bulk.csv'),row.names=FALSE)

# Human ortholog map for mouse bulk xCell/EPIC.
map <- babelgene::orthologs(rownames(expr),species='mouse',human=FALSE,top=TRUE) %>% filter(!is.na(human_symbol),symbol %in% rownames(expr)) %>% arrange(desc(support_n)) %>% distinct(human_symbol,.keep_all=TRUE)
map <- map[!duplicated(map$symbol),]
hm <- expr[match(map$symbol,rownames(expr)),,drop=FALSE]; rownames(hm) <- map$human_symbol
x <- xCell::xCellAnalysis(hm,rnaseq=TRUE,parallel.sz=1,parallel.type='SOCK')
x_tbl <- as.data.frame(t(x),check.names=FALSE) %>% rownames_to_column('sample')
e <- EPIC::EPIC(hm,withOtherCells=TRUE)
e_tbl <- as.data.frame(e$cellFractions,check.names=FALSE) %>% rownames_to_column('sample')
if(!is.null(e$otherCells)) e_tbl$EPIC_otherCells <- as.numeric(e$otherCells[match(e_tbl$sample,rownames(e$otherCells)),1])
x_get <- function(nm,alts=nm) {h<-intersect(alts,names(x_tbl)); if(length(h)) as.numeric(x_tbl[[h[1]]]) else rep(NA_real_,nrow(x_tbl))}
e_get <- function(nm,alts=nm) {h<-intersect(alts,names(e_tbl)); if(length(h)) as.numeric(e_tbl[[h[1]]]) else rep(NA_real_,nrow(e_tbl))}
deconv <- tibble(sample=colnames(expr),
  xCell_Macrophages=x_get('xCell_Macrophages',c('Macrophages')),
  xCell_Adipocytes=x_get('xCell_Adipocytes',c('Adipocytes')),
  xCell_Fibroblasts=x_get('xCell_Fibroblasts',c('Fibroblasts')),
  xCell_Epithelial=x_get('xCell_Epithelial',c('Epithelial cells','Epithelial')),
  xCell_ImmuneScore=x_get('xCell_ImmuneScore',c('ImmuneScore')),
  xCell_StromaScore=x_get('xCell_StromaScore',c('StromaScore')),
  EPIC_Macrophages=e_get('EPIC_Macrophages',c('Macrophages')),
  EPIC_CAFs=e_get('EPIC_CAFs',c('CAFs')),
  EPIC_Bcells=e_get('EPIC_Bcells',c('Bcells')),
  EPIC_Endothelial=e_get('EPIC_Endothelial',c('Endothelial')),
  EPIC_otherCells=e_get('EPIC_otherCells',c('EPIC_otherCells','otherCells')))
write.csv(deconv,file.path(out,'mouse_bulk_xCell_EPIC_estimates.csv'),row.names=FALSE)
write.csv(map,file.path(out,'mouse_to_human_ortholog_map.csv'),row.names=FALSE)
# xCell scores are enrichment scores; EPIC columns are fractions.
writeLines(c('Mouse bulk deconvolution used babelgene orthologs(mouse -> human) with top=TRUE, one human symbol per mapped gene.',
  'xCell output is an enrichment score, not a calibrated cell fraction.',
  'EPIC output is a constrained cell-fraction estimate; otherCells is residual fraction.',
  'These are cross-species ortholog-mapped proxies and are used as covariates, not as direct ground-truth cell counts.'),file.path(out,'deconvolution_method_note.txt'))

bulk_merge <- left_join(bulk_scores,deconv,by='sample')
write.csv(bulk_merge,file.path(out,'mouse_bulk_scores_with_xCell_EPIC.csv'),row.names=FALSE)
panels <- list(
  xCell_core=c('xCell_Macrophages','xCell_Adipocytes','xCell_Fibroblasts','xCell_Epithelial'),
  EPIC_core=c('EPIC_Macrophages','EPIC_CAFs','EPIC_otherCells'),
  reduced_combined=c('xCell_Macrophages','xCell_Fibroblasts','EPIC_Macrophages','EPIC_CAFs')
)
bulk_partial <- bind_rows(lapply(names(panels),function(pn)bind_rows(lapply(c('UCell','AddModuleScore'),function(sm)bind_rows(lapply(spec$score,function(lp){
  psn <- paste0(sm,'_PSPS_core'); lpn <- paste0(sm,'_',lp)
  partial_multi(bulk_merge[[psn]],bulk_merge[[lpn]],bulk_merge[,panels[[pn]],drop=FALSE],pn,sm,'PSPS_core',lp)
}))))))
bulk_partial$q_BH <- p.adjust(bulk_partial$p_value,'BH')
write.csv(bulk_partial,file.path(out,'mouse_bulk_psps_lipid_partial_deconvolution.csv'),row.names=FALSE)

# --------------------------- GSE154778 single cell ---------------------------
sc_dir <- file.path(root,'results','sc_location_GSE154778')
sc_counts <- as.matrix(fread(file.path(root,'data','raw','GSE154778','GSE154778_dgeMtx.csv.gz'),data.table=FALSE,check.names=FALSE,showProgress=FALSE))
sc_gene_col <- colnames(sc_counts)[1]
rownames(sc_counts) <- toupper(sub('\\..*$','',sc_counts[,sc_gene_col]))
sc_counts <- sc_counts[,-1,drop=FALSE]; mode(sc_counts) <- 'numeric'
if(anyDuplicated(rownames(sc_counts))) sc_counts <- rowsum(sc_counts,group=rownames(sc_counts),reorder=FALSE)
sc0 <- fread(file.path(sc_dir,'per_cell_scores.csv'),data.table=FALSE,check.names=FALSE)
sc_id_col <- if('Unnamed: 0' %in% names(sc0)) 'Unnamed: 0' else names(sc0)[1]
sc_meta <- sc0; sc_meta$cell_id <- as.character(sc_meta[[sc_id_col]]); sc_meta <- sc_meta %>% select(-all_of(sc_id_col)); rownames(sc_meta) <- sc_meta$cell_id
common <- intersect(colnames(sc_counts),sc_meta$cell_id); sc_counts <- sc_counts[,common,drop=FALSE]; sc_meta <- sc_meta[common,,drop=FALSE]
sc_meta$cell_type <- as.character(sc_meta$cell_type); sc_meta$sample_id <- as.character(sc_meta$sample_id)
sc_meta$sample_class <- factor(ifelse(grepl('^MET',sc_meta$sample_id),'metastatic','primary'),levels=c('primary','metastatic'))
sets_human <- sets_human[set_order]
sc_cov <- coverage_table(sets_human,sc_counts,'GSE154778_single_cell_any_detected')
# UCell cell scores
sc_uc <- score_matrix(sc_counts,sets_human,'UCell')
rownames(sc_uc) <- colnames(sc_counts)
# AUCell sensitivity
sc_auc <- AUCell_buildRankings(sc_counts,plotStats=FALSE,nCores=1,verbose=FALSE)
sc_auc_obj <- AUCell_calcAUC(sets_human,sc_auc,aucMaxRank=ceiling(0.05*nrow(sc_counts)),nCores=1,verbose=FALSE)
sc_auc_mat <- as.data.frame(t(getAUC(sc_auc_obj)),check.names=FALSE)
sc_auc_mat <- sc_auc_mat[colnames(sc_counts),,drop=FALSE]
names(sc_auc_mat) <- paste0('AUCell_',names(sc_auc_mat))
# AddModuleScore-compatible on logCPM
sc_lib <- colSums(sc_counts); sc_logcpm <- log1p(sweep(sc_counts,2,pmax(sc_lib,1),'/')*1e4)
sc_add <- as.data.frame(lapply(sets_human,function(g)add_module_score(sc_logcpm,g)),check.names=FALSE)
names(sc_add) <- paste0('AddModuleScore_',names(sets_human))
sc_scores <- tibble(cell_id=colnames(sc_counts)) %>% bind_cols(sc_meta[match(colnames(sc_counts),sc_meta$cell_id),c('sample_id','sample_class','cell_type','assignment_margin','low_information')]) %>% bind_cols(sc_uc) %>% bind_cols(sc_add) %>% bind_cols(sc_auc_mat)
write.csv(sc_scores,file.path(out,'sc_per_cell_scores_UCell_AddModuleScore_AUCell.csv'),row.names=FALSE)
sc_cov2 <- coverage_table(sets_human,sc_counts,'GSE154778_single_cell')
write.csv(bind_rows(bulk_cov,sc_cov,sc_cov2),file.path(out,'gene_set_coverage_all.csv'),row.names=FALSE)
# Cell-type heatmap summaries (major types; beta excluded from main heatmap).
major_types <- c('macrophage','acinar','stellate','ductal','epithelial_tumor')
uc_cols <- paste0('UCell_',set_order)
add_cols <- paste0('AddModuleScore_',set_order)
summary_uc <- sc_scores %>% filter(cell_type %in% major_types) %>% group_by(cell_type) %>% summarise(across(all_of(uc_cols),list(median=~median(.x,na.rm=TRUE),positive_rate=~mean(.x>0,na.rm=TRUE))),n_cells=n(),n_samples=n_distinct(sample_id),.groups='drop')
summary_add <- sc_scores %>% filter(cell_type %in% major_types) %>% group_by(cell_type) %>% summarise(across(all_of(add_cols),list(median=~median(.x,na.rm=TRUE),positive_rate=~mean(.x>0,na.rm=TRUE))),n_cells=n(),n_samples=n_distinct(sample_id),.groups='drop')
auc_cols <- paste0('AUCell_',set_order)
summary_auc <- sc_scores %>% filter(cell_type %in% major_types) %>% group_by(cell_type) %>% summarise(across(all_of(auc_cols),list(median=~median(.x,na.rm=TRUE),positive_rate=~mean(.x>0,na.rm=TRUE))),n_cells=n(),n_samples=n_distinct(sample_id),.groups='drop')
write.csv(summary_uc,file.path(out,'sc_cell_type_summary_UCell.csv'),row.names=FALSE); write.csv(summary_add,file.path(out,'sc_cell_type_summary_AddModuleScore.csv'),row.names=FALSE); write.csv(summary_auc,file.path(out,'sc_cell_type_summary_AUCell.csv'),row.names=FALSE)
mat_med <- as.matrix(summary_uc[,paste0(rep(uc_cols,each=1),'_median'),drop=FALSE]); rownames(mat_med)<-summary_uc$cell_type; colnames(mat_med)<-sub('UCell_','',uc_cols)
mat_pos <- as.matrix(summary_uc[,paste0(uc_cols,'_positive_rate'),drop=FALSE]); rownames(mat_pos)<-summary_uc$cell_type; colnames(mat_pos)<-sub('UCell_','',uc_cols)
if(nrow(mat_med)>0) pheatmap(mat_med,cluster_rows=FALSE,cluster_cols=FALSE,filename=file.path(fig,'sc_celltype_score_heatmap_UCell_median.pdf'),main='GSE154778 UCell median scores')
if(nrow(mat_pos)>0) pheatmap(mat_pos,cluster_rows=FALSE,cluster_cols=FALSE,filename=file.path(fig,'sc_celltype_score_heatmap_UCell_positive_rate.pdf'),main='GSE154778 UCell positive rate')
make_heatmaps <- function(summary_df, cols, prefix, title_prefix) {
  med <- as.matrix(summary_df[,paste0(cols,'_median'),drop=FALSE]); rownames(med) <- summary_df$cell_type; colnames(med) <- sub(paste0('^',prefix,'_'),' ',cols)
  pos <- as.matrix(summary_df[,paste0(cols,'_positive_rate'),drop=FALSE]); rownames(pos) <- summary_df$cell_type; colnames(pos) <- sub(paste0('^',prefix,'_'),' ',cols)
  if(nrow(med)>0) pheatmap(med,cluster_rows=FALSE,cluster_cols=FALSE,filename=file.path(fig,paste0('sc_celltype_score_heatmap_',prefix,'_median.pdf')),main=paste(title_prefix,'median scores'))
  if(nrow(pos)>0) pheatmap(pos,cluster_rows=FALSE,cluster_cols=FALSE,filename=file.path(fig,paste0('sc_celltype_score_heatmap_',prefix,'_positive_rate.pdf')),main=paste(title_prefix,'positive rate'))
}
make_heatmaps(summary_add,add_cols,'AddModuleScore','GSE154778 AddModuleScore')
make_heatmaps(summary_auc,auc_cols,'AUCell','GSE154778 AUCell')
# Mixed models. Cell-level results are support only; sample x cell-type pseudobulk is primary.
score_pairs <- spec$score
fit_lmer <- function(dat,ps_col,lp_col,scope) {
  dat <- dat %>% filter(is.finite(.data[[ps_col]]),is.finite(.data[[lp_col]]))
  if(nrow(dat)<50 || n_distinct(dat$sample_id)<3) return(tibble())
  dat$sample_id <- factor(dat$sample_id); dat$sample_class <- factor(dat$sample_class)
  f <- if(scope=='all_major_types_unadjusted') as.formula(paste(ps_col,'~',lp_col,'+(1|sample_id)')) else if(n_distinct(dat$cell_type)>1 && n_distinct(dat$sample_class)>1) as.formula(paste(ps_col,'~',lp_col,'+cell_type+sample_class+(1|sample_id)')) else if(n_distinct(dat$cell_type)>1) as.formula(paste(ps_col,'~',lp_col,'+cell_type+(1|sample_id)')) else if(n_distinct(dat$sample_class)>1) as.formula(paste(ps_col,'~',lp_col,'+sample_class+(1|sample_id)')) else as.formula(paste(ps_col,'~',lp_col,'+(1|sample_id)'))
  fit <- tryCatch(lmer(f,data=dat,REML=TRUE),error=function(e)NULL); if(is.null(fit)) return(tibble())
  co <- summary(fit)$coefficients; if(!(lp_col %in% rownames(co))) return(tibble())
  est <- co[lp_col,'Estimate']; se <- co[lp_col,'Std. Error']; p <- 2*pnorm(abs(est/se),lower.tail=FALSE)
  tibble(scope=scope,psps_method=sub('_PSPS_core','',ps_col),lipid_method=sub(paste0('_',lp_col), '', lp_col),psps_score=ps_col,lipid_score=lp_col,n_cells=nrow(dat),n_samples=n_distinct(dat$sample_id),estimate=est,SE=se,CI95_low=est-1.96*se,CI95_high=est+1.96*se,p_value=p)
}
sc_model <- bind_rows(lapply(c('UCell','AddModuleScore'),function(sm)bind_rows(lapply(score_pairs,function(lp){psn<-paste0(sm,'_PSPS_core'); lpn<-paste0(sm,'_',lp); bind_rows(
  fit_lmer(sc_scores %>% filter(cell_type %in% major_types),psn,lpn,'all_major_types'),
  fit_lmer(sc_scores %>% filter(cell_type %in% major_types),psn,lpn,'all_major_types_unadjusted'),
  bind_rows(lapply(major_types,function(ct)fit_lmer(sc_scores %>% filter(cell_type==ct),psn,lpn,ct)))
)})))) %>% group_by(psps_method,scope) %>% mutate(q_BH=p.adjust(p_value,'BH')) %>% ungroup()
write.csv(sc_model,file.path(out,'sc_cell_level_mixed_models.csv'),row.names=FALSE)

# Pseudobulk raw counts by sample x cell type, >=10 cells.
sc_unit_counts <- sc_meta %>% count(sample_id,cell_type,name='n_cells') %>% filter(n_cells>=10)
sc_units <- paste(sc_unit_counts$sample_id,sc_unit_counts$cell_type,sep='|||')
sc_key <- paste(sc_meta$sample_id,sc_meta$cell_type,sep='|||')
keep <- sc_key %in% sc_units
pb_counts <- t(rowsum(t(sc_counts[,keep,drop=FALSE]),group=sc_key[keep],reorder=FALSE))
pb_unit_key <- colnames(pb_counts)
pb_meta <- tibble(unit_key=pb_unit_key) %>% separate(unit_key,c('sample_id','cell_type'),sep='\\|\\|\\|',remove=FALSE) %>% left_join(sc_unit_counts,by=c('sample_id','cell_type')) %>% mutate(sample_class=factor(ifelse(grepl('^MET',sample_id),'metastatic','primary'),levels=c('primary','metastatic')))
pb_uc <- score_matrix(pb_counts,sets_human,'UCell'); rownames(pb_uc)<-pb_unit_key
pb_lib <- colSums(pb_counts); pb_logcpm <- log1p(sweep(pb_counts,2,pmax(pb_lib,1),'/')*1e4)
pb_add <- as.data.frame(lapply(sets_human,function(g)add_module_score(pb_logcpm,g)),check.names=FALSE); names(pb_add)<-paste0('AddModuleScore_',names(sets_human)); rownames(pb_add)<-pb_unit_key
pb_scores <- pb_meta %>% bind_cols(as.data.frame(pb_uc[pb_unit_key,,drop=FALSE])) %>% bind_cols(pb_add[pb_unit_key,,drop=FALSE])
write.csv(pb_scores,file.path(out,'sc_sample_celltype_pseudobulk_scores.csv'),row.names=FALSE)

# Composition decomposition: weighted W, equal U, and fixed mean-composition F.
decomp_types <- major_types
pb_major <- pb_scores %>% filter(cell_type %in% decomp_types)
prop_tbl <- pb_major %>% group_by(sample_id) %>% mutate(prop=n_cells/sum(n_cells)) %>% ungroup()
mean_prop <- prop_tbl %>% group_by(cell_type) %>% summarise(ref_prop=mean(prop),.groups='drop')
composition_scores <- bind_rows(lapply(c('UCell','AddModuleScore'),function(sm) {
  bind_rows(lapply(set_order,function(sc) {
    col <- paste0(sm,'_',sc)
    if(!(col %in% names(prop_tbl))) return(NULL)
    w <- prop_tbl %>% group_by(sample_id) %>% summarise(W=sum(prop*.data[[col]],na.rm=TRUE),U=mean(.data[[col]],na.rm=TRUE),.groups='drop')
    f <- prop_tbl %>% left_join(mean_prop,by='cell_type') %>% group_by(sample_id) %>% summarise(F=sum(ref_prop*.data[[col]],na.rm=TRUE)/sum(ref_prop),.groups='drop')
    w %>% full_join(f,by='sample_id') %>% mutate(score_method=sm,score_name=sc)
  }))
}))
write.csv(composition_scores,file.path(out,'sc_sample_composition_W_U_F_scores.csv'),row.names=FALSE)

decomp_results <- bind_rows(lapply(c('UCell','AddModuleScore'),function(sm)bind_rows(lapply(spec$score,function(lp){
  for_scope <- function(scope) {
    d <- composition_scores %>% filter(score_method==sm,score_name %in% c('PSPS_core',lp))
  }
  ps <- composition_scores %>% filter(score_method==sm,score_name=='PSPS_core') %>% select(sample_id, W_psps=W,U_psps=U,F_psps=F)
  li <- composition_scores %>% filter(score_method==sm,score_name==lp) %>% select(sample_id,W_lipid=W,U_lipid=U,F_lipid=F)
  d <- left_join(ps,li,by='sample_id') %>% left_join(sc_meta %>% distinct(sample_id,sample_class),by='sample_id')
  bind_rows(lapply(c('pooled','primary','metastatic'),function(sc) {
    z <- if(sc=='pooled') d else d %>% filter(sample_class==sc)
    outz <- bind_rows(lapply(c('W','U','F'),function(mix) {
      x <- z[[paste0(mix,'_psps')]]; y <- z[[paste0(mix,'_lipid')]]
      q <- spearman_row(x,y,paste0(sm,'_',mix),'')
      tibble(score_method=sm,lipid_score=lp,scope=sc,mixing=mix,n=q$n,rho=q$rho,p_value=q$p_value,CI95_low=q$CI95_low,CI95_high=q$CI95_high)
    }))
    rw <- outz$rho[outz$mixing=='W']; ru <- outz$rho[outz$mixing=='U']
    outz %>% mutate(delta_abs_rho=ifelse(mixing=='W',abs(rw)-abs(ru),NA_real_))
  }))
}))))
decomp_results$q_BH <- p.adjust(decomp_results$p_value,'BH')
write.csv(decomp_results,file.path(out,'sc_composition_decomposition_correlations.csv'),row.names=FALSE)

# Partial correlation of weighted pseudobulk scores adjusted for the five major cell proportions.
prop_wide <- prop_tbl %>% select(sample_id,cell_type,prop) %>% pivot_wider(names_from=cell_type,values_from=prop,values_fill=0)
names(prop_wide)[-1] <- paste0('prop_',names(prop_wide)[-1])
sc_partial <- bind_rows(lapply(c('UCell','AddModuleScore'),function(sm)bind_rows(lapply(spec$score,function(lp){
  ps <- composition_scores %>% filter(score_method==sm,score_name=='PSPS_core') %>% select(sample_id,psps=W)
  li <- composition_scores %>% filter(score_method==sm,score_name==lp) %>% select(sample_id,lipid=W)
  d <- ps %>% inner_join(li,by='sample_id') %>% inner_join(prop_wide,by='sample_id') %>% inner_join(sc_meta %>% distinct(sample_id,sample_class),by='sample_id')
  cov <- d %>% select(starts_with('prop_'))
  base_res <- partial_multi(d$psps,d$lipid,cov,'major5_cell_proportions',sm,'PSPS_core',lp)
  site_cov <- bind_cols(cov,site_metastatic=as.numeric(d$sample_class=='metastatic'))
  site_res <- partial_multi(d$psps,d$lipid,site_cov,'major5_cell_proportions_plus_site',sm,'PSPS_core',lp)
  bind_rows(base_res,site_res)
}))))
sc_partial$q_BH <- p.adjust(sc_partial$p_value,'BH')
write.csv(sc_partial,file.path(out,'sc_weighted_pseudobulk_partial_cell_proportions.csv'),row.names=FALSE)

# Composition-only simulation: random macrophage fractions in [0,0.40],
# preserving each sample's non-macrophage relative composition.
set.seed(292815)
base_prop <- prop_tbl %>% select(sample_id,cell_type,prop) %>% pivot_wider(names_from=cell_type,values_from=prop,values_fill=0)
base_prop$macrophage <- if('macrophage' %in% names(base_prop)) base_prop$macrophage else 0
score_wide <- pb_major %>% select(sample_id,cell_type,starts_with('UCell_')) %>% pivot_wider(names_from=cell_type,values_from=starts_with('UCell_'))
# Use UCell PSPS and the four lipid scores for the simulation.
sim_mat <- pb_major %>% select(sample_id,cell_type,all_of(paste0('UCell_',set_order))) %>% pivot_wider(names_from=cell_type,values_from=all_of(paste0('UCell_',set_order)))
sim_out <- bind_rows(lapply(seq_len(500),function(i) {
  mp <- runif(nrow(base_prop),0,0.40)
  res <- base_prop; old <- res$macrophage; non <- setdiff(names(res),c('sample_id','macrophage'))
  scale <- ifelse(1-old>0,(1-mp)/(1-old),0)
  res[non] <- res[non]*scale; res$macrophage <- mp
  sc <- sim_mat %>% left_join(res,by='sample_id')
  calc <- function(nm) rowSums(as.matrix(sc[,paste0(nm,'.x')]),na.rm=TRUE)
  ps <- numeric(); li <- numeric()
  for(j in c('UCell_PSPS_core','UCell_Hallmark_FA_Metabolism','UCell_Hallmark_Adipogenesis','UCell_Reactome_FA_Metabolism','UCell_PPAR_Pathway')) {
    if(!(j %in% names(sim_mat))) next
  }
  # score columns in sim_mat are one per cell type after pivoting with list names.
  ps <- rep(0,nrow(res)); fa <- rep(0,nrow(res))
  for(ct in decomp_types) {
    psct <- sc[[paste0('UCell_PSPS_core_',ct)]]; fact <- sc[[paste0('UCell_Hallmark_FA_Metabolism_',ct)]]
    psct[!is.finite(psct)] <- 0; fact[!is.finite(fact)] <- 0
    ps <- ps + res[[ct]]*psct
    fa <- fa + res[[ct]]*fact
  }
  c <- suppressWarnings(cor.test(ps,fa,method='spearman',exact=FALSE))
  tibble(iteration=i,rho=unname(c$estimate),p_value=c$p.value)
}))
write.csv(sim_out,file.path(out,'sc_composition_only_simulation_UCell.csv'),row.names=FALSE)
sim_summary <- tibble(metric=c('rho_median','rho_q025','rho_q975','fraction_rho_le_minus_0.5','fraction_rho_le_minus_0.7'),value=c(median(sim_out$rho),quantile(sim_out$rho,.025),quantile(sim_out$rho,.975),mean(sim_out$rho<=-0.5),mean(sim_out$rho<=-0.7)))
write.csv(sim_summary,file.path(out,'sc_composition_only_simulation_summary.csv'),row.names=FALSE)

# Compact report.
bulk_partial_top <- bulk_partial %>% filter(psps_score=='PSPS_core') %>% arrange(panel,score_method,p_value) %>% slice_head(n=12)
bulk_report <- bulk_partial %>% filter(panel %in% c('xCell_core','EPIC_core'),score_method %in% c('UCell','AddModuleScore')) %>% arrange(panel,score_method,p_value) %>% mutate(txt=paste0(panel,'/',score_method,'/',lipid_score,': rho=',sprintf('%.3f',rho),', P=',format.pval(p_value,digits=3,eps=1e-3))) %>% pull(txt)
sc_report <- sc_partial %>% filter(panel=='major5_cell_proportions_plus_site') %>% mutate(txt=paste0(score_method,'/',lipid_score,': rho=',sprintf('%.3f',rho),', P=',format.pval(p_value,digits=3,eps=1e-3))) %>% pull(txt)
decomp_report <- decomp_results %>% filter(score_method=='UCell',scope=='pooled',mixing %in% c('W','U')) %>% group_by(lipid_score) %>% summarise(rw=rho[mixing=='W'][1],ru=rho[mixing=='U'][1],.groups='drop') %>% mutate(txt=paste0(lipid_score,': W rho=',sprintf('%.3f',rw),', U rho=',sprintf('%.3f',ru))) %>% pull(txt)
sim_report <- sim_summary %>% transmute(txt=paste0(metric,'=',sprintf('%.3f',value))) %>% pull(txt)
writeLines(c(
  '# bulk + 单细胞组成分解与脂质程序整合',
  '',
  paste0('小鼠 bulk：',ncol(expr),' 个样本；GSE154778：',ncol(sc_counts),' 个细胞、',n_distinct(sc_meta$sample_id),' 个样本。'),
  '基因集：HALLMARK_FATTY_ACID_METABOLISM、HALLMARK_ADIPOGENESIS、REACTOME_FATTY_ACID_METABOLISM、KEGG_PPAR_SIGNALING_PATHWAY，以及 PSPS_core 六基因。',
  paste0('同一 MSigDB 版本：',unique(gene_prov$db_version[!is.na(gene_prov$db_version)])[1], '；小鼠 bulk 通过 msigdbr 的 human MSigDB -> Mus musculus ortholog mapping，GSE154778 使用 human symbols。'),
  '评分方法：主分析 UCell；并行报告 Seurat AddModuleScore 兼容实现；AUCell 作为单细胞敏感性分析。人鼠分数只在各自物种内部比较符号与相对排序，不比较绝对数值。',
  '细胞过滤：β 细胞 7 个不纳入主要统计；α 细胞 11 个同样不作稳定来源判断；主分析使用 macrophage、acinar、stellate、ductal、epithelial_tumor。',
  '',
  'Bulk 去卷积说明：xCell 为富集分数、EPIC 为约束细胞比例；两者均先用 babelgene 将小鼠基因映射到人类同源符号，因此是跨物种 proxy 协变量。',
  '单细胞主推：sample×cell_type 原始 counts pseudobulk；细胞级 mixed model 仅作支持性分析，避免把 14,926 个细胞当独立生物学重复。',
  '组成分解：W=实测组成加权，U=细胞类型等权平均，F=全样本平均组成固定加权；δ=|rho_W|-|rho_U|。偏相关使用五个主要细胞类型比例，并另行加入 site（10 个 primary、6 个 metastatic）联合校正。',
  '组成模拟：将每个样本巨噬细胞比例随机扰动到 0–40%，保持其余细胞类型相对比例，用于评估单纯组成变化能否产生 bulk-like 负相关。',
  '',
  '关键数值（完整结果见 CSV）：',
  paste0('bulk 偏相关：',paste(bulk_report,collapse='；')),
  paste0('单细胞比例校正偏相关：',paste(sc_report,collapse='；')),
  paste0('UCell 组成分解（pooled）：',paste(decomp_report,collapse='；')),
  paste0('组成模拟：',paste(sim_report,collapse='；')),
  '',
  '解释边界：这些结果检验“焦亡评分与脂质程序的关系是否可由细胞组成稀释解释”，不能证明脂肪胰、细胞间因果或功能性焦亡。'
),file.path(out,'report_zh.md'),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('Completed bulk + single-cell composition/lipid integration: ',out)
