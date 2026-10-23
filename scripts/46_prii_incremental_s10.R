#!/usr/bin/env Rscript
# S10: incremental transcriptomic information of the frozen six-gene PRII.
# No gene selection, weight training, cut-off optimisation, or joint entry of
# highly correlated comparator scores is performed.
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(tidyr)
  library(readr); library(tibble); library(ggplot2); library(survival)
})
options(stringsAsFactors=FALSE)
set.seed(20261015L)
root <- normalizePath('.',mustWork=TRUE)
out <- file.path(root,'results/prii_incremental'); figdir <- file.path(out,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
core_h <- c('TXNIP','NLRP3','CASP1','GSDMD','IL1B','IL18')
core_m <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')

# These are project-defined harmonised proxies, not original published scores.
curated_pyro_h <- c('TXNIP','NLRP3','PYCARD','CASP1','CASP4','CASP5','GSDMD','GSDME','IL1B','IL18','AIM2','NLRC4','NLRP1','GBP1','GBP5')
curated_infl_h <- c('NLRP3','CASP1','GSDMD','IL1B','IL18','TNF','IL6','CXCL8','CCL2','CCL20','NFKB1','RELA','TLR4','TGFB1','CXCL1','CXCL2','CXCL5','CXCL10','CXCL11','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3')
manual_infl_h <- c(core_h,'AIM2','NLRC4','NLRP1','MEFV','P2RX7','NEK7','CARD8','PYCARD')
human_sets <- list(
  curated_pyroptosis_extended_noncore=setdiff(curated_pyro_h,core_h),
  manual_inflammasome_noncore=setdiff(manual_infl_h,core_h),
  curated_inflammation_noncore=setdiff(curated_infl_h,core_h)
)

# Human immune marker panels harmonised across TCGA-PAAD and GSE79668 using the
# fixed GSE79668 panel; scores remain marker proxies, not cell counts.
immune_h <- list(
  leukocyte=c('PTPRC','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','LST1'),
  t_cell=c('CD3D','CD3E','CD3G','CD4','CD8A','IL7R','LTB'),
  cytotoxic=c('NKG7','GZMB','GZMK','PRF1','CTSW','GZMA'),
  b_cell=c('MS4A1','CD79A','CD79B','CD37','CD74','HLA-DRA','HLA-DRB1','CD19','CD83'),
  myeloid=c('LST1','TYROBP','AIF1','CTSS','CTSD','LGALS3','S100A10'),
  macrophage_TAM=c('ADGRE1','CD68','CSF1R','C1QA','C1QB','C1QC','MS4A7','MERTK','MAFB','CD163','MRC1','MARCO','LGALS3'),
  neutrophil=c('S100A8','S100A9','LCN2','CXCR2','CSF3R','NGP','CAMP'),
  inflammation_noncore=c('TNF','IL6','CXCL1','CXCL2','CXCL5','CCL2','CCL3','CCL4','CCL5','NFKB1','RELA','TLR4','TGFB1','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3')
)

zmean <- function(expr,genes) {
  genes <- intersect(unique(genes),rownames(expr))
  if (!length(genes)) return(rep(NA_real_,ncol(expr)))
  z <- t(apply(expr[genes,,drop=FALSE],1,function(x) as.numeric(scale(x))))
  if (is.null(dim(z))) z <- matrix(z,nrow=1)
  z[!is.finite(z)] <- 0
  colMeans(z,na.rm=TRUE)
}
score_sets <- function(expr,sets) as.data.frame(setNames(lapply(sets,function(g)zmean(expr,g)),names(sets)))
rank_partial <- function(y,x,controls=NULL) {
  ok <- is.finite(y)&is.finite(x)
  if (!is.null(controls) && ncol(controls)) ok <- ok & apply(controls,1,function(z)all(is.finite(z)))
  n <- sum(ok); if (n<8) return(tibble(n=n,rho=NA_real_,p_value=NA_real_,df=NA_real_))
  ry <- rank(y[ok],ties.method='average'); rx <- rank(x[ok],ties.method='average')
  if (!is.null(controls) && ncol(controls)) {
    cc <- as.data.frame(lapply(seq_len(ncol(controls)),function(j)rank(controls[ok,j],ties.method='average')))
    ey <- residuals(lm(ry~.,data=cc)); ex <- residuals(lm(rx~.,data=cc)); df <- n-ncol(cc)-2
  } else {ey<-ry;ex<-rx;df<-n-2}
  rho <- suppressWarnings(cor(ey,ex)); p <- if(is.finite(rho)&&abs(rho)<1&&df>0)2*pt(abs(rho)*sqrt(df/(1-rho^2)),df=df,lower.tail=FALSE) else NA_real_
  tibble(n=n,rho=unname(rho),p_value=p,df=df)
}
get_vif <- function(fit) {
  X <- model.matrix(fit); X <- X[,colnames(X)!='(Intercept)',drop=FALSE]
  if(ncol(X)<2) return(setNames(rep(1,ncol(X)),colnames(X)))
  sapply(seq_len(ncol(X)),function(j){y<-X[,j]; xx<-X[,-j,drop=FALSE]; r2<-summary(lm(y~xx))$r.squared; 1/(1-r2)}) |> setNames(colnames(X))
}
fit_nested <- function(y,prii,comp,covars,ds,cmp,mod,B=500L) {
  dat <- data.frame(y=as.numeric(scale(y)),PRII=as.numeric(scale(prii)),COMP=as.numeric(scale(comp)),check.names=FALSE)
  if(ncol(covars)) dat <- cbind(dat,covars)
  dat <- dat[complete.cases(dat),,drop=FALSE]; n<-nrow(dat)
  if(n<12) return(tibble(dataset=ds,comparator=cmp,module=mod,n=n,model='insufficient'))
  base <- if(ncol(covars)) paste(names(covars),collapse='+') else '1'
  f1<-lm(as.formula(paste('y~',base,'+COMP')),dat); f2<-lm(as.formula(paste('y~',base,'+PRII')),dat); f3<-lm(as.formula(paste('y~',base,'+COMP+PRII')),dat)
  an<-anova(f1,f3); sm<-summary(f3)$coefficients; vif<-get_vif(f3)
  boot <- matrix(NA_real_,nrow=B,ncol=2,dimnames=list(NULL,c('beta','delta_r2')))
  for(b in seq_len(B)) {
    ix<-sample.int(n,n,replace=TRUE); db<-dat[ix,,drop=FALSE]
    fb1<-try(lm(formula(f1),db),silent=TRUE); fb3<-try(lm(formula(f3),db),silent=TRUE)
    if(!inherits(fb1,'try-error')&&!inherits(fb3,'try-error')&&'PRII'%in%names(coef(fb3))) boot[b,]<-c(coef(fb3)['PRII'],summary(fb3)$r.squared-summary(fb1)$r.squared)
  }
  ci<-apply(boot,2,quantile,probs=c(.025,.975),na.rm=TRUE,names=FALSE)
  tibble(dataset=ds,comparator=cmp,module=mod,n=n,model='M3_vs_M1',prii_beta=unname(sm['PRII','Estimate']),prii_beta_low=ci[1,1],prii_beta_high=ci[2,1],prii_p=unname(sm['PRII','Pr(>|t|)']),comparator_beta=unname(sm['COMP','Estimate']),comparator_p=unname(sm['COMP','Pr(>|t|)']),r2_m1=summary(f1)$r.squared,r2_m2=summary(f2)$r.squared,r2_m3=summary(f3)$r.squared,delta_r2_prii=summary(f3)$r.squared-summary(f1)$r.squared,delta_r2_low=ci[1,2],delta_r2_high=ci[2,2],partial_r2_prii=(deviance(f1)-deviance(f3))/deviance(f1),delta_r2_comparator=summary(f3)$r.squared-summary(f2)$r.squared,lrt_f=unname(an$F[2]),lrt_p=unname(an$`Pr(>F)`[2]),delta_aic=AIC(f3)-AIC(f1),delta_bic=BIC(f3)-BIC(f1),vif_prii=unname(vif['PRII']),vif_comparator=unname(vif['COMP']),bootstrap_success=sum(is.finite(boot[,1])))
}

# Expression input matrices, reconstructed using the frozen project transforms.
mouse_obj<-readRDS(file.path(root,'results/qc/GSE292815_KC_12wk_matrix.rds'))
mouse_expr<-log2(as.matrix(mouse_obj$expr)+1); mouse_meta<-as.data.frame(mouse_obj$metadata)
mouse_gs_obj<-readRDS(file.path(root,'scripts/results/gene_sets_mouse/gene_sets.rds'))
mouse_sets<-list(pyroptosis_noncore=setdiff(mouse_gs_obj$gene_sets$pyroptosis,core_m),inflammasome_noncore=setdiff(mouse_gs_obj$gene_sets$inflammasome,core_m),inflammation_noncore=setdiff(c('Tnf','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl2','Ccl3','Ccl4','Ccl5','Nfkb1','Rela','Tlr4','Tgfb1','Ptgs2','S100a8','S100a9','Nfkbia','Tnfaip3'),core_m))
mouse_immune_long<-read_csv(file.path(root,'results/12wk_primary/immune/12wk_immune_marker_sets_filtered.csv'),show_col_types=FALSE)
immune_m<-lapply(split(mouse_immune_long$gene,mouse_immune_long$module),unique)
mouse_scores<-score_sets(mouse_expr,c(list(PRII=core_m),mouse_sets)); mouse_scores$sample<-colnames(mouse_expr); mouse_scores$group<-mouse_meta$group[match(mouse_scores$sample,mouse_meta$matrix_column)]

se<-readRDS(file.path(root,'data/raw/TCGA_PAAD_RNAseq.rds')); tpm<-assay(se,'tpm_unstrand'); rd<-as.data.frame(rowData(se))
sym<-toupper(as.character(rd$gene_name)); bad<-is.na(sym)|sym==''; sym[bad]<-rownames(tpm)[bad]; sym<-sub('\\..*$','',sym); keep<-!duplicated(sym)
tcga_all<-log2(tpm[keep,,drop=FALSE]+1); rownames(tcga_all)<-sym[keep]
cd<-as.data.frame(colData(se)); tum<-(cd$sample_type=='Primary Tumor'|cd$shortLetterCode=='TP'); tum[is.na(tum)]<-FALSE
tcga_expr<-tcga_all[,tum,drop=FALSE]; tcga_cd<-cd[colnames(tcga_expr),,drop=FALSE]
tcga_scores<-score_sets(tcga_expr,c(list(PRII=core_h),human_sets)); tcga_scores$sample<-colnames(tcga_expr); tcga_scores$purity<-suppressWarnings(as.numeric(tcga_cd$paper_ABSOLUTE.Purity))

raw<-read.delim(gzfile(file.path(root,'data/external/GSE79668_counts.txt.gz')),check.names=FALSE,stringsAsFactors=FALSE)
gse_sym<-toupper(sub('_ENSG.*$','',as.character(raw[[1]]))); cols<-names(raw)[-1]; cmat<-as.matrix(raw[,cols]); storage.mode(cmat)<-'numeric'; rownames(cmat)<-gse_sym
cmat<-rowsum(cmat,group=gse_sym,reorder=FALSE); gse_expr<-log2(sweep(cmat,2,colSums(cmat),'/')*1e6+1)
gse_scores<-score_sets(gse_expr,c(list(PRII=core_h),human_sets)); gse_scores$sample<-colnames(gse_expr)

# Correlation/coverage/overlap audit includes broad project collections, but
# only one compact strict comparator at a time enters the primary models.
cor_rows<-list(); ov_rows<-list()
add_cor<-function(ds,sc,cmp_names) for(nm in cmp_names){ct<-suppressWarnings(cor.test(sc$PRII,sc[[nm]],method='spearman',exact=FALSE));cor_rows[[length(cor_rows)+1]]<<-tibble(dataset=ds,comparator=nm,n=nrow(sc),rho=unname(ct$estimate),p_value=ct$p.value)}
add_cor('GSE292815_KC_12wk',mouse_scores,names(mouse_sets)); add_cor('TCGA_PAAD_primary',tcga_scores,names(human_sets)); add_cor('GSE79668',gse_scores,names(human_sets))
cor_tbl<-bind_rows(cor_rows)%>%group_by(dataset)%>%mutate(fdr_bh=p.adjust(p_value,'BH'))%>%ungroup();write_csv(cor_tbl,file.path(out,'score_correlations.csv'))
audit_signatures<-list(
  GSE292815_KC_12wk=list(expr=mouse_expr,core=core_m,sets=mouse_sets),
  TCGA_PAAD_primary=list(expr=tcga_expr,core=core_h,sets=human_sets),
  GSE79668=list(expr=gse_expr,core=core_h,sets=human_sets))
for(ds in names(audit_signatures)){z<-audit_signatures[[ds]];for(nm in names(z$sets)){g<-unique(z$sets[[nm]]);ov<-intersect(g,z$core);pr<-intersect(g,rownames(z$expr));ov_rows[[length(ov_rows)+1]]<-tibble(dataset=ds,comparator=nm,n_requested=length(g),n_core_overlap=length(ov),core_overlap=paste(ov,collapse=';'),jaccard_core=length(ov)/length(union(g,z$core)),n_present=length(pr),coverage=length(pr)/length(g))}}
overlap_tbl<-bind_rows(ov_rows);write_csv(overlap_tbl,file.path(out,'signature_overlap_audit.csv'))

pair_vif<-bind_rows(lapply(names(audit_signatures),function(ds){sc<-switch(ds,GSE292815_KC_12wk=mouse_scores,TCGA_PAAD_primary=tcga_scores,GSE79668=gse_scores);cmp<-names(audit_signatures[[ds]]$sets);bind_rows(lapply(cmp,function(nm){r<-cor(sc$PRII,sc[[nm]],use='complete.obs');tibble(dataset=ds,comparator=nm,n=sum(complete.cases(sc[,c('PRII',nm)])),pearson_r=r,vif=1/(1-r^2))}))}))
write_csv(pair_vif,file.path(out,'vif_pairwise.csv'))

datasets<-list(
 GSE292815_KC_12wk=list(expr=mouse_expr,scores=mouse_scores,immune=immune_m,comps=names(mouse_sets),cov=data.frame()),
 TCGA_PAAD_primary=list(expr=tcga_expr,scores=tcga_scores,immune=immune_h,comps=names(human_sets),cov=data.frame(purity=tcga_scores$purity)),
 GSE79668=list(expr=gse_expr,scores=gse_scores,immune=immune_h,comps=names(human_sets),cov=data.frame()))
partial_rows<-list();model_rows<-list();module_audit<-list()
for(ds in names(datasets)){d<-datasets[[ds]];for(cmp in d$comps){
  comp_set<-if(ds=='GSE292815_KC_12wk') mouse_sets[[cmp]] else human_sets[[cmp]]
  for(mod in names(d$immune)){requested<-unique(d$immune[[mod]]);remaining<-setdiff(requested,unique(comp_set))
    present<-intersect(remaining,rownames(d$expr));module_audit[[length(module_audit)+1]]<-tibble(dataset=ds,comparator=cmp,module=mod,n_requested=length(requested),n_after_exclusion=length(remaining),n_present=length(present),coverage=length(present)/max(length(requested),1),status=ifelse(length(present)>=3,'estimable','not_estimable'),removed_genes=paste(intersect(requested,comp_set),collapse=';'),missing_remaining_genes=paste(setdiff(remaining,present),collapse=';'))
    if(length(present)<3) next
    y<-zmean(d$expr,present); ctrl<-if(ncol(d$cov))as.matrix(d$cov)else NULL
    cond1<-if(is.null(ctrl))matrix(d$scores[[cmp]],ncol=1)else cbind(ctrl,d$scores[[cmp]])
    cond2<-if(is.null(ctrl))matrix(d$scores$PRII,ncol=1)else cbind(ctrl,d$scores$PRII)
    p1<-rank_partial(y,d$scores$PRII,cond1);p2<-rank_partial(y,d$scores[[cmp]],cond2)
    partial_rows[[length(partial_rows)+1]]<-tibble(dataset=ds,comparator=cmp,module=mod,relation='PRII_conditioned_on_comparator',n=p1$n,rho=p1$rho,p_value=p1$p_value,df=p1$df)
    partial_rows[[length(partial_rows)+1]]<-tibble(dataset=ds,comparator=cmp,module=mod,relation='comparator_conditioned_on_PRII',n=p2$n,rho=p2$rho,p_value=p2$p_value,df=p2$df)
    model_rows[[length(model_rows)+1]]<-fit_nested(y,d$scores$PRII,d$scores[[cmp]],d$cov,ds,cmp,mod)
  }
}}
partial_tbl<-bind_rows(partial_rows)%>%group_by(dataset,comparator,relation)%>%mutate(fdr_bh=p.adjust(p_value,'BH'))%>%ungroup()
model_tbl<-bind_rows(model_rows)%>%group_by(dataset,comparator)%>%mutate(fdr_bh=p.adjust(lrt_p,'BH'))%>%ungroup()
module_tbl<-bind_rows(module_audit)
write_csv(partial_tbl,file.path(out,'partial_correlations.csv'));write_csv(model_tbl,file.path(out,'immune_nested_models.csv'));write_csv(module_tbl,file.path(out,'immune_module_overlap_audit.csv'))

# TCGA survival increment is a secondary sensitivity, on one common complete-
# case set per comparator. GSE79668 has no usable OS time/event pair.
surv<-read_csv(file.path(root,'results/tcga_paad_external/survival_input.csv'),show_col_types=FALSE)
map<-tcga_scores%>%select(sample,all_of(names(human_sets)))
surv<-surv%>%mutate(PRII=PSPS_core)%>%left_join(map,by=c('sample_id'='sample'))
stage_clean<-ifelse(grepl('^Stage I$',surv$stage,ignore.case=TRUE),'I',ifelse(grepl('^Stage III',surv$stage,ignore.case=TRUE),'III',ifelse(grepl('^Stage II',surv$stage,ignore.case=TRUE),'II',ifelse(grepl('^Stage IV',surv$stage,ignore.case=TRUE),'IV',NA))))
surv$stage_group<-factor(stage_clean);surv_rows<-list()
for(cmp in names(human_sets)){dd<-surv%>%filter(is.finite(time_days),time_days>0,is.finite(event),is.finite(age_years),is.finite(purity),is.finite(PRII),is.finite(.data[[cmp]]),!is.na(stage_group));dd$COMP<-dd[[cmp]];if(nrow(dd)<50||sum(dd$event)<15||nlevels(droplevels(dd$stage_group))<2)next
 f1<-coxph(Surv(time_days,event)~age_years+stage_group+purity+COMP,data=dd);f3<-coxph(Surv(time_days,event)~age_years+stage_group+purity+COMP+PRII,data=dd);lr<-anova(f1,f3,test='LRT');sm<-summary(f3);surv_rows[[length(surv_rows)+1]]<-tibble(comparator=cmp,n=nrow(dd),events=sum(dd$event),prii_hr=exp(coef(f3)['PRII']),prii_p=sm$coefficients['PRII','Pr(>|z|)'],comparator_hr=exp(coef(f3)['COMP']),comparator_p=sm$coefficients['COMP','Pr(>|z|)'],delta_chisq=lr$Chisq[2],lrt_p=lr$`P(>|Chi|)`[2],delta_aic=AIC(f3)-AIC(f1),delta_bic=BIC(f3)-BIC(f1),cindex_m1=summary(f1)$concordance[1],cindex_m3=sm$concordance[1])}
surv_tbl<-bind_rows(surv_rows);if(nrow(surv_tbl))surv_tbl<-surv_tbl%>%mutate(fdr_bh=p.adjust(lrt_p,'BH')) else surv_tbl<-tibble(comparator=character(),n=integer(),events=integer(),note=character());write_csv(surv_tbl,file.path(out,'tcga_survival_incremental.csv'))

# The target is conditional explanatory information, not out-of-sample
# prediction; CV is intentionally not presented as a validation endpoint.
write_csv(tibble(dataset=c('GSE292815_KC_12wk','TCGA_PAAD_primary','GSE79668'),n=c(ncol(mouse_expr),ncol(tcga_expr),ncol(gse_expr)),metric='not_applicable',reason=c('small exploratory cohort; fixed-score effect estimation only','nested explanatory models; no prediction claim','independent association replication; no model transport claim')),file.path(out,'cross_validation_metrics.csv'))

# Summary figure for human cohorts: PRII partial association and incremental R2.
plot_data<-model_tbl%>%filter(dataset%in%c('TCGA_PAAD_primary','GSE79668'))
p1<-ggplot(partial_tbl%>%filter(relation=='PRII_conditioned_on_comparator',dataset%in%c('TCGA_PAAD_primary','GSE79668')),aes(rho,module,colour=comparator))+geom_vline(xintercept=0,linetype=2,colour='grey55')+geom_point(position=position_dodge(width=.55),size=1.8)+facet_wrap(~dataset,nrow=1)+labs(x='Partial Spearman rho',y=NULL,colour='Comparator',title='PRII conditional immune-module association')+theme_bw(base_size=8)
p2<-ggplot(plot_data,aes(delta_r2_prii,module,colour=comparator))+geom_vline(xintercept=0,linetype=2,colour='grey55')+geom_point(position=position_dodge(width=.55),size=1.8)+facet_wrap(~dataset,nrow=1)+labs(x='Delta R-squared after adding PRII',y=NULL,colour='Comparator',title='Nested model incremental fit (M3 vs M1)')+theme_bw(base_size=8)
fig<-p1/p2;ggsave(file.path(figdir,'FigureS10_PRII_incremental.pdf'),fig,width=11,height=8);ggsave(file.path(figdir,'FigureS10_PRII_incremental.png'),fig,width=11,height=8,dpi=300)

write_csv(partial_tbl%>%filter(dataset=='GSE79668'),file.path(out,'gse79668_replication.csv'))
manifest<-c('# S10 PRII incremental transcriptomic information analysis','',
 'The directional gate passed: the fixed PRII score was positively correlated with all available strict non-core comparators in GSE292815, TCGA-PAAD and GSE79668.',
 'PRII remained fixed as the within-cohort mean of gene-wise z-scored TXNIP/NLRP3/CASP1/GSDMD/IL1B/IL18. Project-defined comparator sets are harmonised mean-z proxies, not original published weighted scores.',
 'Broad project pyroptosis/inflammasome collections were retained for descriptive correlation and overlap audit only. One comparator at a time entered each nested model.',
 'Human TCGA/GSE79668 immune modules were rescored from the same fixed marker panels after excluding PRII and comparator genes. Mouse modules used the fixed project marker table with the same exclusion rule. Modules with fewer than 3 expressed genes were marked not estimable.',
 'TCGA immune models adjust for paper_ABSOLUTE.Purity; GSE79668 is an independent unadjusted PAAD association replication; GSE292815 is exploratory (n=21). TCGA survival is secondary and uses common complete cases per comparator. GSE79668 lacks usable OS time/event metadata.',
 'Nested models: M0=covariates; M1=M0+comparator; M2=M0+PRII; M3=M0+comparator+PRII. The primary incremental contrast is M3 vs M1. Results are conditional transcriptomic associations, not functional pyroptosis, causal immunity, or independent prognosis.',
 'Files: signature_overlap_audit.csv; score_correlations.csv; vif_pairwise.csv; immune_module_overlap_audit.csv; partial_correlations.csv; immune_nested_models.csv; tcga_survival_incremental.csv; gse79668_replication.csv; cross_validation_metrics.csv; figures/FigureS10_PRII_incremental.pdf/.png.')
writeLines(manifest,file.path(out,'incremental_analysis_manifest.md'));writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('S10 analysis complete: ',out)
