#!/usr/bin/env Rscript
# Reproducible S10 comparison of PRII with project-defined and published-signature proxies.
# The Zuo comparator is a within-cohort standardized directional proxy, not an exact
# reproduction of the published risk score. Xie et al. is descriptive only.
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(dplyr); library(readr); library(tidyr)
  library(tibble); library(ggplot2); library(patchwork)
})
options(stringsAsFactors=FALSE)
root <- normalizePath('.', mustWork=TRUE)
out <- file.path(root,'results/prii_incremental')
sigout <- file.path(out,'signatures')
auditout <- file.path(out,'s10_signature_audit')
figout <- file.path(auditout,'figures')
dir.create(sigout,recursive=TRUE,showWarnings=FALSE); dir.create(figout,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,'results/prii_incremental/s10_signature_audit.R'))
core_h <- c('TXNIP','NLRP3','CASP1','GSDMD','IL1B','IL18')
zuo_beta <- c(ACTA2=-0.135041,C1QTNF9=-0.063701,DNAH8=-0.039539,GATM=-0.097678,LBP=-0.028119,NGF=-0.018908)
zuo_genes <- names(zuo_beta)
xie_genes <- c('ARSG','SLC26A11','MET')
xie52 <- read_csv(file.path(sigout,'Xie2023_candidate_PRG_52.csv'),show_col_types=FALSE)$gene

# Score genes using within-cohort z scores, matching the frozen PRII implementation.
score_z_expr <- function(mat, genes) {
  genes <- intersect(unique(genes),rownames(mat))
  if (!length(genes)) return(rep(NA_real_,ncol(mat)))
  z <- t(scale(t(mat[genes,,drop=FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}
score_weighted_z <- function(mat,beta) {
  genes <- intersect(names(beta),rownames(mat))
  if (!length(genes)) return(rep(NA_real_,ncol(mat)))
  z <- t(scale(t(mat[genes,,drop=FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(crossprod(beta[genes],z))
}
prii_tcga <- score_z_expr(expr,core_h)
zuo_tcga <- score_weighted_z(expr,zuo_beta)
xie_tcga <- score_z_expr(expr,xie_genes)
prii_gse <- score_z_expr(gexpr,core_h)
zuo_gse <- score_weighted_z(gexpr,zuo_beta)
xie_gse <- score_z_expr(gexpr,xie_genes)

# Recreate prespecified immune modules after removing every gene in PRII and the
# comparator being tested. This prevents arithmetic overlap from inflating associations.
tcga_sets <- list(
  immune_infiltration=c('PTPRC','HLA-DRA','HLA-DRB1','CD3D','CD3E','CD8A','MS4A1','NKG7','LST1','TYROBP','FCER1G','CD68'),
  t_cell=c('CD3D','CD3E','CD3G','TRBC1','TRBC2','CD4','CD8A','CD8B','IL7R','LTB'),
  cytotoxic=c('NKG7','GNLY','GZMB','GZMH','PRF1','GZMK','CTSW'),
  b_cell=c('MS4A1','CD79A','CD37','CD74','HLA-DRA','CD79B','CD19'),
  myeloid=c('LST1','TYROBP','FCER1G','AIF1','CTSS','LILRB1','CTSD','LGALS3','S100A10'),
  macrophage_TAM=c('CD68','LST1','TYROBP','FCER1G','AIF1','C1QA','C1QB','C1QC','MS4A7','CTSS','CD163','MRC1','CSF1R','MARCO','VSIG4','LGALS3'),
  CAF=c('COL1A1','COL1A2','COL3A1','COL5A1','COL6A1','COL6A2','DCN','LUM','SPARC','FAP','PDGFRA','PDGFRB','ACTA2','TAGLN','POSTN','COL11A1'),
  EMT=c('VIM','CDH2','SNAI1','SNAI2','TWIST1','TWIST2','ZEB1','ZEB2','MMP2','MMP9','MMP14','ITGA5','ITGB1','TGFB1','TGFBR1','TGFBR2','SERPINE1','FN1','COL1A1','COL3A1'),
  lipid_metabolism=c('ACACA','ACACB','FASN','SCD','ELOVL6','DGAT1','DGAT2','PLIN2','PLIN3','FABP4','FABP5','CPT1A','ACOX1','CD36','LPL','LDLR','SREBF1','SREBF2','HMGCR','HMGCS1','FDFT1','SQLE','ABCA1','APOE','PPARG','PPARA','CYP27A1','LIPE','PNPLA2'),
  inflammation_noncore=c('TNF','IL6','CXCL8','CCL2','CCL20','NFKB1','RELA','TLR4','TGFB1','CXCL1','CXCL2','CXCL5','CXCL10','CXCL11','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3'))
gse_sets <- list(
  leukocyte=c('PTPRC','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','LST1'),
  t_cell=c('CD3D','CD3E','CD3G','CD4','CD8A','IL7R','LTB'),
  cytotoxic=c('NKG7','GZMB','GZMK','PRF1','CTSW','GZMA'),
  b_cell=c('MS4A1','CD79A','CD79B','CD37','CD74','HLA-DRA','HLA-DRB1','CD19','CD83'),
  myeloid=c('LST1','TYROBP','AIF1','CTSS','CTSD','LGALS3','S100A10'),
  macrophage_TAM=c('ADGRE1','CD68','CSF1R','C1QA','C1QB','C1QC','MS4A7','MERTK','MAFB','CD163','MRC1','MARCO','LGALS3'),
  neutrophil=c('S100A8','S100A9','LCN2','CXCR2','CSF3R','NGP','CAMP'),
  inflammation_noncore=c('TNF','IL6','CXCL1','CXCL2','CXCL5','CCL2','CCL3','CCL4','CCL5','NFKB1','RELA','TLR4','TGFB1','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3'))

# Predictor sets and cohort-level proxy scores.
project_sets <- list(
  curated_pyroptosis_extended_noncore=setdiff(curated_pyro_h,core_h),
  manual_inflammasome_noncore=setdiff(manual_infl_h,core_h),
  curated_inflammation_noncore=curated_infl_noncore_h,
  project_pyroptosis_noncore=setdiff(gs_h$pyroptosis,core_h),
  project_inflammasome_noncore=setdiff(gs_h$inflammasome,core_h))
comp_sets <- c(project_sets,list(Zuo2022_standardized_directional_proxy=zuo_genes,
                                 Xie2023_descriptive_3gene_panel=xie_genes))
comp_tcga <- as.data.frame(setNames(lapply(project_sets,function(gs)score_z_expr(expr,gs)),names(project_sets)))
comp_tcga$Zuo2022_standardized_directional_proxy <- zuo_tcga
comp_tcga$Xie2023_descriptive_3gene_panel <- xie_tcga
comp_tcga$PRII <- prii_tcga
comp_gse <- as.data.frame(setNames(lapply(project_sets,function(gs)score_z_expr(gexpr,gs)),names(project_sets)))
comp_gse$Zuo2022_standardized_directional_proxy <- zuo_gse
comp_gse$Xie2023_descriptive_3gene_panel <- xie_gse
comp_gse$PRII <- prii_gse

# Correlation and VIF: one comparator at a time; never enter broad scores together.
dataset_scores <- list(TCGA_PAAD_primary=comp_tcga,GSE79668=comp_gse)
cor_tbl <- bind_rows(lapply(names(dataset_scores),function(ds) {
  s <- dataset_scores[[ds]]
  bind_rows(lapply(setdiff(names(s),'PRII'),function(cmp) {
    ok <- is.finite(s$PRII) & is.finite(s[[cmp]])
    ct <- suppressWarnings(cor.test(s$PRII[ok],s[[cmp]][ok],method='spearman',exact=FALSE))
    tibble(dataset=ds,comparator=cmp,n=sum(ok),spearman_rho=unname(ct$estimate),p_value=ct$p.value)
  }))
})) %>% group_by(dataset) %>% mutate(fdr_bh=p.adjust(p_value,'BH')) %>% ungroup()
write_csv(cor_tbl,file.path(sigout,'published_proxy_correlations.csv'))
write_csv(cor_tbl,file.path(auditout,'S10_published_proxy_correlations.csv'))

pub_vif <- bind_rows(lapply(names(dataset_scores),function(ds) {
  s <- dataset_scores[[ds]]
  bind_rows(lapply(setdiff(names(s),'PRII'),function(cmp) {
    ok <- complete.cases(s[,c('PRII',cmp)])
    r <- cor(s$PRII[ok],s[[cmp]][ok])
    tibble(dataset=ds,comparator=cmp,pearson_r=r,vif=1/(1-r^2),n=sum(ok))
  }))
}))
project_vif <- read_csv(file.path(auditout,'vif_prii_plus_single_comparator.csv'),show_col_types=FALSE)
write_csv(bind_rows(project_vif,pub_vif),file.path(auditout,'S10_vif_prii_single_signature.csv'))

# Gene overlap and expression coverage. Coverage is calculated for the non-core
# genes used in strict conditional models; Xie 52-gene candidate set is descriptive.
expr_by_ds <- list(TCGA_PAAD_primary=rownames(expr),GSE79668=rownames(gexpr))
gene_catalog <- bind_rows(
  bind_rows(lapply(names(project_sets),function(nm)tibble(comparator=nm,gene=project_sets[[nm]],evidence_class='project_defined'))),
  tibble(comparator='Zuo2022_standardized_directional_proxy',gene=zuo_genes,evidence_class='published_weighted_proxy'),
  tibble(comparator='Xie2023_descriptive_3gene_panel',gene=xie_genes,evidence_class='published_descriptive_only'),
  tibble(comparator='Xie2023_candidate_52gene_set',gene=unique(xie52),evidence_class='published_candidate_gene_set_only'))
write_csv(gene_catalog,file.path(auditout,'S10_comparator_gene_catalog.csv'))
full_sets <- list(
  curated_pyroptosis_extended_noncore=curated_pyro_h,
  manual_inflammasome_noncore=manual_infl_h,
  curated_inflammation_noncore=curated_infl_h,
  project_pyroptosis_noncore=gs_h$pyroptosis,
  project_inflammasome_noncore=gs_h$inflammasome,
  Zuo2022_standardized_directional_proxy=zuo_genes,
  Xie2023_descriptive_3gene_panel=xie_genes,
  Xie2023_candidate_52gene_set=unique(xie52))
overlap_table <- bind_rows(lapply(names(full_sets),function(nm) {
  genes <- unique(full_sets[[nm]]); core_overlap <- intersect(genes,core_h); strict <- setdiff(genes,core_h)
  bind_rows(lapply(names(expr_by_ds),function(ds) {
    expressed <- intersect(strict,expr_by_ds[[ds]])
    tibble(dataset=ds,comparator=nm,
      evidence_class=if(grepl('Zuo2022',nm))'published_weighted_proxy' else if(grepl('Xie2023',nm))'published_descriptive_only' else 'project_defined',
      n_signature_genes=length(genes),n_PRII_core_overlap=length(core_overlap),PRII_overlap_genes=paste(core_overlap,collapse=';'),
      n_strict_noncore_genes=length(strict),n_expressed_noncore=length(expressed),coverage_noncore=if(length(strict))length(expressed)/length(strict) else NA_real_,
      expressed_noncore_genes=paste(expressed,collapse=';'))
  }))
}))
write_csv(overlap_table,file.path(auditout,'S10_signature_gene_overlap_coverage.csv'))

# Nested OLS models. Comparisons are M1: purity + comparator, M2: M1 + PRII.
# `anova(m1,m2)` is the partial F test, not a likelihood-ratio test.
partial_rank <- function(d,y,x,covars=character()) {
  vars <- unique(c(y,x,covars)); v <- d[,vars,drop=FALSE]; v <- v[complete.cases(v),,drop=FALSE]
  if(nrow(v)<8) return(c(rho=NA_real_,p=NA_real_,n=nrow(v)))
  rx <- rank(v[[x]]); ry <- rank(v[[y]])
  if(length(covars)) {
    z <- as.data.frame(lapply(v[,covars,drop=FALSE],rank))
    rx <- residuals(lm(rx~.,data=z)); ry <- residuals(lm(ry~.,data=z))
  }
  ct <- suppressWarnings(cor.test(rx,ry,method='pearson'))
  c(rho=unname(ct$estimate),p=ct$p.value,n=nrow(v))
}
fit_one <- function(d,dataset,outcome,cmp,covars=character(),n_outcome_genes=NA_integer_,n_removed=NA_integer_) {
  vars <- c(outcome,cmp,'PRII',covars); v <- d[,vars,drop=FALSE]; v <- v[complete.cases(v),,drop=FALSE]
  if(nrow(v)<12) return(tibble())
  f0 <- as.formula(paste(outcome,'~',if(length(covars))paste(covars,collapse='+') else '1'))
  f_prii <- as.formula(paste(outcome,'~',paste(c(covars,'PRII'),collapse='+')))
  f1 <- as.formula(paste(outcome,'~',paste(c(covars,cmp),collapse='+')))
  f2 <- as.formula(paste(outcome,'~',paste(c(covars,cmp,'PRII'),collapse='+')))
  m0 <- lm(f0,v); m_prii <- lm(f_prii,v); m1 <- lm(f1,v); m2 <- lm(f2,v)
  an <- anova(m1,m2); an_cmp <- anova(m_prii,m2)
  z <- as.data.frame(lapply(v[,c(outcome,c(covars,cmp,'PRII')),drop=FALSE],function(x)as.numeric(scale(x))))
  names(z) <- c(outcome,c(covars,cmp,'PRII'))
  mz <- lm(f2,z); coef <- summary(mz)$coefficients
  X <- model.matrix(mz); e <- residuals(mz); h <- hatvalues(mz); bread <- solve(crossprod(X))
  meat <- crossprod(X * as.numeric(e/(1-h)))
  Vhc3 <- bread %*% meat %*% bread
  se_hc3 <- sqrt(diag(Vhc3)); df_hc3 <- df.residual(mz); tcrit <- qt(.975,df_hc3)
  beta_hc3 <- coef['PRII','Estimate']; se_beta_hc3 <- se_hc3[match('PRII',names(se_hc3))]
  p_hc3 <- 2*pt(-abs(beta_hc3/se_beta_hc3),df=df_hc3)
  pc <- partial_rank(v,outcome,'PRII',c(covars,cmp))
  cc <- partial_rank(v,outcome,cmp,c(covars,'PRII'))
  tibble(dataset=dataset,outcome=outcome,comparator=cmp,n=nrow(v),n_outcome_genes=n_outcome_genes,
    n_comparator_genes=length(intersect(comp_sets[[cmp]],if(dataset=='TCGA_PAAD')rownames(expr) else rownames(gexpr))),
    n_removed_from_outcome=n_removed,r2_M1=summary(m1)$r.squared,r2_M2=summary(m2)$r.squared,
    delta_r2_PRII=summary(m2)$r.squared-summary(m1)$r.squared,
    partial_r2_PRII=(deviance(m1)-deviance(m2))/deviance(m1),p_nested_F_PRII=an$`Pr(>F)`[2],
    delta_r2_comparator_after_PRII=summary(m2)$r.squared-summary(m_prii)$r.squared,
    partial_r2_comparator_after_PRII=(deviance(m_prii)-deviance(m2))/deviance(m_prii),p_nested_F_comparator_after_PRII=an_cmp$`Pr(>F)`[2],
    delta_AIC_PRII=AIC(m2)-AIC(m1),delta_BIC_PRII=BIC(m2)-BIC(m1),
    beta_PRII=coef['PRII','Estimate'],p_beta_PRII=coef['PRII','Pr(>|t|)'],beta_PRII_HC3=beta_hc3,beta_PRII_HC3_SE=se_beta_hc3,
    beta_PRII_HC3_CI_low=beta_hc3-tcrit*se_beta_hc3,beta_PRII_HC3_CI_high=beta_hc3+tcrit*se_beta_hc3,p_beta_PRII_HC3=p_hc3,
    beta_comparator=coef[cmp,'Estimate'],p_beta_comparator=coef[cmp,'Pr(>|t|)'],
    partial_rho_PRII=pc['rho'],partial_p_PRII=pc['p'],partial_rho_comparator=cc['rho'],partial_p_comparator=cc['p'],
    covariates=paste(covars,collapse='+'),outcome_scoring='strict_nonoverlap')
}

# Read exact sample-level TCGA and GSE metadata, then recompute each outcome from expression.
tcga_mod <- read_csv(file.path(root,'results/immune_validation/tcga_sample_scores_with_purity.csv'),show_col_types=FALSE) %>%
  filter(sample_id %in% colnames(expr)) %>% arrange(match(sample_id,colnames(expr)))
stopifnot(!anyDuplicated(tcga_mod$sample_id),all(tcga_mod$sample_id %in% colnames(expr)))
tcga_dat <- tibble(sample_id=tcga_mod$sample_id,purity=as.numeric(tcga_mod$purity),PRII=prii_tcga[match(tcga_mod$sample_id,colnames(expr))])
for(cmp in names(comp_sets)) tcga_dat[[cmp]] <- dataset_scores$TCGA_PAAD_primary[[cmp]][match(tcga_mod$sample_id,colnames(expr))]
gse_mod <- read_csv(file.path(root,'results/external_validation/GSE79668_PSPS_immune_scores.csv'),show_col_types=FALSE) %>%
  filter(sample %in% colnames(gexpr)) %>% arrange(match(sample,colnames(gexpr)))
stopifnot(!anyDuplicated(gse_mod$sample),all(gse_mod$sample %in% colnames(gexpr)))
gse_dat <- tibble(sample=gse_mod$sample,PRII=prii_gse[match(gse_mod$sample,colnames(gexpr))])
for(cmp in names(comp_sets)) gse_dat[[cmp]] <- dataset_scores$GSE79668[[cmp]][match(gse_mod$sample,colnames(gexpr))]

coverage <- list(); published_models <- list()
for(cmp in c('Zuo2022_standardized_directional_proxy','Xie2023_descriptive_3gene_panel')) {
  if(cmp=='curated_inflammation_noncore') {tcga_out <- tcga_sets[setdiff(names(tcga_sets),'inflammation_noncore')]; gse_out <- gse_sets[setdiff(names(gse_sets),'inflammation_noncore')]}
  else {tcga_out <- tcga_sets; gse_out <- gse_sets}
  for(outcome in names(tcga_out)) {
    req <- tcga_out[[outcome]]; present <- intersect(setdiff(req,unique(c(core_h,comp_sets[[cmp]]))),rownames(expr))
    coverage[[length(coverage)+1]] <- tibble(dataset='TCGA_PAAD',outcome=outcome,comparator=cmp,n_requested=length(req),n_after_core_comparator_removal=length(setdiff(req,unique(c(core_h,comp_sets[[cmp]])))),n_expressed_after_removal=length(present),coverage_after_removal=length(present)/length(req),removed_genes=paste(setdiff(req,setdiff(req,unique(c(core_h,comp_sets[[cmp]])))),collapse=';'),remaining_genes=paste(present,collapse=';'),modelable=length(present)>=3)
    if(length(present)>=3) {
      one <- tcga_dat; one$outcome_score <- score_z_expr(expr,present)[match(one$sample_id,colnames(expr))]
      published_models[[length(published_models)+1]] <- fit_one(one,'TCGA_PAAD','outcome_score',cmp,'purity',length(present),length(setdiff(req,setdiff(req,unique(c(core_h,comp_sets[[cmp]])))))) %>% mutate(outcome=.env$outcome)
    }
  }
  for(outcome in names(gse_out)) {
    req <- gse_out[[outcome]]; present <- intersect(setdiff(req,unique(c(core_h,comp_sets[[cmp]]))),rownames(gexpr))
    coverage[[length(coverage)+1]] <- tibble(dataset='GSE79668',outcome=outcome,comparator=cmp,n_requested=length(req),n_after_core_comparator_removal=length(setdiff(req,unique(c(core_h,comp_sets[[cmp]])))),n_expressed_after_removal=length(present),coverage_after_removal=length(present)/length(req),removed_genes=paste(setdiff(req,setdiff(req,unique(c(core_h,comp_sets[[cmp]])))),collapse=';'),remaining_genes=paste(present,collapse=';'),modelable=length(present)>=3)
    if(length(present)>=3) {
      one <- gse_dat; one$outcome_score <- score_z_expr(gexpr,present)[match(one$sample,colnames(gexpr))]
      published_models[[length(published_models)+1]] <- fit_one(one,'GSE79668','outcome_score',cmp,character(),length(present),length(setdiff(req,setdiff(req,unique(c(core_h,comp_sets[[cmp]])))))) %>% mutate(outcome=.env$outcome)
    }
  }
}
pub_cov <- bind_rows(coverage)
write_csv(pub_cov,file.path(auditout,'S10_published_strict_immune_module_coverage.csv'))
pub_models <- bind_rows(published_models) %>% select(dataset,outcome,comparator,everything())
write_csv(pub_models,file.path(auditout,'S10_published_proxy_nested_models.csv'))

# Combine with the pre-existing strict project-defined comparator models and
# recompute BH families across all S10 comparator/outcome tests within each cohort.
project_models <- read_csv(file.path(auditout,'S10_nested_model_results.csv'),show_col_types=FALSE) %>%
  mutate(n_comparator_genes=ifelse(is.na(n_comparator_genes),NA_integer_,as.integer(n_comparator_genes))) %>%
  select(-any_of(c('lrt_p_PRII','fdr_lrt')))
all_models <- bind_rows(project_models,pub_models) %>%
  mutate(fdr_family=ifelse(comparator=='Xie2023_descriptive_3gene_panel','Xie_descriptive_sensitivity','project_sets_plus_Zuo_proxy')) %>%
  group_by(dataset,fdr_family) %>%
  mutate(fdr_nested_F=p.adjust(p_nested_F_PRII,'BH'),fdr_partial=p.adjust(partial_p_PRII,'BH'),fdr_beta=p.adjust(p_beta_PRII,'BH'),fdr_beta_HC3=p.adjust(p_beta_PRII_HC3,'BH'),fdr_nested_F_comparator=p.adjust(p_nested_F_comparator_after_PRII,'BH')) %>% ungroup()
# Alias the primary family for concise downstream use; Xie rows remain separate exploratory tests.
all_models <- all_models %>% mutate(fdr_primary_by_dataset=ifelse(fdr_family=='project_sets_plus_Zuo_proxy',fdr_nested_F,NA_real_),fdr_beta_HC3_primary_by_dataset=ifelse(fdr_family=='project_sets_plus_Zuo_proxy',fdr_beta_HC3,NA_real_),fdr_nested_F_comparator_primary_by_dataset=ifelse(fdr_family=='project_sets_plus_Zuo_proxy',fdr_nested_F_comparator,NA_real_),fdr_Xie_sensitivity=ifelse(fdr_family=='Xie_descriptive_sensitivity',fdr_nested_F,NA_real_))
write_csv(all_models,file.path(auditout,'S10_nested_models_all_comparators.csv'))

# External direction replication is evaluated only by signed conditional effects;
# nested delta-R2 is an improvement magnitude and is never treated as direction.
direction <- all_models %>% filter(dataset %in% c('TCGA_PAAD','GSE79668')) %>%
  select(dataset,outcome,comparator,beta_PRII,partial_rho_PRII,delta_r2_PRII) %>%
  pivot_wider(names_from=dataset,values_from=c(beta_PRII,partial_rho_PRII,delta_r2_PRII),names_sep='__') %>%
  filter(!is.na(beta_PRII__TCGA_PAAD),!is.na(beta_PRII__GSE79668)) %>%
  mutate(beta_same_direction=sign(beta_PRII__TCGA_PAAD)==sign(beta_PRII__GSE79668),
         partial_rho_same_direction=sign(partial_rho_PRII__TCGA_PAAD)==sign(partial_rho_PRII__GSE79668),
         delta_r2_TCGA=delta_r2_PRII__TCGA_PAAD,delta_r2_GSE79668=delta_r2_PRII__GSE79668)
write_csv(direction,file.path(auditout,'S10_direction_consistency_all_comparators.csv'))

# Model-complete-case VIF includes purity for TCGA and is calculated for the actual M2 predictor set.
model_vif <- bind_rows(lapply(names(comp_sets),function(cmp) {
  bind_rows(lapply(names(dataset_scores),function(ds) {
    d <- if(ds=='TCGA_PAAD_primary') tcga_dat else gse_dat
    vars <- c('PRII',cmp,if(ds=='TCGA_PAAD_primary')'purity')
    d <- d[complete.cases(d[,vars,drop=FALSE]),,drop=FALSE]
    bind_rows(lapply(vars,function(term) {
      others <- setdiff(vars,term); mod <- lm(as.formula(paste(term,'~',paste(others,collapse='+'))),data=d)
      r2 <- summary(mod)$r.squared
      tibble(dataset=if(ds=='TCGA_PAAD_primary')'TCGA_PAAD' else 'GSE79668',comparator=cmp,term=term,n=nrow(d),vif=1/(1-r2))
    }))
  }))
}))
write_csv(model_vif,file.path(auditout,'S10_vif_model_complete_cases.csv'))

# Fixed-composition sensitivity: use the same non-overlapping marker subset for each
# selected comparator within a cohort, allowing comparable model fit estimates.
common_comps <- c('curated_pyroptosis_extended_noncore','manual_inflammasome_noncore','curated_inflammation_noncore','Zuo2022_standardized_directional_proxy')
common_drop <- unique(c(core_h,unlist(comp_sets[common_comps])))
common_outcomes <- c('t_cell','cytotoxic','b_cell','myeloid','macrophage_TAM')
common_models <- list(); common_coverage <- list()
for(cmp in common_comps) {
  for(outcome in common_outcomes) {
    for(ds in c('TCGA_PAAD','GSE79668')) {
      req <- if(ds=='TCGA_PAAD')tcga_sets[[outcome]] else gse_sets[[outcome]]
      mat <- if(ds=='TCGA_PAAD')expr else gexpr
      present <- intersect(setdiff(req,common_drop),rownames(mat))
      common_coverage[[length(common_coverage)+1]] <- tibble(dataset=ds,outcome=outcome,comparator=cmp,n_requested=length(req),n_common_remaining=length(setdiff(req,common_drop)),n_expressed=length(present),coverage=length(present)/length(req),genes=paste(present,collapse=';'))
      if(length(present)<3) next
      d <- if(ds=='TCGA_PAAD')tcga_dat else gse_dat
      d$outcome_score <- score_z_expr(mat,present)[match(if(ds=='TCGA_PAAD')d$sample_id else d$sample,colnames(mat))]
      common_models[[length(common_models)+1]] <- fit_one(d,ds,'outcome_score',cmp,if(ds=='TCGA_PAAD')'purity' else character(),length(present),length(req)-length(setdiff(req,common_drop))) %>% mutate(outcome=.env$outcome)
    }
  }
}
common_cov <- bind_rows(common_coverage)
common_res <- bind_rows(common_models) %>% group_by(dataset) %>% mutate(fdr_nested_F=p.adjust(p_nested_F_PRII,'BH'),fdr_partial=p.adjust(partial_p_PRII,'BH'),fdr_beta_HC3=p.adjust(p_beta_PRII_HC3,'BH'),fdr_nested_F_comparator=p.adjust(p_nested_F_comparator_after_PRII,'BH')) %>% ungroup()
write_csv(common_cov,file.path(auditout,'S10_common_outcome_coverage.csv'))
write_csv(common_res,file.path(auditout,'S10_common_composition_nested_models.csv'))

# Unified correlation output with the previously frozen project-defined audit.
project_cor <- read_csv(file.path(auditout,'prii_comparator_correlations.csv'),show_col_types=FALSE)
write_csv(bind_rows(project_cor %>% mutate(signature_class='project_defined'),
                    cor_tbl %>% mutate(signature_class=ifelse(comparator=='Zuo2022_standardized_directional_proxy','published_weighted_proxy','published_descriptive_only'))),
          file.path(auditout,'S10_PRII_signature_correlations_all.csv'))

# Main Figure S10: A) signature overlap/coverage; B) PRII-comparator correlations;
# C) strict conditional partial rho; D) incremental R2. Signed quantities are shown
# for replication; delta-R2 is shown as magnitude only.
show_comps <- c('curated_pyroptosis_extended_noncore','manual_inflammasome_noncore','curated_inflammation_noncore','Zuo2022_standardized_directional_proxy')
labels <- c(curated_pyroptosis_extended_noncore='Curated pyroptosis',manual_inflammasome_noncore='Manual inflammasome',curated_inflammation_noncore='Inflammation baseline',Zuo2022_standardized_directional_proxy='Zuo 2022 proxy')
plot_overlap <- bind_rows(lapply(show_comps,function(cmp) {
  genes <- full_sets[[cmp]]; core_n <- length(intersect(genes,core_h)); noncore <- setdiff(genes,core_h)
  bind_rows(lapply(names(expr_by_ds),function(ds)tibble(dataset=if(ds=='TCGA_PAAD_primary')'TCGA-PAAD' else ds,comparator=cmp,core_overlap=core_n,coverage=length(intersect(noncore,expr_by_ds[[ds]]))/length(noncore),n_noncore=length(noncore))))
})) %>% mutate(comparator_label=factor(labels[comparator],levels=rev(unname(labels[show_comps]))))
pA <- ggplot(plot_overlap,aes(x=coverage,y=comparator_label,colour=dataset)) + geom_point(aes(size=n_noncore),position=position_dodge(width=.55)) +
  geom_text(aes(label=paste0(round(coverage*100),'% (',core_overlap,' core overlap)')),position=position_dodge(width=.55),hjust=-.08,size=3,show.legend=FALSE) +
  scale_x_continuous(limits=c(0,1.35),labels=scales::percent_format(accuracy=1)) + labs(title='A  Full signature overlap and non-core gene coverage',x='Expressed non-core genes / non-core signature genes',y=NULL,colour=NULL,size='Non-core genes') + theme_bw(base_size=10) + theme(legend.position='bottom')
plot_cor <- cor_tbl %>% filter(comparator %in% show_comps) %>% mutate(comparator_label=factor(labels[comparator],levels=rev(unname(labels[show_comps]))),q_label=ifelse(fdr_bh<.001,'q<0.001',paste0('q=',formatC(fdr_bh,format='f',digits=3))),dataset=recode(dataset,TCGA_PAAD_primary='TCGA-PAAD (n=178)',GSE79668='GSE79668 (n=51)'))
pB <- ggplot(plot_cor,aes(x=spearman_rho,y=comparator_label,colour=dataset)) + geom_vline(xintercept=0,lty=2,colour='grey65') + geom_point(size=2.6,position=position_dodge(width=.60)) +
  geom_text(aes(label=sprintf('ρ=%.2f\n%s',spearman_rho,q_label)),position=position_dodge(width=.60),vjust=-.55,size=2.5,show.legend=FALSE) +
  scale_x_continuous(limits=c(-1,1)) + labs(title='B  PRII–comparator correlation (BH q)',x='Spearman ρ',y=NULL,colour=NULL) + theme_bw(base_size=10) + theme(legend.position='bottom')
outcome_labels <- c(t_cell='T cell',cytotoxic='Cytotoxic',b_cell='B cell',myeloid='Myeloid',macrophage_TAM='Macrophage/TAM')
plotmods <- common_res %>% mutate(comparator_label=factor(labels[comparator],levels=rev(unname(labels[common_comps]))),outcome_label=factor(outcome_labels[outcome],levels=rev(unname(outcome_labels[common_outcomes]))),dataset=recode(dataset,TCGA_PAAD='TCGA-PAAD (n=149; purity-adjusted)',GSE79668='GSE79668 (n=51; unadjusted)'))
pC <- ggplot(plotmods,aes(x=outcome_label,y=comparator_label,fill=partial_rho_PRII)) + geom_tile(colour='white',linewidth=.35) + geom_text(aes(label=ifelse(is.na(partial_rho_PRII),'',sprintf('%.2f\n(n=%d genes)',partial_rho_PRII,n_outcome_genes))),size=2.2) +
  facet_wrap(~dataset,nrow=1) + scale_fill_gradient2(low='#3B6FB6',mid='white',high='#B34A4A',midpoint=0,limits=c(-1,1),na.value='grey90') +
  labs(title='C  Strict conditional partial Spearman ρ for PRII',x=NULL,y=NULL,fill='Partial ρ') + theme_bw(base_size=9) + theme(axis.text.x=element_text(angle=35,hjust=1),strip.background=element_rect(fill='grey93'))
plotmods_long <- plotmods %>% pivot_longer(c(delta_r2_PRII,delta_r2_comparator_after_PRII),names_to='increment_term',values_to='incremental_r2') %>% mutate(increment_term=recode(increment_term,delta_r2_PRII='PRII added after comparator',delta_r2_comparator_after_PRII='Comparator added after PRII'))
pD <- ggplot(plotmods_long,aes(x=outcome_label,y=comparator_label,fill=incremental_r2)) + geom_tile(colour='white',linewidth=.35) + geom_text(aes(label=ifelse(is.na(incremental_r2),'',sprintf('%.03f\n(n=%d genes)',incremental_r2,n_outcome_genes))),size=2.0) +
  facet_grid(increment_term~dataset) + scale_fill_gradientn(colours=c('#fff7ec','#fee8c8','#fdbb84','#e34a33','#b30000'),limits=c(0,max(plotmods_long$incremental_r2,na.rm=TRUE)),na.value='grey90') +
  labs(title='D  In-sample incremental R² in either model order',x=NULL,y=NULL,fill='Incremental R²') + theme_bw(base_size=8.5) + theme(axis.text.x=element_text(angle=35,hjust=1),strip.background=element_rect(fill='grey93'))
fig <- (pA | pB) / (pC / pD) + plot_layout(heights=c(1,2.1)) + plot_annotation(theme=theme(plot.title=element_text(face='bold',size=10)))
ggsave(file.path(figout,'Figure_S10_PRII_incremental_association.pdf'),fig,width=15,height=13,device=cairo_pdf)
ggsave(file.path(figout,'Figure_S10_PRII_incremental_association.png'),fig,width=15,height=13,dpi=400,bg='white')
ggsave(file.path(figout,'Figure_S10_PRII_incremental_association.tiff'),fig,width=15,height=13,dpi=300,bg='white')

# Structured run metadata and cautious, bilingual interpretation.
writeLines(capture.output(sessionInfo()),file.path(auditout,'sessionInfo_S10.txt'))
writeLines(c(
  '# S10: PRII 与既有焦亡相关评分的增量关联 / Incremental association of PRII against existing pyroptosis-related scores','',
  '## 方法口径','',
  '- PRII 为 TXNIP/NLRP3/CASP1/GSDMD/IL1B/IL18 六基因队列内逐基因 z 分数均值。比较评分每次单独建模。','- Zuo et al. 2022 的六个公开 LASSO-Cox 系数用于计算 `sum(beta × within-cohort gene z-score)`。原文未充分说明表达单位/预处理，因此这是标准化方向性代理评分，不是原始风险评分复现。Xie et al. 2023 的 ARSG/SLC26A11/MET 仅为无权重描述性面板，不能称为复现的 Xie 风险评分。','- 条件免疫模块按每一组 comparator 重新计分，并剔除 PRII 六基因及当前 comparator 的全部基因。TCGA 模型调整纯度；GSE79668 无有据可用的纯度协变量，报告未调整的条件关联。','- M1 = strict immune module ~ comparator (+ purity in TCGA); M2 = M1 + PRII。并对调顺序计算 comparator 在 PRII 之后的增量拟合。`p_nested_F_PRII` 与 `p_nested_F_comparator_after_PRII` 均是嵌套线性模型的部分 F 检验，不是似然比检验。模型同时报告 HC3 异方差稳健的标准化 PRII 系数区间与 P 值。主 BH 家族包括项目定义比较项与 Zuo 2022 proxy（TCGA 58 项、GSE79668 45 项）；Xie 三基因描述性面板单独作为探索性敏感性家族校正（TCGA 10 项、GSE79668 8 项）。','- 外部方向支持以 partial Spearman ρ 和标准化 β 的符号为准。ΔR² 仅表示样本内拟合增量，不是预测性能，也不作为方向性复制证据。C/D 面板在每个队列内部对比较项固定了 marker set，但 TCGA 与 GSE79668 的对应模块定义及 marker 数仍有差异，因此属于同类模块的方向支持而非完全同构复现。','',
  '## 主要限制','',
  '- 这是批量转录组层面的关联/模型增量分析，不证明功能性焦亡、炎症小体活化、因果免疫调控或独立临床预后价值。','- Zuo proxy 受输入表达尺度未报告的限制；Xie comparator 无法恢复原始系数和公式。外部队列模型仅为 GSE79668 的条件方向复制，不代表预后模型验证。','- 项目定义的广泛通路集与 PRII/免疫模块可能共享生物学成分；严格去重只能降低基因算术重叠，不能消除共同细胞组成/共同转录程序造成的混杂。','',
  '## 输出','',
  '- `S10_nested_models_all_comparators.csv`：项目定义比较评分与可用已发表评分代理的统一嵌套模型，含 PRII-after-comparator 与 comparator-after-PRII 两个增量方向；`fdr_primary_by_dataset` 为项目比较项 + Zuo proxy 的主 BH 家族，Xie 描述性敏感性使用独立 FDR 列。','- `S10_direction_consistency_all_comparators.csv`：TCGA 与 GSE79668 的 conditional effect 方向对照。','- `S10_PRII_signature_correlations_all.csv`、`S10_vif_model_complete_cases.csv`、`S10_signature_gene_overlap_coverage.csv`：相关性、实际模型完整病例 VIF、基因重叠和覆盖。','- `S10_published_strict_immune_module_coverage.csv`、`S10_published_proxy_nested_models.csv`：已发表评分代理的严格条件模型及 marker 覆盖。`S10_common_composition_nested_models.csv` 对选定比较项使用相同的非重叠 marker 子集，供跨比较项公平比较。','- `figures/Figure_S10_PRII_incremental_association.pdf/.png/.tiff`：S10 主图。','- `signatures/published_proxy_correlations.csv` 已替换为按样本计算、每个队列每个评分仅一行的纠正版。','',
  '## English interpretation','',
  'PRII was evaluated as a compact pyroptosis-associated transcriptional proxy against project-defined non-core comparator sets and a standardized directional proxy derived from the published Zuo et al. coefficients. Immune modules were rescored after removing genes shared with PRII and the comparator. Nested-model comparisons quantify the additional transcriptomic association contributed by PRII after the comparator (and purity in TCGA) is included; the reverse comparator increment was also calculated. The primary BH family comprised project-defined comparators and the Zuo proxy, while the Xie three-gene descriptive panel was assessed in a separate exploratory family. GSE79668 provides direction-of-association support without purity adjustment; corresponding modules are analogous but not fully harmonized across cohorts. In-sample ΔR² is an association-fit increment, not predictive performance. These analyses do not establish functional pyroptosis, causal immune regulation, or independent prognostic value. The Zuo proxy is not an exact reproduction of the published score, and the Xie three-gene panel is descriptive only because its coefficients and formula are unavailable.'
),file.path(auditout,'S10_methods_results_zh_en.md'))
# Replace the pre-analysis audit report regenerated by the sourced audit script
# with the final corrected methods/interpretation manifest.
manifest_text <- readLines(file.path(auditout,'S10_methods_results_zh_en.md'),warn=FALSE)
writeLines(c(manifest_text,'','For manuscript-ready bilingual Results, see `S10_results_writeup_zh_en.md`; for the caption, see `Figure_S10_legend_zh_en.md`.'),file.path(auditout,'S10_signature_audit_report_zh_en.md'))
message('Completed strict S10 published comparator audit; rows: ',nrow(pub_models),'; combined models: ',nrow(all_models))
