#!/usr/bin/env Rscript
# S10 conditional immune association and nested model analysis.
# Comparator lists are frozen by s10_signature_audit.R; this script does not
# search or reweight genes and does not claim functional pyroptosis.
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(tibble)})
options(stringsAsFactors=FALSE)
root <- normalizePath('.',mustWork=TRUE)
source(file.path(root,'results/prii_incremental/s10_signature_audit.R'))
outdir <- file.path(root,'results/prii_incremental/s10_signature_audit')

comp_names <- c('curated_pyroptosis_extended_noncore','manual_inflammasome_noncore','curated_inflammation_noncore','project_pyroptosis_noncore','project_inflammasome_noncore')
outcomes_tcga <- c('immune_infiltration','t_cell','cytotoxic','b_cell','myeloid','macrophage_TAM','CAF','EMT','lipid_metabolism','inflammation_noncore')
outcomes_gse <- c('leukocyte','t_cell','cytotoxic','b_cell','myeloid','macrophage_TAM','neutrophil','inflammation_noncore')
score_z_expr <- function(mat, genes) {
  g <- intersect(unique(genes),rownames(mat)); if (!length(g)) return(rep(NA_real_,ncol(mat)))
  zz <- t(apply(mat[g,,drop=FALSE],1,function(x)as.numeric(scale(x)))); zz[!is.finite(zz)]<-0; colMeans(zz,na.rm=TRUE)
}
tcga_out_sets <- list(
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
gse_out_sets <- list(
  leukocyte=c('PTPRC','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','LST1'),
  t_cell=c('CD3D','CD3E','CD3G','CD4','CD8A','IL7R','LTB'),
  cytotoxic=c('NKG7','GZMB','GZMK','PRF1','CTSW','GZMA'),
  b_cell=c('MS4A1','CD79A','CD79B','CD37','CD74','HLA-DRA','HLA-DRB1','CD19','CD83'),
  myeloid=c('LST1','TYROBP','AIF1','CTSS','CTSD','LGALS3','S100A10'),
  macrophage_TAM=c('ADGRE1','CD68','CSF1R','C1QA','C1QB','C1QC','MS4A7','MERTK','MAFB','CD163','MRC1','MARCO','LGALS3'),
  neutrophil=c('S100A8','S100A9','LCN2','CXCR2','CSF3R','NGP','CAMP'),
  inflammation_noncore=c('TNF','IL6','CXCL1','CXCL2','CXCL5','CCL2','CCL3','CCL4','CCL5','NFKB1','RELA','TLR4','TGFB1','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3'))

standard_beta <- function(dat, y, vars) {
  z <- as.data.frame(lapply(dat[,c(y,vars),drop=FALSE], function(x) as.numeric(scale(x))))
  names(z) <- c(y,vars)
  f <- as.formula(paste(y,'~',paste(vars,collapse='+')))
  cc <- summary(lm(f,data=z))$coefficients
  tibble(term=rownames(cc),beta=cc[,1],beta_se=cc[,2],beta_p=cc[,4]) %>% filter(term %in% vars)
}
partial_rank <- function(dat,y,x,covars=character()) {
  v <- dat[,c(y,x,covars),drop=FALSE]; v <- v[complete.cases(v),,drop=FALSE]
  if (nrow(v)<8) return(tibble(n=nrow(v),rho=NA_real_,p_value=NA_real_))
  rx <- rank(v[[x]]); ry <- rank(v[[y]])
  if (length(covars)) {z <- as.data.frame(lapply(v[,covars,drop=FALSE],rank)); rx <- residuals(lm(rx~.,data=z)); ry <- residuals(lm(ry~.,data=z))}
  ct <- suppressWarnings(cor.test(rx,ry,method='pearson'))
  tibble(n=nrow(v),rho=unname(ct$estimate),p_value=ct$p.value)
}
nested_one <- function(dat, dataset, outcome, comparator, covars=character(), outcome_label=outcome, n_outcome_genes=NA_integer_, n_comparator_genes=NA_integer_, n_removed=NA_integer_) {
  vars <- c(outcome,comparator,'PRII',covars); d <- dat[,vars,drop=FALSE]; d <- d[complete.cases(d),,drop=FALSE]
  if (nrow(d)<12) return(tibble())
  f0 <- as.formula(paste(outcome,'~',ifelse(length(covars),paste(covars,collapse='+'),'1')))
  f1 <- as.formula(paste(outcome,'~',paste(c(covars,comparator),collapse='+')))
  f_prii <- as.formula(paste(outcome,'~',paste(c(covars,'PRII'),collapse='+')))
  f2 <- as.formula(paste(outcome,'~',paste(c(covars,comparator,'PRII'),collapse='+')))
  m0<-lm(f0,d); m1<-lm(f1,d); m_prii<-lm(f_prii,d); m2<-lm(f2,d)
  rss1<-deviance(m1); rss2<-deviance(m2); r2_1<-summary(m1)$r.squared; r2_2<-summary(m2)$r.squared
  f_prii_test <- anova(m1,m2); f_comp_test <- anova(m_prii,m2)
  delta_r2 <- r2_2-r2_1; partial_r2 <- (rss1-rss2)/rss1
  sb <- standard_beta(d,outcome,c(covars,comparator,'PRII'))
  zdata <- as.data.frame(lapply(d[,c(outcome,c(covars,comparator,'PRII')),drop=FALSE],function(x)as.numeric(scale(x))))
  names(zdata) <- c(outcome,c(covars,comparator,'PRII'))
  mz <- lm(f2,zdata); X <- model.matrix(mz); e <- residuals(mz); h <- hatvalues(mz)
  bread <- solve(crossprod(X)); meat <- crossprod(X * as.numeric(e/(1-h))); Vhc3 <- bread %*% meat %*% bread
  se_hc3 <- sqrt(diag(Vhc3)); beta_hc3 <- coef(mz)['PRII']; se_prii_hc3 <- se_hc3[match('PRII',names(se_hc3))]
  df_hc3 <- df.residual(mz); p_prii_hc3 <- 2*pt(-abs(beta_hc3/se_prii_hc3),df=df_hc3); tcrit <- qt(.975,df_hc3)
  bcomp_hc3 <- coef(mz)[comparator]; se_comp_hc3 <- se_hc3[match(comparator,names(se_hc3))]
  p_comp_hc3 <- 2*pt(-abs(bcomp_hc3/se_comp_hc3),df=df_hc3)
  bp <- sb$beta[match('PRII',sb$term)]; bcomp<-sb$beta[match(comparator,sb$term)]; pp<-sb$beta_p[match('PRII',sb$term)]; pc<-sb$beta_p[match(comparator,sb$term)]
  pr <- partial_rank(d,outcome,'PRII',c(covars,comparator)); cr <- partial_rank(d,outcome,comparator,c(covars,'PRII'))
  tibble(dataset=dataset,outcome=outcome_label,comparator=comparator,n=nrow(d),n_outcome_genes=n_outcome_genes,n_comparator_genes=n_comparator_genes,n_removed_from_outcome=n_removed,r2_M1=r2_1,r2_M2=r2_2,delta_r2_PRII=delta_r2,partial_r2_PRII=partial_r2,p_nested_F_PRII=f_prii_test$`Pr(>F)`[2],delta_r2_comparator_after_PRII=summary(m2)$r.squared-summary(m_prii)$r.squared,partial_r2_comparator_after_PRII=(deviance(m_prii)-deviance(m2))/deviance(m_prii),p_nested_F_comparator_after_PRII=f_comp_test$`Pr(>F)`[2],delta_AIC_PRII=AIC(m2)-AIC(m1),delta_BIC_PRII=BIC(m2)-BIC(m1),beta_PRII=bp,p_beta_PRII=pp,beta_PRII_HC3=beta_hc3,beta_PRII_HC3_SE=se_prii_hc3,beta_PRII_HC3_CI_low=beta_hc3-tcrit*se_prii_hc3,beta_PRII_HC3_CI_high=beta_hc3+tcrit*se_prii_hc3,p_beta_PRII_HC3=p_prii_hc3,beta_comparator=bcomp,p_beta_comparator=pc,beta_comparator_HC3=bcomp_hc3,beta_comparator_HC3_SE=se_comp_hc3,p_beta_comparator_HC3=p_comp_hc3,partial_rho_PRII=pr$rho,partial_p_PRII=pr$p_value,partial_rho_comparator=cr$rho,partial_p_comparator=cr$p_value,covariates=paste(covars,collapse='+'),outcome_scoring='strict_nonoverlap')
}

# TCGA bulk module scores and exact sample barcode matching.
tcga_mod <- read_csv(file.path(root,'results/immune_validation/tcga_sample_scores_with_purity.csv'),show_col_types=FALSE)
tcga_mod <- tcga_mod %>% filter(sample_id %in% rownames(tcga_comp)) %>% arrange(match(sample_id,rownames(tcga_comp)))
tcga_dat <- bind_cols(tcga_mod,as.data.frame(tcga_comp[tcga_mod$sample_id,,drop=FALSE]))
stopifnot(all(tcga_dat$sample_id==rownames(tcga_comp)[match(tcga_dat$sample_id,rownames(tcga_comp))]))
comp_gene_sets <- list(
  curated_pyroptosis_extended_noncore=setdiff(curated_pyro_h,core_h),
  manual_inflammasome_noncore=setdiff(manual_infl_h,core_h),
  curated_inflammation_noncore=curated_infl_noncore_h,
  project_pyroptosis_noncore=setdiff(gs_h$pyroptosis,core_h),
  project_inflammasome_noncore=setdiff(gs_h$inflammasome,core_h))
tcga_rows <- list()
coverage_rows <- list()
for (cmp in comp_names) for (out in outcomes_tcga) {
  # The same marker set cannot serve as both predictor and outcome.
  if (cmp=='curated_inflammation_noncore' && out=='inflammation_noncore') next
  req <- unique(tcga_out_sets[[out]])
  strict <- setdiff(req,unique(c(core_h,comp_gene_sets[[cmp]])))
  present <- intersect(strict,rownames(expr))
  coverage_rows[[length(coverage_rows)+1]] <- tibble(dataset='TCGA_PAAD',outcome=out,comparator=cmp,n_requested=length(req),n_after_core_comparator_removal=length(strict),n_expressed_after_removal=length(present),coverage_after_removal=length(present)/length(req),removed_genes=paste(setdiff(req,strict),collapse=';'),remaining_genes=paste(present,collapse=';'),modelable=length(present)>=3)
  if (length(present)<3) next
  one <- tcga_dat
  one$outcome_score <- score_z_expr(expr,present)[match(one$sample_id,colnames(expr))]
  tcga_rows[[length(tcga_rows)+1]] <- nested_one(one,'TCGA_PAAD','outcome_score',cmp,covars='purity',outcome_label=out,n_outcome_genes=length(present),n_comparator_genes=length(intersect(comp_gene_sets[[cmp]],rownames(expr))),n_removed=length(setdiff(req,strict)))
}
tcga_res <- bind_rows(tcga_rows) %>% mutate(fdr_nested_F=p.adjust(p_nested_F_PRII,'BH'),fdr_partial=p.adjust(partial_p_PRII,'BH'),fdr_beta=p.adjust(p_beta_PRII,'BH'),fdr_beta_HC3=p.adjust(p_beta_PRII_HC3,'BH'),fdr_nested_F_comparator=p.adjust(p_nested_F_comparator_after_PRII,'BH'))

# External PAAD cohort: no documented purity covariate, therefore unadjusted
# conditional associations are the prespecified replication endpoint.
gse_mod <- read_csv(file.path(root,'results/external_validation/GSE79668_PSPS_immune_scores.csv'),show_col_types=FALSE)
gse_mod <- gse_mod %>% filter(sample %in% rownames(g_comp)) %>% arrange(match(sample,rownames(g_comp)))
gse_dat <- bind_cols(gse_mod,as.data.frame(g_comp[gse_mod$sample,,drop=FALSE]))
gse_rows <- list()
for (cmp in comp_names) for (out in outcomes_gse) {
  if (cmp=='curated_inflammation_noncore' && out=='inflammation_noncore') next
  req <- unique(gse_out_sets[[out]])
  strict <- setdiff(req,unique(c(core_h,comp_gene_sets[[cmp]])))
  present <- intersect(strict,rownames(gexpr))
  coverage_rows[[length(coverage_rows)+1]] <- tibble(dataset='GSE79668',outcome=out,comparator=cmp,n_requested=length(req),n_after_core_comparator_removal=length(strict),n_expressed_after_removal=length(present),coverage_after_removal=length(present)/length(req),removed_genes=paste(setdiff(req,strict),collapse=';'),remaining_genes=paste(present,collapse=';'),modelable=length(present)>=3)
  if (length(present)<3) next
  one <- gse_dat
  one$outcome_score <- score_z_expr(gexpr,present)[match(one$sample,colnames(gexpr))]
  gse_rows[[length(gse_rows)+1]] <- nested_one(one,'GSE79668','outcome_score',cmp,covars=character(),outcome_label=out,n_outcome_genes=length(present),n_comparator_genes=length(intersect(comp_gene_sets[[cmp]],rownames(gexpr))),n_removed=length(setdiff(req,strict)))
}
write_csv(bind_rows(coverage_rows),file.path(outdir,'S10_strict_immune_module_coverage.csv'))
gse_res <- bind_rows(gse_rows) %>% mutate(fdr_nested_F=p.adjust(p_nested_F_PRII,'BH'),fdr_partial=p.adjust(partial_p_PRII,'BH'),fdr_beta=p.adjust(p_beta_PRII,'BH'),fdr_beta_HC3=p.adjust(p_beta_PRII_HC3,'BH'),fdr_nested_F_comparator=p.adjust(p_nested_F_comparator_after_PRII,'BH'))
all_res <- bind_rows(tcga_res,gse_res)
write_csv(all_res,file.path(outdir,'S10_nested_model_results.csv'))

# Directional replication uses signed effects only; nested delta-R2 is a
# non-negative improvement magnitude by construction and is not a direction test.
direction <- all_res %>% select(dataset,outcome,comparator,beta_PRII,partial_rho_PRII,delta_r2_PRII,partial_r2_PRII) %>%
  pivot_wider(names_from=dataset,values_from=c(beta_PRII,partial_rho_PRII,delta_r2_PRII,partial_r2_PRII),names_sep='__') %>%
  mutate(beta_same_direction=sign(beta_PRII__TCGA_PAAD)==sign(beta_PRII__GSE79668),partial_rho_same_direction=sign(partial_rho_PRII__TCGA_PAAD)==sign(partial_rho_PRII__GSE79668),direction_consistent=beta_same_direction & partial_rho_same_direction)
write_csv(direction,file.path(outdir,'S10_direction_consistency.csv'))

manifest <- c(
  '# S10 conditional association and nested models',
  '',
  '- PRII and comparator genes were frozen by `s10_signature_audit.R`; comparator scores are within-cohort gene-wise z-score means.',
  '- TCGA-PAAD outcomes: 10 marker modules from the existing immune validation matrix. Models use the exact primary-tumour complete cases for each outcome/comparator and purity as the prespecified covariate.',
  '- GSE79668 outcomes: 8 prespecified immune modules. No purity or survival covariate was available; external replication is therefore unadjusted conditional association only.',
  '- To avoid arithmetic overlap, each immune outcome is rescored after removing all six PRII genes and all genes in the comparator being tested. Coverage and removed genes are recorded in `S10_strict_immune_module_coverage.csv`; model rows with fewer than three retained genes are omitted.',
  '- M1 = outcome ~ covariates + comparator; M2 = M1 + PRII. The reverse nested increment of the comparator after PRII is also reported. `anova(m1,m2)` is a partial F test, not a likelihood-ratio test. HC3 robust standard errors and confidence intervals are supplied for standardized coefficients; partial Spearman is computed from rank residuals.',
  '- Comparators are fitted one at a time. The same `inflammation_noncore` marker set is not used simultaneously as predictor and outcome.',
  '- ΔR² is an incremental fit magnitude and is always nonnegative in nested OLS; it is not evidence of directional replication. This analysis tests transcriptomic association only; it cannot establish functional pyroptosis, inflammasome activation, causal immune communication or independent prognosis.',
  '',
  'Outputs: `S10_nested_model_results.csv`; `S10_direction_consistency.csv`.'
)
writeLines(manifest,file.path(outdir,'S10_nested_model_manifest_zh_en.md'))
message('Completed S10 nested models: ',nrow(all_res),' rows')
