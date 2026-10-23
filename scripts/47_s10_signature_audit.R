#!/usr/bin/env Rscript
# S10 pre-analysis audit: freeze and compare PRII with available
# pyroptosis/inflammasome comparators.  This audit does not edit manuscripts.
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
outdir <- file.path(root, 'results/prii_incremental/s10_signature_audit')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

core_h <- c('TXNIP','NLRP3','CASP1','GSDMD','IL1B','IL18')
core_m <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
human_long <- read_csv(file.path(root,'scripts/results/gene_sets_human/gene_sets_long.csv'), show_col_types=FALSE)
mouse_long <- read_csv(file.path(root,'scripts/results/gene_sets_mouse/gene_sets_long.csv'), show_col_types=FALSE)
gs_h <- lapply(split(human_long$gene_symbol, human_long$gene_set), unique)
gs_m <- lapply(split(mouse_long$gene_symbol, mouse_long$gene_set), unique)

# Explicit curated comparator already used in TCGA validation.  This is a
# project-defined comparator, not a published weighted signature.
curated_pyro_h <- c('TXNIP','NLRP3','PYCARD','CASP1','CASP4','CASP5','GSDMD','GSDME','IL1B','IL18','AIM2','NLRC4','NLRP1','GBP1','GBP5')
curated_infl_h <- c('NLRP3','CASP1','GSDMD','IL1B','IL18','TNF','IL6','CXCL8','CCL2','CCL20','NFKB1','RELA','TLR4','TGFB1','CXCL1','CXCL2','CXCL5','CXCL10','CXCL11','PTGS2','S100A8','S100A9','NFKBIA','TNFAIP3')
curated_infl_noncore_h <- setdiff(curated_infl_h, core_h)
manual_infl_h <- c(core_h,'AIM2','NLRC4','NLRP1','MEFV','P2RX7','NEK7','CARD8','PYCARD')

zmean <- function(expr, genes) {
  genes <- intersect(unique(genes), rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  zz <- t(vapply(seq_len(nrow(expr[genes,,drop=FALSE])), function(i) {
    x <- as.numeric(expr[genes,,drop=FALSE][i,]); as.numeric(scale(x))
  }, numeric(ncol(expr))))
  if (is.null(dim(zz))) zz <- matrix(zz, nrow=1)
  colMeans(zz, na.rm=TRUE)
}
safe_cor <- function(x,y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok)<5 || sd(x[ok])==0 || sd(y[ok])==0) return(c(rho=NA_real_,p=NA_real_,n=sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok],y[ok], method='spearman', exact=FALSE))
  c(rho=unname(ct$estimate),p=ct$p.value,n=sum(ok))
}
vif_manual <- function(df) {
  nm <- names(df); out <- lapply(nm, function(x) {
    oth <- setdiff(nm,x); ok <- complete.cases(df[,c(x,oth),drop=FALSE])
    if (length(oth)==0 || sum(ok)<5) return(tibble(variable=x,vif=NA_real_,n=sum(ok),r2=NA_real_))
    fit <- lm(df[[x]][ok] ~ ., data=as.data.frame(df[ok,oth,drop=FALSE]))
    r2 <- summary(fit)$r.squared
    tibble(variable=x,vif=ifelse(r2>=1,Inf,1/(1-r2)),n=sum(ok),r2=r2)
  })
  bind_rows(out)
}
overlap_row <- function(dataset, comparator, genes, core, expressed=NULL) {
  genes <- unique(genes); inter <- intersect(genes,core); present <- if (is.null(expressed)) character() else intersect(genes,expressed)
  tibble(dataset=dataset, comparator=comparator, n_requested=length(genes), n_core_overlap=length(inter),
         core_overlap=paste(inter,collapse=';'), n_strict=length(setdiff(genes,core)),
         jaccard_core=length(inter)/length(unique(c(genes,core))), n_present=length(present),
         coverage=ifelse(is.null(expressed),NA_real_,length(present)/length(genes)),
         present_genes=ifelse(length(genes)<=50,paste(present,collapse=';'),NA_character_))
}
score_cor_table <- function(dataset, scores, score_names, prii='PRII') {
  bind_rows(lapply(setdiff(score_names,prii), function(nm) {
    z <- safe_cor(scores[[prii]], scores[[nm]])
    tibble(dataset=dataset, prii=prii, comparator=nm, n=z[['n']], spearman_rho=z[['rho']], p_value=z[['p']])
  })) %>% mutate(fdr_bh=p.adjust(p_value,'BH'))
}

all_overlap <- list(); all_cor <- list(); all_vif <- list(); all_scores <- list()
all_pair_vif <- list()

# ---- Mouse discovery: reuse frozen scores for strict non-core comparators. ----
m_scores <- read_csv(file.path(root,'results/12wk_primary/robustness/12wk_noncore_scores.csv'), show_col_types=FALSE)
m_ps <- read_csv(file.path(root,'results/12wk_primary/scores/12wk_sample_scores.csv'), show_col_types=FALSE)
m_scores <- m_scores %>% left_join(m_ps %>% select(sample,PSPS_core), by='sample')
m_comp <- as.data.frame(m_scores %>% transmute(PRII=PSPS_core, pyroptosis_noncore=pyroptosis_noncore_score,
                                  inflammasome_noncore=inflammasome_noncore_score, inflammation_noncore=inflammation_noncore_score))
rownames(m_comp) <- m_scores$sample
all_scores[['GSE292815_KC_12wk']] <- m_comp
all_cor[['GSE292815_KC_12wk']] <- score_cor_table('GSE292815_KC_12wk',m_comp,names(m_comp))
all_vif[['GSE292815_KC_12wk']] <- vif_manual(m_comp)
all_pair_vif[['GSE292815_KC_12wk']] <- bind_rows(lapply(setdiff(names(m_comp),'PRII'),function(nm){
  x<-m_comp$PRII; y<-m_comp[[nm]]; r<-cor(x,y,use='complete.obs'); tibble(dataset='GSE292815_KC_12wk',comparator=nm,pearson_r=r,vif=1/(1-r^2),n=sum(is.finite(x)&is.finite(y)))
}))
mouse_obj <- readRDS(file.path(root,'results/qc/GSE292815_KC_12wk_matrix.rds'))
mouse_genes <- rownames(mouse_obj$expr)
all_overlap[['GSE292815_KC_12wk']] <- bind_rows(
  overlap_row('GSE292815_KC_12wk','project_pyroptosis_full',gs_m$pyroptosis,core_m,mouse_genes),
  overlap_row('GSE292815_KC_12wk','project_pyroptosis_noncore',setdiff(gs_m$pyroptosis,core_m),core_m,mouse_genes),
  overlap_row('GSE292815_KC_12wk','project_inflammasome_full',gs_m$inflammasome,core_m,mouse_genes),
  overlap_row('GSE292815_KC_12wk','project_inflammasome_noncore',setdiff(gs_m$inflammasome,core_m),core_m,mouse_genes),
  overlap_row('GSE292815_KC_12wk','inflammation_noncore',c('Tnf','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl2','Ccl3','Ccl4','Ccl5','Nfkb1','Rela','Tlr4','Tgfb1','Ptgs2','S100a8','S100a9','Nfkbia','Tnfaip3'),core_m,mouse_genes))

# ---- TCGA-PAAD: recompute project and curated strict scores on exact primary tumors. ----
se <- readRDS(file.path(root,'data/raw/TCGA_PAAD_RNAseq.rds'))
expr0 <- assay(se,'tpm_unstrand'); rd <- as.data.frame(rowData(se)); gene <- toupper(as.character(rd$gene_name)); gene[is.na(gene)|gene==''] <- rownames(expr0)[is.na(gene)|gene=='']; gene <- sub('\\..*$','',gene)
keep <- !duplicated(gene); expr <- log2(expr0[keep,,drop=FALSE]+1); rownames(expr) <- gene[keep]
cd <- as.data.frame(colData(se)); tumor <- (cd$sample_type=='Primary Tumor' | cd$shortLetterCode=='TP'); tumor[is.na(tumor)] <- FALSE; expr <- expr[,tumor,drop=FALSE]
tcga_sets <- list(
  project_pyroptosis_noncore=setdiff(gs_h$pyroptosis,core_h),
  project_inflammasome_noncore=setdiff(gs_h$inflammasome,core_h),
  curated_pyroptosis_extended_noncore=setdiff(curated_pyro_h,core_h),
  curated_inflammation_noncore=curated_infl_noncore_h,
  manual_inflammasome_noncore=setdiff(manual_infl_h,core_h))
tcga_comp <- as.data.frame(setNames(lapply(tcga_sets,function(g)zmean(expr,g)),names(tcga_sets)))
tcga_comp$PRII <- zmean(expr,core_h); tcga_comp <- tcga_comp %>% select(PRII,everything())
rownames(tcga_comp) <- colnames(expr)
all_scores[['TCGA_PAAD_primary']] <- tcga_comp
all_cor[['TCGA_PAAD_primary']] <- score_cor_table('TCGA_PAAD_primary',tcga_comp,names(tcga_comp))
all_vif[['TCGA_PAAD_primary']] <- vif_manual(tcga_comp)
all_pair_vif[['TCGA_PAAD_primary']] <- bind_rows(lapply(setdiff(names(tcga_comp),'PRII'),function(nm){
  x<-tcga_comp$PRII; y<-tcga_comp[[nm]]; r<-cor(x,y,use='complete.obs'); tibble(dataset='TCGA_PAAD_primary',comparator=nm,pearson_r=r,vif=1/(1-r^2),n=sum(is.finite(x)&is.finite(y)))
}))
all_overlap[['TCGA_PAAD_primary']] <- bind_rows(
  overlap_row('TCGA_PAAD_primary','project_pyroptosis_full',gs_h$pyroptosis,core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','project_pyroptosis_noncore',setdiff(gs_h$pyroptosis,core_h),core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','project_inflammasome_full',gs_h$inflammasome,core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','project_inflammasome_noncore',setdiff(gs_h$inflammasome,core_h),core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','curated_pyroptosis_extended_full',curated_pyro_h,core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','curated_pyroptosis_extended_noncore',setdiff(curated_pyro_h,core_h),core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','curated_inflammation_noncore',curated_infl_noncore_h,core_h,rownames(expr)),
  overlap_row('TCGA_PAAD_primary','manual_inflammasome_noncore',setdiff(manual_infl_h,core_h),core_h,rownames(expr)))

# ---- GSE79668: recompute using same human gene lists and independent log2(CPM+1). ----
raw <- read.delim(gzfile(file.path(root,'data/external/GSE79668_counts.txt.gz')),check.names=FALSE,stringsAsFactors=FALSE)
sym <- toupper(sub('_ENSG.*$','',as.character(raw[[1]]))); sample_cols <- names(raw)[-1]
mat <- as.matrix(raw[,sample_cols]); storage.mode(mat)<-'numeric'; rownames(mat)<-sym; mat <- rowsum(mat,group=sym,reorder=FALSE)
lib <- colSums(mat); gexpr <- log2(sweep(mat,2,lib,'/')*1e6+1)
g_sets <- list(project_pyroptosis_noncore=setdiff(gs_h$pyroptosis,core_h), project_inflammasome_noncore=setdiff(gs_h$inflammasome,core_h), curated_pyroptosis_extended_noncore=setdiff(curated_pyro_h,core_h), curated_inflammation_noncore=curated_infl_noncore_h, manual_inflammasome_noncore=setdiff(manual_infl_h,core_h))
g_comp <- as.data.frame(setNames(lapply(g_sets,function(g)zmean(gexpr,g)),names(g_sets))); g_comp$PRII <- zmean(gexpr,core_h); g_comp <- g_comp %>% select(PRII,everything())
rownames(g_comp) <- colnames(gexpr)
all_scores[['GSE79668']] <- g_comp
all_cor[['GSE79668']] <- score_cor_table('GSE79668',g_comp,names(g_comp))
all_vif[['GSE79668']] <- vif_manual(g_comp)
all_pair_vif[['GSE79668']] <- bind_rows(lapply(setdiff(names(g_comp),'PRII'),function(nm){
  x<-g_comp$PRII; y<-g_comp[[nm]]; r<-cor(x,y,use='complete.obs'); tibble(dataset='GSE79668',comparator=nm,pearson_r=r,vif=1/(1-r^2),n=sum(is.finite(x)&is.finite(y)))
}))
all_overlap[['GSE79668']] <- bind_rows(
  overlap_row('GSE79668','project_pyroptosis_full',gs_h$pyroptosis,core_h,rownames(gexpr)),
  overlap_row('GSE79668','project_pyroptosis_noncore',setdiff(gs_h$pyroptosis,core_h),core_h,rownames(gexpr)),
  overlap_row('GSE79668','project_inflammasome_full',gs_h$inflammasome,core_h,rownames(gexpr)),
  overlap_row('GSE79668','project_inflammasome_noncore',setdiff(gs_h$inflammasome,core_h),core_h,rownames(gexpr)),
  overlap_row('GSE79668','curated_pyroptosis_extended_full',curated_pyro_h,core_h,rownames(gexpr)),
  overlap_row('GSE79668','curated_pyroptosis_extended_noncore',setdiff(curated_pyro_h,core_h),core_h,rownames(gexpr)),
  overlap_row('GSE79668','curated_inflammation_noncore',curated_infl_noncore_h,core_h,rownames(gexpr)),
  overlap_row('GSE79668','manual_inflammasome_noncore',setdiff(manual_infl_h,core_h),core_h,rownames(gexpr)))

ov <- bind_rows(all_overlap); cor_tbl <- bind_rows(all_cor); vif_tbl <- bind_rows(lapply(names(all_vif),function(nm) all_vif[[nm]] %>% mutate(dataset=nm,.before=1)))
pair_vif_tbl <- bind_rows(all_pair_vif)
write_csv(ov,file.path(outdir,'signature_overlap_coverage.csv')); write_csv(cor_tbl,file.path(outdir,'prii_comparator_correlations.csv')); write_csv(vif_tbl %>% mutate(vif_scope='all_candidate_scores_entered_together'),file.path(outdir,'vif_all_candidate_scores.csv')); write_csv(pair_vif_tbl,file.path(outdir,'vif_prii_plus_single_comparator.csv'))
set_overlap <- function(dataset, a_name, a, b_name, b) {
  a <- unique(a); b <- unique(b); i <- intersect(a,b)
  tibble(dataset=dataset,set_a=a_name,set_b=b_name,n_a=length(a),n_b=length(b),n_intersection=length(i),jaccard=length(i)/length(union(a,b)),intersection_genes=ifelse(length(i)<=100,paste(i,collapse=';'),NA_character_))
}
set_ov <- bind_rows(
  set_overlap('human','project_pyroptosis',gs_h$pyroptosis,'project_inflammasome',gs_h$inflammasome),
  set_overlap('mouse','project_pyroptosis',gs_m$pyroptosis,'project_inflammasome',gs_m$inflammasome),
  set_overlap('human','project_pyroptosis',gs_h$pyroptosis,'curated_pyroptosis_extended',curated_pyro_h),
  set_overlap('human','project_inflammasome',gs_h$inflammasome,'manual_inflammasome',manual_infl_h),
  set_overlap('human','project_pyroptosis',gs_h$pyroptosis,'curated_inflammation',curated_infl_h))
write_csv(set_ov,file.path(outdir,'comparator_set_overlap.csv'))
frozen_lists <- bind_rows(
  tibble(signature='project_pyroptosis_full_human',role='descriptive_broad',gene=gs_h$pyroptosis),
  tibble(signature='project_pyroptosis_noncore_human',role='sensitivity_only',gene=setdiff(gs_h$pyroptosis,core_h)),
  tibble(signature='project_inflammasome_full_human',role='descriptive_broad',gene=gs_h$inflammasome),
  tibble(signature='project_inflammasome_noncore_human',role='sensitivity_only',gene=setdiff(gs_h$inflammasome,core_h)),
  tibble(signature='curated_pyroptosis_extended_full_human',role='sensitivity_full_overlap',gene=curated_pyro_h),
  tibble(signature='curated_pyroptosis_extended_noncore_human',role='primary_candidate_project_defined',gene=setdiff(curated_pyro_h,core_h)),
  tibble(signature='manual_inflammasome_full_human',role='sensitivity_full_overlap',gene=manual_infl_h),
  tibble(signature='manual_inflammasome_noncore_human',role='primary_candidate_project_defined',gene=setdiff(manual_infl_h,core_h)),
  tibble(signature='curated_inflammation_full_human',role='secondary_baseline_full_overlap',gene=curated_infl_h),
  tibble(signature='curated_inflammation_noncore_human',role='secondary_general_inflammation_baseline',gene=curated_infl_noncore_h),
  tibble(signature='project_pyroptosis_noncore_mouse',role='sensitivity_only',gene=setdiff(gs_m$pyroptosis,core_m)),
  tibble(signature='project_inflammasome_noncore_mouse',role='sensitivity_only',gene=setdiff(gs_m$inflammasome,core_m)))
write_csv(frozen_lists,file.path(outdir,'frozen_comparator_gene_lists.csv'))

cor_wide <- cor_tbl %>% select(dataset,comparator,spearman_rho,p_value,fdr_bh) %>% tidyr::pivot_wider(names_from=comparator,values_from=c(spearman_rho,p_value,fdr_bh))
audit_lines <- c(
  '# S10 comparator pre-analysis audit',
  '',
  '## 中文结论',
  '- PRII固定为队列内六基因逐基因z分数均值：TXNIP/NLRP3/CASP1/GSDMD/IL1B/IL18；小鼠队列使用对应小鼠基因符号。',
  '- 项目现有文件中未找到带 DOI、并按原始基因/权重冻结的已发表焦亡或炎症小体评分。当前比较对象均应称为“项目定义的基因集/评分”，不可称作外部已发表签名。',
  '- 三个队列中可用比较评分与PRII的相关方向全部为正：GSE292815 n=21，rho=.651至.892；TCGA-PAAD n=178，rho=.627至.724；GSE79668 n=51，rho=.788至.843。因此可以继续S10条件关联、TCGA嵌套模型和GSE79668外部方向复制。',
  '- PRII与单个比较评分组成二变量预测部分时，pairwise VIF为1.62–4.02（多数<4）；但多个候选评分同时入模时，广义pyroptosis/inflammasome评分的VIF约55–61，不能把所有评分并列放入同一模型。',
  '- 主模型候选：项目定义的15基因 extended-pyroptosis评分去除PRII六基因后剩余9基因（GSE79668覆盖8/9）；项目定义的manual inflammasome去核心评分8基因（两个PAAD队列均8/8覆盖）。一般炎症noncore评分作为次级基线；极宽MSigDB/manual综合集合只作敏感性/描述性分析。',
  '- 项目自建mouse集合中pyroptosis与inflammasome高度重叠：mouse为44/46，human为98/98 inflammasome基因均在broad pyroptosis内；两者不应被解释为独立通路复现。',
  '',
  '中文输出文件包括 `signature_overlap_coverage.csv`、`comparator_set_overlap.csv`、`frozen_comparator_gene_lists.csv`、`prii_comparator_correlations.csv`、`vif_prii_plus_single_comparator.csv` 和 `vif_all_candidate_scores.csv`。',
  '',
  '## English summary',
  '## Scope and naming',
  '- PRII is frozen as the six-gene, within-cohort gene-wise z-score mean: TXNIP/NLRP3/CASP1/GSDMD/IL1B/IL18 (human symbols; orthologous mouse symbols in GSE292815).',
  '- Current project files contain no DOI-linked, literature-frozen, weighted pyroptosis or inflammasome signature. The MSigDB/manual composites and the 15-gene TCGA extended set are project-defined comparator sets and should not be described as published signatures.',
  '- Strict comparators remove all six PRII genes. Full signatures are sensitivity analyses only.',
  '',
  '## Directional audit',
  capture.output(print(cor_tbl %>% mutate(spearman_rho=round(spearman_rho,3),p_value=format.pval(p_value,digits=3,eps=1e-300),fdr_bh=format.pval(fdr_bh,digits=3,eps=1e-300)),row.names=FALSE)),
  '',
  'All available comparator scores correlated positively with PRII in the discovery cohort, TCGA-PAAD, and GSE79668. This supports proceeding to S10 conditional association and nested model analysis, but it does not establish independent biology.',
  '',
  '## Collinearity',
  capture.output(print(pair_vif_tbl %>% mutate(pearson_r=round(pearson_r,3),vif=round(vif,2)),row.names=FALSE)),
  '',
  'Pairwise PRII-plus-one-comparator VIF values are generally below 5 (the inflammation_noncore comparator is 4.02 in the mouse discovery cohort). The all-candidate-score VIF audit is much higher for broad composite sets (approximately 55–61 in TCGA and GSE79668), so broad pyroptosis/inflammasome composites must not be entered simultaneously.',
  '',
  '## Comparator selection for S10',
  '1. Primary: curated_pyroptosis_extended_noncore (9 genes; project-defined TCGA comparator) and manual_inflammasome_noncore (8 genes; project-defined comparator), one at a time.',
  '2. Secondary: curated_inflammation_noncore (19 genes) as a general inflammatory baseline.',
  '3. Sensitivity only: project_pyroptosis_noncore and project_inflammasome_noncore (very broad MSigDB/manual composites; 801/40 mouse genes, 1438/92 human genes).',
  '4. Do not include all comparator scores in one model; do not call any current comparator an external published signature until a DOI and exact gene list/weights are frozen.',
  '',
  'Generated files: signature_overlap_coverage.csv; comparator_set_overlap.csv; frozen_comparator_gene_lists.csv; prii_comparator_correlations.csv; vif_prii_plus_single_comparator.csv; vif_all_candidate_scores.csv.'
)
writeLines(audit_lines,file.path(outdir,'S10_signature_audit_report_zh_en.md'))
write_csv(bind_rows(lapply(names(all_scores),function(nm) {x<-all_scores[[nm]]; x$sample<-rownames(x); x$dataset<-nm; x})),file.path(outdir,'score_matrix_audit_long.csv'))
writeLines(c(
  'S10 pre-analysis signature audit',
  'Frozen PRII: mean of within-cohort gene-wise z scores for TXNIP/NLRP3/CASP1/GSDMD/IL1B/IL18 (human) or orthologous mouse symbols.',
  'Available comparator classes: project-defined MSigDB/manual composite pyroptosis and inflammasome sets; TCGA curated extended pyroptosis set; project-defined inflammation_noncore.',
  'No DOI-linked published weighted pyroptosis or inflammasome signature was found in the current project files; the curated 15-gene extended set is project-defined and must not be called an external published signature.',
  'Primary strict comparator removes all six PRII genes. Full signatures are sensitivity analyses only.',
  'Scores are recalculated independently within each cohort. VIF is 1/(1-R2) from regressing each score on all other scores.',
  paste0('TCGA primary tumor n=',ncol(expr),'; GSE79668 n=',ncol(gexpr),'; GSE292815 12-week mouse n=',nrow(m_comp)),
  'This audit is descriptive and does not prove functional pyroptosis, inflammasome activation, causal immunity, or independent prognosis.'
),file.path(outdir,'README.md'))
message('Wrote S10 signature audit to ',outdir)
