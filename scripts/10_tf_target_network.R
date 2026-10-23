#!/usr/bin/env Rscript
# Prespecified TF--target regulon activity and network analysis.
#
# This analysis uses the mouse/human DoRothEA prior regulons supplied by the
# dorothea package. TFs were prespecified from the study hypothesis (NF-kB,
# STAT, NRF2, ATF3, lipid-regulatory and hypoxia programs), rather than chosen
# after inspecting the data. PSPS core genes are excluded from activity scoring
# and from target edges to avoid arithmetic overlap with PSPS_core.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr);
  library(limma); library(ggplot2); library(igraph);
  library(SummarizedExperiment); library(dorothea); library(decoupleR)
})
options(stringsAsFactors = FALSE)
set.seed(292815)
root <- normalizePath('.', mustWork = TRUE)
outdir <- file.path(root, 'results/network_tf')
figdir <- file.path(outdir, 'figures')
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

tf_set <- c('Rela','Nfkb1','Stat1','Stat3','Nfe2l2','Atf3','Ppara','Pparg',
            'Srebf1','Srebf2','Cebpb','Hif1a')
tf_labels <- c(Rela='RELA/NF-kB', Nfkb1='NFKB1/NF-kB', Stat1='STAT1',
               Stat3='STAT3', Nfe2l2='NRF2', Atf3='ATF3', Ppara='PPARalpha',
               Pparg='PPARgamma', Srebf1='SREBP1', Srebf2='SREBP2',
               Cebpb='CEBPB', Hif1a='HIF1A')
core_mouse <- c('Txnip','Nlrp3','Casp1','Gsdmd','Il1b','Il18')
core_human <- toupper(core_mouse)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(rho=NA_real_, p=NA_real_, n=sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method='spearman', exact=FALSE))
  c(rho=unname(ct$estimate), p=ct$p.value, n=sum(ok))
}
partial_spearman <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 10L) return(c(rho=NA_real_, p=NA_real_, n=sum(ok)))
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  ex <- resid(lm(rx ~ rz)); ey <- resid(lm(ry ~ rz))
  ct <- suppressWarnings(cor.test(ex, ey, method='pearson'))
  c(rho=unname(ct$estimate), p=ct$p.value, n=sum(ok))
}

writeLines(c(
  paste0('TF set: ', paste(tf_set, collapse=', ')),
  'Prior regulon: dorothea_mm/dorothea_hs, confidence A/B/C.',
  paste0('PSPS core genes excluded from TF activity and network edges: ', paste(core_mouse, collapse=', ')),
  'Activity estimator: decoupleR run_mlm on log2(normalized expression + 1).',
  'Group model: limma with five pre-existing 12-week diet groups.',
  'Interpretation: TF activities and edges are computationally inferred associations; they do not prove TF binding or causal regulation.'
), file.path(outdir, 'analysis_manifest.txt'))

# ---- Mouse 12-week bulk ----------------------------------------------------
obj <- readRDS(file.path(root, 'results/qc/GSE292815_KC_12wk_matrix.rds'))
expr <- obj$expr; meta <- obj$metadata
stopifnot(ncol(expr) == nrow(meta), all(colnames(expr) == meta$matrix_column))
logexpr <- log2(expr + 1)
ps <- read_csv(file.path(root, 'results/12wk_primary/scores/12wk_sample_scores.csv'), show_col_types=FALSE)
ps <- ps[match(colnames(expr), ps$sample),,drop=FALSE]
stopifnot(all(ps$sample == colnames(expr)))

net0 <- dorothea_mm %>%
  filter(tf %in% tf_set, confidence %in% c('A','B','C')) %>%
  transmute(tf=as.character(tf), target=as.character(target), mor=as.numeric(mor),
            confidence=as.character(confidence)) %>% distinct()
net_audit <- net0 %>%
  group_by(tf) %>% summarise(n_regulon_total=n(), n_targets_matrix=sum(target %in% rownames(logexpr)),
                             n_core_targets=sum(target %in% core_mouse),
                             n_targets_used=sum(target %in% rownames(logexpr) & !(target %in% core_mouse)),
                             confidence=paste(sort(unique(confidence)), collapse=';'), .groups='drop')
write_csv(net_audit, file.path(outdir, 'mouse_regulon_audit.csv'))
net <- net0 %>% filter(target %in% rownames(logexpr), !(target %in% core_mouse)) %>%
  transmute(source=tf,target,mor)
if (!nrow(net)) stop('No mouse regulon targets overlap expression matrix after core exclusion.')

act_long <- decoupleR::run_mlm(logexpr, net, minsize=5) %>%
  transmute(tf=as.character(source), sample=as.character(condition), activity=as.numeric(score),
            p_value=as.numeric(p_value))
act <- act_long %>% select(tf,sample,activity) %>% pivot_wider(names_from=tf, values_from=activity)
act <- act %>% right_join(tibble(sample=colnames(expr)), by='sample') %>%
  left_join(meta %>% transmute(sample=matrix_column, group), by='sample') %>%
  left_join(ps %>% select(sample,PSPS_core), by='sample')
write_csv(act, file.path(outdir, 'mouse_tf_activity_scores.csv'))
write_csv(act_long, file.path(outdir, 'mouse_tf_activity_mlm_long.csv'))

tf_cor <- bind_rows(lapply(tf_set, function(tf) {
  z <- safe_cor(act$PSPS_core, act[[tf]])
  tibble(tf=tf, tf_label=unname(tf_labels[tf]), n=unname(z['n']), spearman_rho=unname(z['rho']), p_value=unname(z['p']))
})) %>% mutate(fdr_BH=p.adjust(p_value, 'BH'))
write_csv(tf_cor, file.path(outdir, 'mouse_tf_psps_spearman.csv'))

design <- model.matrix(~0+group, data=meta); colnames(design) <- sub('^group','',colnames(design))
fit <- lmFit(t(as.matrix(act[,tf_set,drop=FALSE])), design)
contrasts <- c(HO_Safflower_vs_Safflower='HO_Safflower_12wk-Safflower_12wk',
               HO_Sunflower_vs_Sunflower='HO_Sunflower_12wk-Sunflower_12wk')
contrast_tbl <- bind_rows(lapply(names(contrasts), function(nm) {
  cf <- eBayes(contrasts.fit(fit, makeContrasts(contrasts=contrasts[[nm]], levels=design)), robust=TRUE)
  topTable(cf, number=Inf, sort.by='none') %>% rownames_to_column('tf') %>%
    mutate(tf_label=unname(tf_labels[tf]), contrast=nm)
}))
write_csv(contrast_tbl, file.path(outdir, 'mouse_tf_group_contrasts.csv'))

# Add observed target-level associations to the prior regulon edge table.
cf <- eBayes(contrasts.fit(lmFit(logexpr, design), makeContrasts(contrasts=contrasts[['HO_Safflower_vs_Safflower']], levels=design)), robust=TRUE)
de_target <- topTable(cf, number=Inf, sort.by='none') %>% rownames_to_column('target') %>%
  transmute(target, target_logFC=logFC, target_adjP=adj.P.Val)
edge_tbl <- net0 %>% filter(target %in% rownames(logexpr), !(target %in% core_mouse)) %>%
  left_join(de_target, by='target') %>%
  rowwise() %>% mutate(target_tf_rho=unname(safe_cor(act[[tf]], logexpr[target,])[['rho']]),
                      target_tf_p=unname(safe_cor(act[[tf]], logexpr[target,])[['p']])) %>% ungroup() %>%
  mutate(tf_label=unname(tf_labels[tf]), target_present=TRUE, core_excluded=FALSE,
         edge_rank_score=abs(target_tf_rho) + 0.15*abs(target_logFC),
         evidence_type='表达关联', evidence_components='表达关联;数据库先验', prior_source='DoRothEA')
write_csv(edge_tbl, file.path(outdir, 'mouse_tf_target_edges.csv'))

top_edges <- edge_tbl %>% group_by(tf) %>% arrange(desc(edge_rank_score), .by_group=TRUE) %>% slice_head(n=12) %>% ungroup()
write_csv(top_edges, file.path(outdir, 'mouse_tf_target_edges_top12.csv'))
nodes <- bind_rows(tibble(name=unique(top_edges$tf), node_type='TF', label=unname(tf_labels[unique(top_edges$tf)])),
                   tibble(name=unique(top_edges$target), node_type='Target', label=unique(top_edges$target))) %>% distinct(name,.keep_all=TRUE)
write_csv(nodes, file.path(outdir, 'mouse_tf_target_nodes_top12.csv'))
g <- graph_from_data_frame(top_edges %>% transmute(from=tf,to=target,mor,evidence_type,evidence_components,prior_source), directed=TRUE, vertices=nodes)
V(g)$color <- ifelse(V(g)$node_type=='TF','#D95F02','#2C7FB8'); V(g)$size <- ifelse(V(g)$node_type=='TF',10,4)
V(g)$frame.color <- 'grey30'; V(g)$label <- ifelse(V(g)$node_type=='TF', V(g)$label, NA)
E(g)$color <- ifelse(E(g)$mor>0, '#1B9E77', '#7570B3'); E(g)$arrow.size <- 0.2; E(g)$width <- 0.8
write_graph(g, file.path(outdir, 'mouse_tf_target_network_top12.graphml'), format='graphml')
pdf(file.path(figdir, 'mouse_tf_target_network_top12.pdf'), width=13, height=10)
set.seed(292815); plot(g, layout=layout_with_fr(g, weights=NULL), vertex.label.cex=.75,
     main='Prespecified TF–target network (12-week mouse bulk; top 12 edges/TF)')
legend('topleft', legend=c('Prespecified TF','Target gene','Activation edge','Inhibition edge'),
       col=c('#D95F02','#2C7FB8','#1B9E77','#7570B3'), pch=c(19,19,15,15), pt.cex=c(1.2,1,1,1), bty='n')
dev.off()

# ---- TCGA-PAAD human external sensitivity analysis -------------------------
se <- readRDS(file.path(root, 'data/raw/TCGA_PAAD_RNAseq.rds'))
expr0 <- assay(se, 'tpm_unstrand'); rd <- as.data.frame(rowData(se)); gene <- toupper(as.character(rd$gene_name));
gene[is.na(gene)|gene==''] <- toupper(rownames(expr0)[is.na(gene)|gene=='']); gene <- sub('\\..*$', '', gene)
keep <- !duplicated(gene); exprh <- expr0[keep,,drop=FALSE]; rownames(exprh) <- gene[keep]
cd <- as.data.frame(colData(se)); tumor <- cd$sample_type=='Primary Tumor' | cd$shortLetterCode=='TP'; tumor[is.na(tumor)] <- FALSE
exprh <- log2(exprh[,tumor,drop=FALSE] + 1); cdh <- cd[colnames(exprh),,drop=FALSE]
tf_human <- toupper(tf_set); names(tf_human) <- tf_set
net_h0 <- dorothea_hs %>% filter(tf %in% tf_human, confidence %in% c('A','B','C')) %>%
  transmute(tf=as.character(tf), target=as.character(target), mor=as.numeric(mor), confidence=as.character(confidence)) %>% distinct()
net_h <- net_h0 %>% filter(target %in% rownames(exprh), !(target %in% core_human)) %>% select(tf,target,mor)
write_csv(net_h0 %>% group_by(tf) %>% summarise(n_regulon_total=n(), n_targets_matrix=sum(target%in%rownames(exprh)), n_core_targets=sum(target%in%core_human), n_targets_used=sum(target%in%rownames(exprh)&!(target%in%core_human)), confidence=paste(sort(unique(confidence)),collapse=';'),.groups='drop'), file.path(outdir, 'tcga_regulon_audit.csv'))
if (nrow(net_h)>0) {
  ah_long <- decoupleR::run_mlm(exprh, net_h %>% transmute(source=tf,target,mor), minsize=5) %>% transmute(tf=as.character(source), sample=as.character(condition), activity=as.numeric(score))
  ah <- ah_long %>% pivot_wider(names_from=tf, values_from=activity)
  # PSPS is recalculated from the same tumor expression after core exclusion only for TF target scoring.
  psps_h <- colMeans(t(scale(t(exprh[intersect(core_human, rownames(exprh)),,drop=FALSE]))), na.rm=TRUE)
  pur <- suppressWarnings(as.numeric(cdh$paper_ABSOLUTE.Purity)); names(pur) <- rownames(cdh)
  ah <- ah %>% mutate(PSPS_core=psps_h[match(sample, names(psps_h))], purity=pur[match(sample,names(pur))])
  # Shared TCGA TF sensitivity cohort: exact primary-tumor RNA-seq barcode
  # matching, finite PSPS, all prespecified TF activities and purity.
  tf_cols <- tf_human
  shared_tf <- ah %>% mutate(psps_complete=is.finite(PSPS_core),
                             tf_complete=if_all(all_of(tf_cols), is.finite),
                             purity_complete=is.finite(purity),
                             shared_complete=psps_complete & tf_complete & purity_complete)
  write_csv(shared_tf %>% select(sample,PSPS_core,purity,psps_complete,tf_complete,purity_complete,shared_complete), file.path(outdir, 'tcga_tf_shared_analysis_cohort.csv'))
  write_csv(tibble(cohort='TCGA-PAAD Primary Tumor',
                   metric=c('primary_tumor_total','PSPS_complete','TF_activities_complete','purity_complete','shared_complete'),
                   n=c(nrow(shared_tf),sum(shared_tf$psps_complete),sum(shared_tf$tf_complete),sum(shared_tf$purity_complete),sum(shared_tf$shared_complete)),
                   definition=c('Exact Primary Tumor RNA-seq barcode','Finite PSPS_core','Finite all prespecified TF activities','Finite paper_ABSOLUTE.Purity matched by exact barcode','PSPS + TF activities + purity complete')),
            file.path(outdir, 'tcga_tf_sample_count_audit.csv'))
  tcga_cor <- bind_rows(lapply(tf_human, function(h) {
    tf <- names(tf_human)[tf_human==h][1]; z <- partial_spearman(ah$PSPS_core, ah[[h]], ah$purity)
    tibble(tf_mouse=tf, tf_human=h, n=unname(z['n']), partial_spearman_rho=unname(z['rho']), p_value=unname(z['p']))
  })) %>% mutate(fdr_BH=p.adjust(p_value,'BH'))
  write_csv(ah, file.path(outdir, 'tcga_tf_activity_scores.csv')); write_csv(tcga_cor, file.path(outdir, 'tcga_tf_psps_partial_spearman_purity.csv'))
  p <- ggplot(tcga_cor, aes(x=reorder(tf_human, partial_spearman_rho), y=partial_spearman_rho, fill=fdr_BH<.05)) + geom_col() + coord_flip() + theme_classic(base_size=11) + scale_fill_manual(values=c('grey70','#D95F02')) + labs(x='TF', y='Partial Spearman rho (PSPS | purity)', title='TCGA-PAAD TF activity association') + theme(legend.position='none')
  ggsave(file.path(figdir, 'tcga_tf_psps_partial_spearman.pdf'), p, width=6.5, height=5.2)
} else { writeLines('No TCGA regulon targets overlapped the expression matrix after PSPS core exclusion.', file.path(outdir, 'tcga_tf_activity_scores.txt')) }

writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
readme <- c(
  'TF–target network analysis',
  '',
  paste0('Prespecified TFs: ', paste(tf_set, collapse=', ')),
  'Mouse prior network: DoRothEA confidence A/B/C; activity inferred by decoupleR run_mlm.',
  'PSPS core genes were audited and removed from TF activity target sets and network edges: Txnip, Nlrp3, Casp1, Gsdmd, Il1b, Il18.',
  'Mouse trend: RELA, STAT1 and NRF2 activity correlated positively with PSPS_core, whereas ATF3 correlated negatively; CEBPB/HIF1A were borderline after BH correction.',
  'Matched HO-Safflower contrast: ATF3 decreased significantly after BH correction; STAT1 and HIF1A increased nominally but were not BH-significant. HO-Sunflower showed no BH-significant TF differences.',
  'TCGA-PAAD purity-adjusted sensitivity: NFKB1, STAT1 and PPARG positive; STAT3 and ATF3 negative associations with PSPS_core; report as external computational sensitivity rather than causal validation.',
  'Network edge evidence_type: 表达关联; evidence_components=表达关联;数据库先验; prior_source=DoRothEA. Edges combine prior TF–target relationships with observed target logFC and TF–target activity correlation; they do not prove binding or causal regulation.',
  'Key files: mouse_tf_activity_scores.csv; mouse_tf_psps_spearman.csv; mouse_tf_group_contrasts.csv; mouse_tf_target_edges.csv; mouse_tf_target_nodes_top12.csv; figures/mouse_tf_target_network_top12.pdf.'
)
# Preserve a manually curated README (including cohort/evidence wording) when it exists.
if (!file.exists(file.path(outdir, 'README.md'))) writeLines(readme, file.path(outdir, 'README.md'))
message('Completed TF network analysis: ', outdir)
