#!/usr/bin/env Rscript
# Figure S6: independent mouse metabolic-context cohorts.
# PRII is the manuscript name for legacy PSPS_core. Scores are cohort-specific
# mean-z scores; only within-cohort contrasts are displayed and no meta-analysis
# is performed.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork); library(grid)
})
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "graph/FigS6_PRII_external_metabolic_context.R"
root <- normalizePath(file.path(dirname(normalizePath(script_path)), ".."), mustWork = TRUE)
gd <- file.path(root, "graph"); vd <- file.path(root, "results", "external_validation")
out <- file.path(root, "results", "supplementary", "S6_external_metabolic_context")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read0 <- function(f) read_csv(file.path(vd, f), show_col_types = FALSE, progress = FALSE)
e234 <- read0("GSE234956_PSPS_sample_scores.csv")
e234_strain <- read0("GSE234956_PSPS_by_strain.csv")
e234_eff <- read0("GSE234956_HFHS_effects_adjusted_strain_sex.csv")
e190 <- read0("GSE190812_PSPS_sample_scores.csv"); e190_test <- read0("GSE190812_PSPS_group_tests.csv")
e249 <- read0("GSE249790_PSPS_sample_scores.csv"); e249_test <- read0("GSE249790_PSPS_group_tests.csv")
e293 <- read0("GSE293790_PSPS_sample_scores.csv"); e293_test <- read0("GSE293790_PSPS_group_tests.csv")
e293_gene <- read0("GSE293790_core_gene_group_means.csv")
stopifnot(nrow(e234)==45L, nrow(e190)==21L, nrow(e249)==7L, nrow(e293)==9L,
          all(c("PSPS_core","diet","strain","sex") %in% names(e234)),
          all(c("PSPS_core","group") %in% names(e190)), all(c("PSPS_core","group") %in% names(e249)),
          all(c("PSPS_core","group") %in% names(e293)))
ink <- "#243746"; teal <- "#2C7C7A"; blue <- "#4B7C9A"; coral <- "#CE6A58"; gold <- "#C08A3E"; muted <- "#6F7F89"; gridc <- "#E6ECEE"
group_cols <- c("10% fat + fiber"=blue,"44% fat + fiber"=coral,RD=blue,HFD=coral,HFD_O304=gold,CD=blue,DIO=coral,CDDIO=gold)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))

# Design/source tables ------------------------------------------------------
design <- tibble(
  cohort=c("GSE234956","GSE190812","GSE249790","GSE293790"),
  tissue_model=c("WT pancreatic islets","CBA pancreatic islet pooled libraries","sorted pancreatic beta cells","female KC/Kras-driven pancreas"),
  comparison=c("44% fat + fiber vs 10% fat + fiber","RD / HFD / HFD-O304","chow (CD) vs 2-day HFD","CD / DIO / DIO→CD"),
  n_per_group=c("23 vs 22","7 / 7 / 7 pools","4 vs 3","3 / 3 / 3"),
  sex_strain=c("C57BL/6J, CAST/EiJ, NZO/HlLtJ; female and male","CBA; male","C57BL/6N; not reported in local table","female KC; LSL-KrasG12D/+; p48Cre/+"),
  exposure=c("duration NR from public sample metadata","~10 weeks; pooled 5 mice/library","HFD for 2 days","11 kcal% fat CD or 60 kcal% fat DIO for 21 weeks; DIO→CD after 8 weeks DIO + 13 weeks CD"),
  score_scope="within-cohort six-gene mean-z PRII",
  interpretation=c("primary metabolic-context contrast","diet plus O304 intervention","early beta-cell metabolic stress","cancer-susceptible KC model"))
write_csv(design, file.path(out,"S6A_cohort_design_metadata.csv"))
e234_out <- e234 %>% transmute(cohort="GSE234956", sample=gsm, PRII=PSPS_core, group=diet, strain, sex, model="WT pancreatic islet")
e190_out <- e190 %>% transmute(cohort="GSE190812", sample, PRII=PSPS_core, group, strain="CBA", sex="male", model="pooled pancreatic islet library")
e249_out <- e249 %>% transmute(cohort="GSE249790", sample, PRII=PSPS_core, group, strain="C57BL/6N", sex="NR", model="sorted beta cell")
e293_out <- e293 %>% transmute(cohort="GSE293790", sample, PRII=PSPS_core, group, strain="KC/Kras", sex="female", model="whole pancreas, tumor-prone")
write_csv(bind_rows(e234_out,e190_out,e249_out,e293_out), file.path(out,"S6B_PRII_sample_scores_all_cohorts.csv"))

gene_order <- c("TXNIP","NLRP3","CASP1","GSDMD","IL1B","IL18")
g234 <- e234_eff %>% filter(feature %in% c("Txnip","Nlrp3","Casp1","Gsdmd","Il1b","Il18")) %>% transmute(cohort="GSE234956", contrast="44% fat + fiber − 10% fat + fiber", gene=toupper(feature), delta=effect_HFHS, p_value, q_value=FDR_BH, estimate_type="strain+sex-adjusted")
g190 <- e190_test %>% filter(feature!="PSPS_core") %>% transmute(cohort="GSE190812", contrast=recode(comparison,HFD_vs_RD="HFD − RD",HFD_O304_vs_HFD="HFD-O304 − HFD"), gene=toupper(feature), delta, p_value=mw_p, q_value=NA_real_, estimate_type="Mann–Whitney nominal")
g249 <- e249_test %>% filter(feature!="PSPS_core") %>% transmute(cohort="GSE249790", contrast="HFD − CD", gene=toupper(feature), delta, p_value=mw_p, q_value=FDR_BH, estimate_type="Mann–Whitney nominal")
g293 <- e293_gene %>% mutate(gene=toupper(gene)) %>% filter(gene %in% gene_order) %>% select(gene,group,mean_log2_norm) %>% pivot_wider(names_from=group, values_from=mean_log2_norm) %>% mutate(cohort="GSE293790", estimate_type="group mean difference") %>% transmute(cohort, contrast="DIO − CD", gene, delta=DIO-CD, p_value=NA_real_, q_value=NA_real_, estimate_type="group mean difference; no gene-level test")
gene_long <- bind_rows(g234,g190,g249,g293) %>% mutate(gene=factor(gene, levels=gene_order))
write_csv(gene_long, file.path(out,"S6C_core_gene_direction_long.csv"))

# Panel A: effect-direction summary. CI is only available for GSE234956.
e234_prii <- e234_eff %>% filter(feature=="PSPS_core")
e190_prii <- e190_test %>% filter(feature=="PSPS_core")
e249_prii <- e249_test %>% filter(feature=="PSPS_core")
e293_dio <- e293_test %>% filter(comparison=="DIO_vs_CD")
e293_switch <- e293_test %>% filter(comparison=="CDDIO_vs_CD")
forest <- tibble(
  label=c("GSE234956\n44% fat + fiber vs 10% fat + fiber","GSE190812\nHFD vs RD","GSE190812\nHFD-O304 vs HFD","GSE249790\nHFD vs CD (2 days)","GSE293790\nDIO vs CD (KC pancreas)","GSE293790\nDIO→CD vs CD (KC pancreas)"),
  effect=c(e234_prii$effect_HFHS,e190_prii$delta[e190_prii$comparison=="HFD_vs_RD"],e190_prii$delta[e190_prii$comparison=="HFD_O304_vs_HFD"],e249_prii$delta,e293_dio$delta,e293_switch$delta),
  lo=c(e234_prii$ci_low,NA,NA,NA,NA,NA), hi=c(e234_prii$ci_high,NA,NA,NA,NA,NA),
  p=c(e234_prii$p_value,e190_prii$mw_p[e190_prii$comparison=="HFD_vs_RD"],e190_prii$mw_p[e190_prii$comparison=="HFD_O304_vs_HFD"],e249_prii$mw_p,e293_dio$mw_p,e293_switch$mw_p),
  n=c(paste(sum(e234$diet=="44% fat + fiber"),"vs",sum(e234$diet=="10% fat + fiber")),"7 vs 7 pools","7 vs 7 pools",paste(e249_prii$n_HFD,"vs",e249_prii$n_CD),paste(e293_dio$n_a,"vs",e293_dio$n_b),paste(e293_switch$n_a,"vs",e293_switch$n_b)),
  q=c(e234_prii$FDR_BH,NA,NA,e249_prii$FDR_BH,NA,NA),
  estimate_type=c("strain+sex-adjusted","nominal MW","nominal MW","nominal MW; BH-adjusted q","nominal MW","nominal MW"),stringsAsFactors=FALSE) %>%
  mutate(label=factor(label,levels=rev(label)),direction=ifelse(effect>=0,"positive","negative"),stat=paste0("Δ=",sprintf("%+.3f",effect),"; P=",fmtp(p),ifelse(is.na(q),"",paste0("\nBH q=",fmtp(q))),"\n",n))
write_csv(forest %>% mutate(label=as.character(label)) %>% select(-stringsAsFactors), file.path(out,"S6A_within_cohort_effect_summary.csv"))
pA <- ggplot(forest,aes(y=label,x=effect,colour=direction)) + geom_vline(xintercept=0,linetype=2,linewidth=.4,colour="#AAB6BC") +
  geom_segment(aes(x=lo,xend=hi,yend=label),linewidth=.6,na.rm=TRUE) + geom_point(size=2.5,shape=21,fill="white",stroke=.7) +
  geom_text(aes(x=1.34,label=stat),hjust=0,size=2.55,colour=ink,lineheight=.92) + scale_colour_manual(values=c(positive=teal,negative=coral),guide="none") +
  scale_x_continuous(limits=c(-1.05,1.9),breaks=c(-1,-.5,0,.5,1,1.5),expand=c(0,0)) +
  labs(tag="A",title="Within-cohort PRII contrasts",subtitle="Effect estimates are not pooled; CI shown only where available",x="PRII difference (higher-fat/intervention − comparator)",y=NULL) +
  theme_classic(base_size=8.2) + theme(plot.title=element_text(face="bold",colour=ink,size=10),plot.subtitle=element_text(colour=muted,size=7.1,margin=margin(b=4)),axis.text.y=element_text(colour=ink,size=7),axis.text.x=element_text(colour=muted,size=7),axis.title.x=element_text(colour=muted,size=7.3,margin=margin(t=4)),axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major.x=element_line(colour=gridc,linewidth=.3),plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),plot.margin=margin(8,70,5,8))

# Panel B: GSE234956 multi-strain pancreatic islets.
strain_lab <- e234_strain %>% mutate(strain_display=paste0(strain,"\nunadj. Δ=",sprintf("%+.3f",delta),"; P=",fmtp(p_value),"\nn(44%/10%)=",n_HFHS,"/",n_control))
e234p <- e234 %>% left_join(strain_lab %>% select(strain,strain_display),by="strain") %>% mutate(diet=factor(diet,levels=c("10% fat + fiber","44% fat + fiber")),strain_display=factor(strain_display,levels=strain_lab$strain_display))
pB <- ggplot(e234p,aes(diet,PSPS_core,colour=diet)) + geom_boxplot(aes(group=diet),width=.5,outlier.shape=NA,alpha=.16,colour=ink,linewidth=.35) + geom_jitter(aes(shape=sex),width=.07,height=0,size=1.75,alpha=.82,stroke=.2) + facet_wrap(~strain_display,nrow=1) +
  scale_colour_manual(values=group_cols[c("10% fat + fiber","44% fat + fiber")],guide="none") + scale_shape_manual(values=c(Female=16,Male=17),name="Sex") +
  scale_y_continuous(limits=c(-1.75,2.7),breaks=c(-1.5,-1,-.5,0,.5,1,1.5,2,2.5),expand=c(0,0)) + labs(tag="B",title="GSE234956: multi-strain pancreatic islets",subtitle=paste0("Adjusted: Δ=",sprintf("%+.3f",e234_prii$effect_HFHS)," (95% CI ",sprintf("%+.3f",e234_prii$ci_low)," to ",sprintf("%+.3f",e234_prii$ci_high),"), P=",fmtp(e234_prii$p_value),", BH q=",fmtp(e234_prii$FDR_BH),"; adjusted for strain and sex. Duration NR."),x=NULL,y="Cohort-specific PRII") +
  theme_classic(base_size=8) + theme(plot.title=element_text(face="bold",colour=ink,size=9.8),plot.subtitle=element_text(colour=muted,size=6.8,margin=margin(b=4)),strip.background=element_rect(fill="#F1F4F5",colour=NA),strip.text=element_text(face="bold",colour=ink,size=5.8,lineheight=.9),axis.text=element_text(colour=muted,size=6.5),axis.title.y=element_text(colour=muted,size=7),axis.text.x=element_text(angle=25,hjust=1,size=6),axis.ticks.x=element_blank(),panel.grid.major.y=element_line(colour=gridc,linewidth=.3),panel.grid.major.x=element_blank(),legend.position="top",legend.title=element_text(size=6.8),legend.text=element_text(size=6.3),plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),plot.margin=margin(8,6,4,7))

# Panel C: GSE190812 pooled libraries; O304 is a pharmacologic arm.
e190p <- e190 %>% mutate(group=factor(group,levels=c("RD","HFD","HFD_O304")))
pC <- ggplot(e190p,aes(group,PSPS_core,colour=group)) + geom_boxplot(aes(group=group),width=.48,outlier.shape=NA,alpha=.16,colour=ink,linewidth=.35) + geom_jitter(width=.08,height=0,size=2,alpha=.85) + scale_colour_manual(values=group_cols[c("RD","HFD","HFD_O304")],labels=c("RD","HFD","HFD + O304"),guide="none") + scale_x_discrete(labels=c(RD="RD",HFD="HFD",HFD_O304="HFD + O304")) +
  annotate("text",x=2,y=1.2,label=paste0("HFD−RD: Δ=",sprintf("%+.3f",e190_prii$delta[e190_prii$comparison=="HFD_vs_RD"]),"; P=",fmtp(e190_prii$mw_p[e190_prii$comparison=="HFD_vs_RD"]),"\nHFD+O304−HFD: Δ=",sprintf("%+.3f",e190_prii$delta[e190_prii$comparison=="HFD_O304_vs_HFD"]),"; P=",fmtp(e190_prii$mw_p[e190_prii$comparison=="HFD_O304_vs_HFD"])),size=2.35,colour=ink,lineheight=.92) + scale_y_continuous(limits=c(-2.2,1.45),breaks=c(-2,-1.5,-1,-.5,0,.5,1),expand=c(0,0)) + labs(tag="C",title="GSE190812: pooled CBA islet libraries",subtitle="n=7/group; each library pooled islets from five mice; O304 is a pharmacologic arm",x=NULL,y="Cohort-specific PRII") +
  theme_classic(base_size=8) + theme(plot.title=element_text(face="bold",colour=ink,size=9.6),plot.subtitle=element_text(colour=muted,size=6.7,margin=margin(b=4)),axis.text=element_text(colour=muted,size=6.8),axis.title.y=element_text(colour=muted,size=7),axis.ticks.x=element_blank(),panel.grid.major.y=element_line(colour=gridc,linewidth=.3),plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),plot.margin=margin(8,7,4,7))

# Panel D: two contexts are separate facets; scores are not comparable across them.
e249p <- e249 %>% mutate(cohort_panel="GSE249790\nβ cells; 2-day HFD",group=factor(group,levels=c("CD","HFD")),group_label=factor(recode(group,CD="CD\nn=4",HFD="HFD\nn=3"),levels=c("CD\nn=4","HFD\nn=3")))
e293p <- e293 %>% mutate(cohort_panel="GSE293790\nKC pancreas; DIO model",group=factor(group,levels=c("CD","DIO","CDDIO")),group_label=factor(recode(group,CD="CD\nn=3",DIO="DIO\nn=3",CDDIO="DIO→CD\nn=3"),levels=c("CD\nn=3","DIO\nn=3","DIO→CD\nn=3")))
dd <- bind_rows(e249p %>% select(sample,PSPS_core,group,group_label,cohort_panel),e293p %>% select(sample,PSPS_core,group,group_label,cohort_panel)); dd$cohort_panel <- factor(dd$cohort_panel,levels=c("GSE249790\nβ cells; 2-day HFD","GSE293790\nKC pancreas; DIO model"))
pD_anno <- tibble(cohort_panel=factor(c("GSE249790\nβ cells; 2-day HFD","GSE293790\nKC pancreas; DIO model"),levels=levels(dd$cohort_panel)),
  x=c(1.5,2), y=c(1.31,1.31), label=c(paste0("HFD−CD: Δ=",sprintf("%+.3f",e249_prii$delta),"; P=",fmtp(e249_prii$mw_p),"\nBH q=",fmtp(e249_prii$FDR_BH)),paste0("DIO−CD: Δ=",sprintf("%+.3f",e293_dio$delta),"; P=",fmtp(e293_dio$mw_p),"\nDIO→CD−CD: Δ=",sprintf("%+.3f",e293_switch$delta),"; P=",fmtp(e293_switch$mw_p))))
pD <- ggplot(dd,aes(group_label,PSPS_core,colour=group)) + geom_boxplot(aes(group=group_label),width=.48,outlier.shape=NA,alpha=.16,colour=ink,linewidth=.35) + geom_jitter(width=.08,height=0,size=2,alpha=.85) + facet_wrap(~cohort_panel,scales="free_x",nrow=1) + scale_colour_manual(values=group_cols,guide="none") +
  geom_text(data=pD_anno,aes(x=x,y=y,label=label),inherit.aes=FALSE,size=2.3,colour=ink,lineheight=.92) +
  scale_y_continuous(limits=c(-1.15,1.52),breaks=c(-1,-.5,0,.5,1,1.5),expand=c(0,0)) + labs(tag="D",title="Context-specific score distributions",subtitle="GSE249790: sorted β cells; GSE293790: female KC/Kras-driven pancreas",x=NULL,y="Cohort-specific PRII") +
  theme_classic(base_size=8) + theme(plot.title=element_text(face="bold",colour=ink,size=9.6),plot.subtitle=element_text(colour=muted,size=6.7,margin=margin(b=4)),strip.background=element_rect(fill="#F1F4F5",colour=NA),strip.text=element_text(face="bold",colour=ink,size=7),axis.text=element_text(colour=muted,size=6.8),axis.title.y=element_text(colour=muted,size=7),axis.ticks.x=element_blank(),panel.grid.major.y=element_line(colour=gridc,linewidth=.3),plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),plot.margin=margin(8,7,4,7))

fig <- pA/(pB|pC)/pD + plot_layout(heights=c(1.32,1.05,1.02)) + plot_annotation(title="Context-dependent validation of PRII across independent metabolic cohorts",subtitle="Four independent mouse datasets; within-cohort contrasts only; absolute PRII values are not cross-cohort comparable",theme=theme(plot.title=element_text(face="bold",colour=ink,size=13,margin=margin(b=2)),plot.subtitle=element_text(colour=muted,size=8,margin=margin(b=6)),plot.margin=margin(8,10,5,10)))
pdf_file <- file.path(gd,"FigS6_PRII_external_metabolic_context.pdf"); tiff_file <- file.path(gd,"FigS6_PRII_external_metabolic_context.tiff"); png_file <- file.path(gd,"FigS6_PRII_external_metabolic_context_preview.png")
ggsave(pdf_file,fig,width=12,height=10.4,device=cairo_pdf,bg="white",limitsize=FALSE)
tiff(tiff_file,width=12,height=10.4,units="in",res=600,compression="lzw",type="cairo",bg="white"); print(fig); dev.off()
ggsave(png_file,fig,width=12,height=10.4,dpi=180,bg="white",limitsize=FALSE)
run_manifest <- c("Figure S6 run manifest",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),"Manuscript score name: PRII; source column PSPS_core is retained as a historical alias.","All scores are cohort-specific means of per-gene z scores and are not compared in absolute magnitude across cohorts.","GSE234956: 45 WT mouse islet samples; 44% fat + fiber versus 10% fat + fiber; model adjusted for strain and sex; 95% CI and BH q are from the frozen adjusted table; exposure duration was not available in public sample metadata.","GSE190812: CBA mouse islet pooled libraries; RD, HFD and HFD-O304, n=7 per arm; each library pooled islets from five mice; reported P values are nominal Mann–Whitney tests.","GSE249790: sorted beta-cell RNA-seq; CD n=4 versus 2-day HFD n=3; low/zero counts for several genes limit interpretation of negative directions; BH q is retained for the six-gene table.","GSE293790: female KC/Kras-driven pancreas; CD (11 kcal% fat), DIO (60 kcal% fat) and DIO→CD, n=3 per arm; cancer-susceptible model, not normal fatty pancreas.","No effect estimates were pooled and no cross-cohort meta-analysis was performed because tissues, models, exposures and score scales differ.","GSE190812, GSE249790 and GSE293790 P values are exploratory/nominal unless explicitly labelled otherwise.","Interpretation is limited to context-dependent transcriptional directionality; the figure does not establish functional pyroptosis, fatty pancreas or causality.","Outputs: 600-dpi LZW TIFF, vector PDF, preview PNG and panel-specific source CSV files.")
writeLines(run_manifest,file.path(out,"S6_run_manifest.txt")); writeLines(capture.output(sessionInfo()),file.path(out,"S6_sessionInfo.txt"))
message("Wrote ",pdf_file," and ",tiff_file)
