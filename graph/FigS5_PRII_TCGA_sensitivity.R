#!/usr/bin/env Rscript
# Figure S5: TCGA-PAAD deconvolution and sensitivity analyses.
# New outputs use PRII; PSPS_core is the legacy source-column name.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(survival); library(grid)
})
args <- commandArgs(trailingOnly = FALSE)
fa <- grep("^--file=", args, value = TRUE)
script <- if (length(fa)) sub("^--file=", "", fa[1]) else "graph/FigS5_PRII_TCGA_sensitivity.R"
root <- normalizePath(file.path(dirname(normalizePath(script)), ".."), mustWork = TRUE)
gd <- file.path(root, "graph")
out <- file.path(root, "results/supplementary/S5_tcga_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dd <- file.path(root, "results/tcga_paad_external/deconvolution")
metrics <- read_csv(file.path(dd, "deconvolution_harmonized_metrics.csv"), show_col_types = FALSE)
assocs <- read_csv(file.path(dd, "deconvolution_psps_associations.csv"), show_col_types = FALSE)
focus <- read_csv(file.path(dd, "deconvolution_focus_associations.csv"), show_col_types = FALSE)
cons <- read_csv(file.path(dd, "deconvolution_method_consistency.csv"), show_col_types = FALSE)
epic <- read_csv(file.path(dd, "epic_convergence_sensitivity.csv"), show_col_types = FALSE)
status <- read_csv(file.path(dd, "deconvolution_method_status.csv"), show_col_types = FALSE)
gof <- read_csv(file.path(dd, "epic_fit_gof.csv"), show_col_types = FALSE)
surv0 <- read_csv(file.path(root, "results/tcga_paad_external/survival_input.csv"), show_col_types = FALSE)
stopifnot(nrow(metrics) == 178L, sum(gof$convergeCode == 0, na.rm=TRUE) == 167L)
stopifnot(nrow(surv0) == 160L)

ink <- "#243746"; teal <- "#327C7A"; blue <- "#4C7E9C"
coral <- "#D16B58"; muted <- "#70818D"; gridc <- "#E5EBEE"
fmtq <- function(q) ifelse(q < 0.001, formatC(q, format="e", digits=1), sprintf("%.3f", q))

# A. Purity-adjusted association forest: types are intentionally labeled by method.
labs <- c(xCell_Macrophages="Macrophages (xCell enrichment)",
          EPIC_Macrophages="Macrophages (EPIC fraction)",
          CIBERSORT_M0="LM22-compatible M0", CIBERSORT_M1="LM22-compatible M1",
          CIBERSORT_M2="LM22-compatible M2", ImmuneScore="ImmuneScore (ESTIMATE)",
          StromalScore="StromalScore (ESTIMATE)")
a <- focus %>% filter(analysis=="partial_spearman_adjusted_for_ABSOLUTE_purity",
                      metric %in% names(labs)) %>%
  mutate(label=unname(labs[metric]), label=factor(label, levels=rev(unname(labs))),
         stat=paste0("n=",n,"; BH q=",fmtq(fdr_BH)),
         sig=ifelse(fdr_BH<.05,"BH q < 0.05","BH q ≥ 0.05"))
stopifnot(nrow(a)==length(labs))
write_csv(a %>% transmute(metric, metric_label=label, n, partial_rho=rho,
                          p_value, df, p_value_method, BH_q=fdr_BH, covariate, analysis),
          file.path(out,"S5A_PRII_deconvolution_purity_adjusted_associations.csv"))
pA <- ggplot(a,aes(y=label))+
  geom_vline(xintercept=0,linetype=2,linewidth=.4,colour="#AAB6BC")+
  geom_segment(aes(x=0,xend=rho,yend=label,colour=sig),linewidth=.6)+
  geom_point(aes(x=rho,fill=sig),shape=21,size=2.4,colour=ink,stroke=.3)+
  geom_text(aes(x=.78,label=stat),hjust=0,size=2.5,colour=ink)+
  annotate("text",x=.78,y=length(labs)+.55,label="n; BH q",hjust=0,
           size=2.6,colour=muted,fontface="bold")+
  scale_colour_manual(values=c("BH q < 0.05"=teal,"BH q ≥ 0.05"="#A8B3B8"),guide="none")+
  scale_fill_manual(values=c("BH q < 0.05"=teal,"BH q ≥ 0.05"="white"),guide="none")+
  scale_x_continuous(limits=c(-.28,1.35),breaks=c(-.2,0,.2,.4,.6),expand=c(0,0))+
  labs(tag="A",title="PRII associations after ABSOLUTE-purity adjustment",
       subtitle="Rank-residual partial correlations; purity-complete cases",
       x="Partial Spearman ρ",y=NULL)+
  theme_classic(base_size=8.4)+
  theme(plot.title=element_text(face="bold",colour=ink,size=10),
        plot.subtitle=element_text(colour=muted,size=7.2,margin=margin(b=4)),
        axis.text.y=element_text(colour=ink,size=7.4),axis.text.x=element_text(colour=muted,size=7),
        axis.title.x=element_text(colour=muted,size=7.4,margin=margin(t=4)),
        axis.line.y=element_blank(),axis.ticks.y=element_blank(),
        panel.grid.major.x=element_line(colour=gridc,linewidth=.3),
        plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),
        plot.margin=margin(8,9,4,7))

# B. Method correlation matrix for five macrophage-related estimates.
mm <- c("xCell_Macrophages","EPIC_Macrophages","CIBERSORT_M0","CIBERSORT_M1","CIBERSORT_M2")
ms <- c(xCell_Macrophages="xCell",EPIC_Macrophages="EPIC",
        CIBERSORT_M0="LM22 M0",CIBERSORT_M1="LM22 M1",CIBERSORT_M2="LM22 M2")
pairs <- cons %>% filter(concept=="macrophage",metric_1%in%mm,metric_2%in%mm)
stopifnot(nrow(pairs)==10L)
b <- bind_rows(pairs %>% transmute(m1=metric_1,m2=metric_2,n,rho,q=fdr_BH),
               pairs %>% transmute(m1=metric_2,m2=metric_1,n,rho,q=fdr_BH),
               tibble(m1=mm,m2=mm,n=178L,rho=1,q=NA_real_)) %>%
  mutate(m1=factor(m1,levels=mm),m2=factor(m2,levels=rev(mm)),
         txt=ifelse(is.na(q),"1.00",paste0(sprintf("%.2f",rho),ifelse(q<.05,"*",""))))
write_csv(pairs %>% mutate(metric_1_label=unname(ms[metric_1]),
                           metric_2_label=unname(ms[metric_2])) %>%
            select(concept,metric_1,metric_1_label,metric_2,metric_2_label,n,
                   spearman_rho=rho,p_value,BH_q=fdr_BH),
          file.path(out,"S5B_macrophage_method_pairwise_correlations.csv"))
pB <- ggplot(b,aes(m1,m2,fill=rho))+geom_tile(colour="white",linewidth=.55)+
  geom_text(aes(label=txt),size=2.55,colour=ink)+
  scale_x_discrete(labels=ms,position="top")+scale_y_discrete(labels=ms)+
  scale_fill_gradient2(low="#B85C54",mid="#F6F7F7",high="#387D89",
                       midpoint=0,limits=c(-1,1),name="Spearman ρ")+
  labs(tag="B",title="Macrophage-related estimates across methods",
       subtitle="Pairwise correlations; * BH q < 0.05 across 35 comparisons")+
  theme_minimal(base_size=7.8)+theme(plot.title=element_text(face="bold",colour=ink,size=9.2),
    plot.subtitle=element_text(colour=muted,size=6.6,margin=margin(b=3)),
    axis.title=element_blank(),axis.text=element_text(colour=ink,size=6.5),
    panel.grid=element_blank(),legend.position="bottom",legend.title=element_text(size=6.4),
    legend.text=element_text(size=6.1),legend.key.width=unit(13,"mm"),
    plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),
    plot.margin=margin(8,5,4,7))

# C. EPIC macrophage estimate before and after removing non-converged fits.
c1 <- focus %>% filter(metric=="EPIC_Macrophages",
  analysis=="partial_spearman_adjusted_for_ABSOLUTE_purity") %>%
  transmute(dataset="All EPIC estimates",n,rho,p_value,df,p_value_method,BH_q=fdr_BH)
c2 <- epic %>% filter(metric=="EPIC_Macrophages",
  analysis=="EPIC_converged_partial_ABSOLUTE_purity") %>%
  transmute(dataset="EPIC converged only",n,rho,p_value,df,p_value_method,BH_q=fdr_BH)
c <- bind_rows(c1,c2) %>% mutate(dataset=factor(dataset,levels=rev(c("All EPIC estimates","EPIC converged only"))),
                                  txt=paste0("n=",n,"; BH q=",fmtq(BH_q)))
write_csv(c,file.path(out,"S5C_EPIC_macrophage_convergence_sensitivity.csv"))
pC <- ggplot(c,aes(y=dataset))+geom_vline(xintercept=0,linetype=2,linewidth=.4,colour="#AAB6BC")+
  geom_point(aes(x=rho,fill=dataset),shape=21,size=2.7,colour=ink,stroke=.3)+
  geom_text(aes(x=.62,label=txt),hjust=0,size=2.5,colour=ink)+
  scale_fill_manual(values=c("All EPIC estimates"=blue,"EPIC converged only"=coral),guide="none")+
  scale_x_continuous(limits=c(0,1.05),breaks=c(0,.25,.5,.75,1),expand=c(0,0))+
  labs(tag="C",title="EPIC macrophage convergence sensitivity",
       subtitle="11/178 samples did not converge; n reflects purity completeness",
       x="Partial Spearman ρ",y=NULL)+
  theme_classic(base_size=8)+theme(plot.title=element_text(face="bold",colour=ink,size=9.1),
    plot.subtitle=element_text(colour=muted,size=6.5,margin=margin(b=3)),
    axis.text=element_text(colour=muted,size=6.7),axis.text.y=element_text(colour=ink,size=7),
    axis.title.x=element_text(colour=muted,size=6.8,margin=margin(t=3)),
    axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    panel.grid.major.x=element_line(colour=gridc,linewidth=.3),
    plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),
    plot.margin=margin(8,7,4,7))

# D. Distinct purity estimators compared directly.
pd <- metrics %>% transmute(sample_id,ABSOLUTE=as.numeric(ABSOLUTE_Purity),
                            ESTIMATE=as.numeric(TumorPurity)) %>%
  filter(is.finite(ABSOLUTE),is.finite(ESTIMATE))
pt <- suppressWarnings(cor.test(pd$ABSOLUTE,pd$ESTIMATE,method="spearman",exact=FALSE))
ps <- tibble(n=nrow(pd),spearman_rho=unname(pt$estimate),p_value=pt$p.value)
write_csv(pd,file.path(out,"S5D_purity_estimator_sample_data.csv"))
write_csv(ps,file.path(out,"S5D_purity_estimator_correlation.csv"))
pD <- ggplot(pd,aes(ABSOLUTE,ESTIMATE))+geom_point(shape=21,size=1.15,alpha=.65,
  fill=blue,colour="white",stroke=.18)+
  geom_smooth(method="lm",formula=y~x,se=TRUE,colour=coral,fill="#F1D5CF",linewidth=.6)+
  annotate("label",x=.06,y=.98,label=sprintf("Spearman ρ=%.3f\nn=%d; P<0.001",
    ps$spearman_rho,ps$n),hjust=0,vjust=1,size=2.45,colour=ink,fill="white",linewidth=.18)+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=c(0,0))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=c(0,0))+
  labs(tag="D",title="Correlation between purity estimators",subtitle="Distinct measures; not interchangeable",
       x="ABSOLUTE tumor purity",y="ESTIMATE TumorPurity")+
  theme_classic(base_size=7.8)+theme(plot.title=element_text(face="bold",colour=ink,size=9.1),
    plot.subtitle=element_text(colour=muted,size=6.5,margin=margin(b=3)),
    axis.text=element_text(colour=muted,size=6.5),axis.title=element_text(colour=muted,size=6.8),
    plot.tag=element_text(face="bold",colour=ink,size=12),plot.tag.position=c(0,1),
    plot.margin=margin(8,6,4,7))

# E. Median-split KM and risk table are descriptive; continuous Cox is primary.
sd <- surv0 %>% transmute(sample_id,PRII=as.numeric(PSPS_core),
                          time_days=as.numeric(time_days),event=as.integer(event)) %>%
  filter(is.finite(PRII),is.finite(time_days),!is.na(event))
med <- median(sd$PRII)
sd <- sd %>% mutate(group=factor(ifelse(PRII>=med,"High","Low"),levels=c("Low","High")))
fit <- survfit(Surv(time_days,event)~group,data=sd)
ss <- summary(fit,censored=TRUE)
curve <- tibble(year=ss$time/365.25,surv=ss$surv,
                group=factor(sub("^group=","",as.character(ss$strata)),levels=c("Low","High")),
                n_censor=ss$n.censor)
curve <- bind_rows(tibble(year=0,surv=1,group=factor(c("Low","High"),levels=c("Low","High")),n_censor=0L),curve) %>%
  arrange(group,year)
maxyear <- min(5,ceiling(max(sd$time_days)/365.25)); times <- 0:maxyear
curve <- filter(curve,year<=maxyear)
risk <- expand_grid(group=factor(c("Low","High"),levels=c("Low","High")),year=times) %>%
  rowwise() %>% mutate(at_risk=sum(sd$group==group & sd$time_days>=year*365.25)) %>% ungroup()
gs <- sd %>% group_by(group) %>% summarise(n=n(),deaths=sum(event==1),median_PRII=median(PRII),.groups="drop")
lr <- survdiff(Surv(time_days,event)~group,data=sd)
lp <- pchisq(lr$chisq,df=length(lr$n)-1,lower.tail=FALSE)
write_csv(sd,file.path(out,"S5E_KM_sample_data.csv"))
write_csv(risk,file.path(out,"S5E_KM_risk_table.csv"))
write_csv(gs %>% mutate(logrank_p=lp),file.path(out,"S5E_KM_group_summary.csv"))
deaths_hi <- gs$deaths[gs$group=="High"]; n_hi <- gs$n[gs$group=="High"]
deaths_lo <- gs$deaths[gs$group=="Low"]; n_lo <- gs$n[gs$group=="Low"]
subE <- paste0("Descriptive median split; log-rank P=",sprintf("%.3f",lp),
               "; deaths ",deaths_hi,"/",n_hi," vs ",deaths_lo,"/",n_lo)
pEcurve <- ggplot(curve,aes(year,surv,colour=group))+geom_step(linewidth=.68,direction="hv")+
  geom_point(data=filter(curve,n_censor>0),shape=3,size=1.1,stroke=.5,show.legend=FALSE)+
  scale_colour_manual(values=c("Low"=blue,"High"=coral),name="PRII group")+
  scale_x_continuous(limits=c(0,maxyear),breaks=times,labels=rep("",length(times)),expand=c(.01,0))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=c(0,0))+
  labs(tag="E",title="Overall survival by median-split PRII",subtitle=subE,x=NULL,y="Survival")+
  theme_classic(base_size=7.7)+theme(plot.title=element_text(face="bold",colour=ink,size=9.1),
    plot.subtitle=element_text(colour=muted,size=6.5,margin=margin(b=3)),
    axis.text=element_text(colour=muted,size=6.5),axis.title.y=element_text(colour=muted,size=6.7),
    legend.position="top",legend.title=element_text(size=6.1),legend.text=element_text(size=6.1),
    legend.key.width=unit(9,"mm"),plot.tag=element_text(face="bold",colour=ink,size=12),
    plot.tag.position=c(0,1),plot.margin=margin(8,6,0,7))
pErisk <- ggplot(risk,aes(year,group,label=at_risk,colour=group))+geom_text(size=2.4,show.legend=FALSE)+
  scale_colour_manual(values=c("Low"=blue,"High"=coral),guide="none")+
  scale_x_continuous(limits=c(0,maxyear),breaks=times,labels=times,expand=c(.01,0))+
  scale_y_discrete(limits=rev(c("Low","High")),labels=c("High","Low"))+
  labs(x="Follow-up (years)",y="At risk")+
  theme_classic(base_size=7)+theme(axis.text=element_text(colour=muted,size=6.2),
    axis.text.y=element_text(colour=ink),axis.title=element_text(colour=muted,size=6.4),
    axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(0,6,4,7))
pE <- pEcurve/pErisk+plot_layout(heights=c(3.1,.8))

fig <- pA/(pB|pC)/(pD|pE)+plot_layout(heights=c(1.15,1,1.12))+
  plot_annotation(title="TCGA-PAAD sensitivity analyses for PRII and immune-context estimates",
    subtitle="178 primary tumors; purity-adjusted analyses n=149; EPIC converged-only sensitivity n=142",
    theme=theme(plot.title=element_text(face="bold",colour=ink,size=12),
      plot.subtitle=element_text(colour=muted,size=8),
      plot.margin=margin(8,10,5,10)))
pdf <- file.path(gd,"FigS5_PRII_TCGA_sensitivity.pdf")
tiffout <- file.path(gd,"FigS5_PRII_TCGA_sensitivity.tiff")
png <- file.path(gd,"FigS5_PRII_TCGA_sensitivity_preview.png")
ggsave(pdf,fig,width=11.2,height=10.3,device=cairo_pdf,bg="white",limitsize=FALSE)
tiff(tiffout,width=11.2,height=10.3,units="in",res=600,compression="lzw",type="cairo",bg="white")
print(fig); dev.off()
ggsave(png,fig,width=11.2,height=10.3,dpi=180,bg="white",limitsize=FALSE)
write_csv(status,file.path(out,"S5_method_status.csv"))
write_csv(assocs %>% mutate(PRII_source_column="PSPS_core"),
          file.path(out,"S5_all_deconvolution_associations_source.csv"))
writeLines(capture.output(sessionInfo()),file.path(out,"S5_sessionInfo.txt"))
writeLines(c(
  "Figure S5 run manifest",
  paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),
  "TCGA-PAAD exact-barcode Primary Tumor cohort: 178 samples.",
  "PRII is the manuscript name for the legacy PSPS_core six-gene score.",
  "xCell outputs are enrichment scores; EPIC outputs are estimated fractions; ESTIMATE outputs are transcriptomic scores/estimates.",
  "LM22-compatible estimates are from local SVR (permutations=0); official CIBERSORTx was not run.",
  "Purity-adjusted analyses use paper_ABSOLUTE.Purity; complete cases are metric-specific, generally n=149.",
  "Partial-correlation P values use the rank-residual t approximation with df=n-3 for one covariate; BH families are inherited from the full metric tables.",
  "Method-to-method macrophage correlations use the original 35-comparison BH family.",
  "EPIC converged in 167/178 samples; 11 non-converged samples were excluded in the sensitivity subset.",
  "Panel E median-split KM is descriptive; continuous PRII Cox is the primary survival analysis.",
  "Purity estimates from ABSOLUTE and ESTIMATE are correlated but not interchangeable.",
  "PDF and 600-dpi LZW TIFF plus panel source CSVs are provided."
),file.path(out,"S5_run_manifest.txt"))
cat("Wrote ",pdf," and ",tiffout,"; n=",nrow(sd),"; log-rank P=",signif(lp,3),"\n",sep="")
