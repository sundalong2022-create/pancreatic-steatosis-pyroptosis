#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(readr); library(dplyr); library(igraph)})
root <- normalizePath(getwd())
out <- file.path(root, "results/network_integrated")
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)

tf <- read_csv(file.path(root, "results/network_tf/mouse_tf_target_edges_top12.csv"), show_col_types=FALSE) %>%
  arrange(edge_rank_score) %>% group_by(tf) %>% slice_head(n=8) %>% ungroup() %>%
  transmute(from=paste0("TF_", toupper(tf)), to=toupper(target), edge_type="TF-target (DoRothEA + expression)", evidence_type="表达关联", evidence_components="表达关联;数据库先验", prior_source="DoRothEA", weight=pmax(edge_rank_score,0.1))
ppi <- read_csv(file.path(root, "results/network_ppi/ppi_edges.csv"), show_col_types=FALSE)
psps <- c("TXNIP","NLRP3","CASP1","GSDMD","IL1B","IL18")
ppi <- ppi %>% mutate(from=toupper(from), to=toupper(to), edge_type="STRING PPI", evidence_type="数据库先验", evidence_components="数据库先验", prior_source="STRING", weight=combined_score)
ppi_keep <- bind_rows(ppi %>% filter(from %in% psps | to %in% psps) %>% arrange(desc(combined_score)) %>% slice_head(n=120),
                      ppi %>% arrange(desc(combined_score)) %>% slice_head(n=60)) %>% distinct(from,to,.keep_all=TRUE) %>% select(from,to,edge_type,evidence_type,evidence_components,prior_source,weight)
lr <- read_csv(file.path(root, "results/network_ligand_receptor/lr_top_edges.csv"), show_col_types=FALSE) %>%
  arrange(desc(prevalence), desc(mean_score)) %>% slice_head(n=25) %>%
  transmute(from=toupper(ligand), to=toupper(receptor), edge_type="scRNA LR expression support", evidence_type="候选配体–受体表达支持", evidence_components="候选配体–受体表达支持", prior_source="OmniPath + sample expression", weight=pmax(mean_score,0.1))
edges <- bind_rows(tf, ppi_keep, lr) %>% distinct(from,to,edge_type,.keep_all=TRUE)
nodes <- tibble(name=unique(c(edges$from, edges$to))) %>% mutate(is_tf=grepl("^TF_",name), symbol=sub("^TF_","",name), category=case_when(symbol %in% psps ~ "PSPS core", is_tf ~ "Prespecified TF", name %in% unique(c(lr$from,lr$to)) ~ "LR ligand/receptor", TRUE ~ "Other gene"))
g <- graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
cols <- c("PSPS core"="#D73027","Prespecified TF"="#4575B4","LR ligand/receptor"="#1A9850","Other gene"="#BDBDBD")
V(g)$color <- cols[V(g)$category]; V(g)$size <- ifelse(V(g)$category=="PSPS core",8,ifelse(V(g)$category=="Prespecified TF",7,4)); V(g)$label <- ifelse(V(g)$is_tf,V(g)$symbol,ifelse(V(g)$category %in% c("PSPS core","LR ligand/receptor"),V(g)$symbol,NA)); V(g)$frame.color <- "white"; V(g)$label.color <- "#222222"; V(g)$label.cex <- .55
E(g)$color <- c("TF-target (DoRothEA + expression)"="#4575B4","STRING PPI"="#969696","scRNA LR expression support"="#1A9850")[E(g)$edge_type]; E(g)$width <- ifelse(E(g)$edge_type=="STRING PPI",.7,1.5); E(g)$arrow.size <- .25; E(g)$curved <- .08
set.seed(13); lay <- layout_with_fr(g, niter=1500)
draw <- function(file, type="pdf") { if(type=="pdf") pdf(file,width=13,height=10,useDingbats=FALSE) else png(file,width=2600,height=2000,res=220); par(mar=c(0,0,4,0)); plot(g,layout=lay,main="Integrated candidate network: TF regulation, scRNA communication and STRING PPI",sub=paste0(vcount(g)," nodes; ",ecount(g)," edges; evidence layers are color-coded"),vertex.label.dist=.4,edge.arrow.size=.25); legend("topleft",legend=c("PSPS core","Prespecified TF","LR ligand/receptor","Other gene","TF-target","STRING PPI","scRNA LR"),col=c("#D73027","#4575B4","#1A9850","#BDBDBD","#4575B4","#969696","#1A9850"),pch=c(19,19,19,19,NA,NA,NA),lty=c(NA,NA,NA,NA,1,1,1),lwd=c(NA,NA,NA,NA,2,2,2),bty="n",cex=.8); dev.off() }
draw(file.path(out,"figures/integrated_network_overview.pdf"),"pdf"); draw(file.path(out,"figures/integrated_network_overview.png"),"png")
write_csv(as_data_frame(g,what="vertices"),file.path(out,"integrated_nodes.csv")); write_csv(as_data_frame(g,what="edges"),file.path(out,"integrated_edges.csv"))
readme_lines <- c("Integrated network figure",paste0("nodes=",vcount(g)),paste0("edges=",ecount(g)),paste0("TF edges=",sum(edges$edge_type=="TF-target (DoRothEA + expression)")),paste0("PPI edges=",sum(edges$edge_type=="STRING PPI")),paste0("LR edges=",sum(edges$edge_type=="scRNA LR expression support")),"Evidence labels: TF edges = evidence_type 表达关联; evidence_components=表达关联;数据库先验; prior_source=DoRothEA. PPI edges = 数据库先验; LR edges = 候选配体–受体表达支持.","Edge widths are not comparable across evidence layers. All edges are hypothesis-generating and do not establish physical binding, secretion, receptor activation, or causal regulation.")
if (!file.exists(file.path(out, "README.md"))) writeLines(readme_lines,file.path(out,"README.md"))
cat("Completed integrated network: ",out,"\n")
