suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2); library(igraph)})
out <- "results/network_ligand_receptor"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
x <- read_csv(file.path(out, "lr_edge_summary.csv"), show_col_types=FALSE) %>%
  mutate(evidence_type="候选配体–受体表达支持")
stopifnot(nrow(x)>0)
x <- x %>% mutate(edge=paste(ligand, receptor, sep="→"), receiver=factor(receiver, levels=c("epithelial_tumor","stellate","t_cell")))
top <- x %>% filter(n_supported_samples>=3, prevalence>0) %>% arrange(desc(prevalence), desc(mean_score)) %>% slice_head(n=30)
write_csv(top, file.path(out, "lr_top30_edges.csv"))
p1 <- ggplot(top %>% mutate(edge=reorder(edge, prevalence)), aes(prevalence, edge, colour=receiver, size=mean_score)) + geom_point(alpha=.85) + scale_x_continuous(limits=c(0,1), labels=scales::percent) + labs(x="Supported-sample prevalence", y=NULL, colour="Receiver", size="Mean interaction score", title="Macrophage candidate ligand–receptor edges") + theme_bw(base_size=10) + theme(panel.grid.minor=element_blank())
ggsave(file.path(out, "lr_top_edges_dotplot.pdf"), p1, width=8.5, height=7)
ggsave(file.path(out, "lr_top_edges_dotplot.png"), p1, width=8.5, height=7, dpi=300)
hm <- x %>% filter(n_supported_samples>=3) %>% group_by(ligand, receiver) %>% summarise(mean_score=max(mean_score), prevalence=max(prevalence), .groups="drop")
p2 <- ggplot(hm, aes(receiver, reorder(ligand, mean_score), fill=mean_score)) + geom_tile(colour="white") + scale_fill_viridis_c() + labs(x=NULL, y="Macrophage ligand", fill="Mean score", title="Candidate macrophage-to-compartment signalling") + theme_bw(base_size=10) + theme(panel.grid=element_blank())
ggsave(file.path(out, "lr_ligand_receiver_heatmap.pdf"), p2, width=7.5, height=8)
ggsave(file.path(out, "lr_ligand_receiver_heatmap.png"), p2, width=7.5, height=8, dpi=300)
nodes <- bind_rows(
  tibble(name="C::macrophage", label="macrophage", type="sender"),
  top %>% transmute(name=paste0("C::", receiver), label=as.character(receiver), type="receiver"),
  top %>% transmute(name=paste0("L::", ligand), label=ligand, type="ligand"),
  top %>% transmute(name=paste0("R::", receptor), label=receptor, type="receptor")
) %>% distinct(name, .keep_all=TRUE)
edge_lr <- top %>% transmute(from=paste0("L::", ligand), to=paste0("R::", receptor), weight=mean_score, receiver=as.character(receiver), relation="ligand_receptor", evidence_type)
edge_send <- top %>% transmute(from="C::macrophage", to=paste0("L::", ligand), weight=mean_score, receiver=as.character(receiver), relation="candidate_sender", evidence_type)
edge_recv <- top %>% transmute(from=paste0("R::", receptor), to=paste0("C::", receiver), weight=mean_score, receiver=as.character(receiver), relation="candidate_receiver", evidence_type)
edges <- bind_rows(edge_send, edge_lr, edge_recv) %>% distinct(from,to,receiver,.keep_all=TRUE)
g <- graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
V(g)$color <- ifelse(V(g)$type=="ligand", "#D55E00", ifelse(V(g)$type=="receptor", "#0072B2", ifelse(V(g)$type=="receiver", "#009E73", "#666666")))
V(g)$size <- ifelse(V(g)$type %in% c("ligand","receptor"), 8, 13)
V(g)$label <- V(g)$label
E(g)$width <- 1 + 3*scales::rescale(E(g)$weight, to=c(0,1)); E(g)$arrow.size <- .25
write_graph(g, file.path(out, "lr_network_top30.graphml"), format="graphml")
pdf(file.path(out, "lr_network_top30.pdf"), width=11, height=8.5); set.seed(11); plot(g, layout=layout_with_fr(g), vertex.label.cex=.58, vertex.label.color="black", edge.curved=.12, main="Focused macrophage ligand–receptor candidate network"); dev.off()
png(file.path(out, "lr_network_top30.png"), width=1980, height=1530, res=180); set.seed(11); plot(g, layout=layout_with_fr(g), vertex.label.cex=.58, vertex.label.color="black", edge.curved=.12, main="Focused macrophage ligand–receptor candidate network"); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_plot.txt"))
