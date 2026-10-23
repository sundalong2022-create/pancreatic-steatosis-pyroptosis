# ============================================================================
# 01_build_gene_sets.R：构建焦亡、炎症小体、脂质代谢和氧化应激/ROS 基因集
# ============================================================================
# 使用方法：在项目根目录运行 source("scripts/01_build_gene_sets.R")。
# 默认物种为人类；如需小鼠，将 organism 改为 "mouse"，脚本会使用小鼠符号。
# 固定核心轴：TXNIP–NLRP3–CASP1–GSDMD–IL1B/IL18。
# TXNIP 是 ROS 与 NLRP3 激活之间的耦联节点；CASP8/CASP3-GSDME 是替代性模块。

# ------------------------- 1. 项目参数 --------------------------------------
organism <- "human"  # 可改为 "mouse"
out_dir <- file.path("results", "gene_sets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- 2. 依赖检查 ----------------------------------------
required_packages <- c("msigdbr", "dplyr", "purrr", "stringr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "\n请先安装缺失的 R 包后再运行脚本。")
}
suppressPackageStartupMessages({
  library(msigdbr); library(dplyr); library(purrr)
  library(stringr); library(readr); library(tibble)
})

# ------------------------- 3. 物种和核心轴 ------------------------------------
if (organism == "human") {
  msig_species <- "Homo sapiens"
  core_axis <- c("TXNIP", "NLRP3", "CASP1", "GSDMD", "IL1B", "IL18")
} else if (organism == "mouse") {
  msig_species <- "Mus musculus"
  core_axis <- c("Txnip", "Nlrp3", "Casp1", "Gsdmd", "Il1b", "Il18")
} else stop("organism 只能是 human 或 mouse")
message("物种：", organism, "；核心轴：", paste(core_axis, collapse = " – "))

# ------------------------- 4. 手工锚定基因 ------------------------------------
# 手工集合用于稳定保留主轴。CASP11 主要是小鼠基因；人类使用 CASP4/CASP5。
# CASP8/CASP3/GSDME 作为替代性焦亡分支记录，不与 CASP1-GSDMD 主轴混淆。
manual_sets <- if (organism == "human") list(
  pyroptosis = c(core_axis, "AIM2","PYCARD","CASP4","CASP5","NLRC4","GBP1","GBP2","GBP4","GBP5","GSDME","CASP8","CASP3"),
  inflammasome = c(core_axis, "AIM2","NLRC4","NLRP1","MEFV","P2RX7","NEK7","CARD8","PYCARD"),
  lipid_metabolism = c("SREBF1","SREBF2","FASN","ACACA","ACACB","SCD","CPT1A","CPT2","PPARA","PPARG","LPL","CD36","HADHA","HADHB","ABCA1","ABCG1","DGAT1","DGAT2","PNPLA2","LIPE"),
  oxidative_stress_ROS = c("TXNIP","NFE2L2","KEAP1","SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX3","TXN","TXNRD1","NOX1","NOX2","CYBB","DUOX1","DUOX2","HMOX1","GCLC","GCLM")
) else list(
  pyroptosis = c(core_axis, "Aim2","Pycard","Casp11","Nlrc4","Gsdme","Casp8","Casp3"),
  inflammasome = c(core_axis, "Aim2","Nlrc4","Nlrp1","Mefv","P2rx7","Nek7","Card8","Pycard"),
  lipid_metabolism = c("Srebf1","Srebf2","Fasn","Acaca","Acacb","Scd1","Cpt1a","Cpt2","Ppara","Pparg","Lpl","Cd36","Hadha","Hadhb","Abca1","Abcg1","Dgat1","Dgat2","Pnpla2","Lipe"),
  oxidative_stress_ROS = c("Txnip","Nfe2l2","Keap1","Sod1","Sod2","Cat","Gpx1","Gpx4","Prdx1","Prdx3","Txn","Txnrd1","Nox1","Cybb","Duoxa1","Duox2","Hmox1","Gclc","Gclm")
)
manual_sets <- map(manual_sets, ~ unique(.x[!is.na(.x) & .x != ""]))

# ------------------------- 5. 兼容各版本 msigdbr ------------------------------
# 不同版本参数名不同：先尝试 collection/subcollection，再回退到 category/subcategory。
get_msigdb <- function(collection, subcollection = NULL, species_name) {
  a <- list(db_species = "HS", species = species_name, collection = collection)
  if (!is.null(subcollection)) a$subcollection <- subcollection
  z <- tryCatch(do.call(msigdbr, a), error = function(e) NULL)
  if (!is.null(z) && nrow(z) > 0) return(z)
  a <- list(species = species_name, category = collection)
  if (!is.null(subcollection)) a$subcategory <- subcollection
  z <- tryCatch(do.call(msigdbr, a), error = function(e) NULL)
  if (is.null(z)) tibble() else z
}
standardise_msigdb_columns <- function(x) {
  if (nrow(x) == 0) return(tibble())
  g <- intersect(c("gene_symbol","human_gene_symbol","mouse_gene_symbol"), names(x))[1]
  s <- intersect(c("gs_name","gene_set"), names(x))[1]
  if (is.na(g) || is.na(s)) stop("msigdbr 返回结果缺少基因符号或基因集名称列")
  x %>% transmute(gs_name = .data[[s]], gene_symbol = .data[[g]])
}
message("读取 MSigDB（需要联网）……")
msig_h <- standardise_msigdb_columns(get_msigdb("H", species_name = msig_species))
msig_reactome <- standardise_msigdb_columns(get_msigdb("C2", "CP:REACTOME", msig_species))
msig_go <- standardise_msigdb_columns(get_msigdb("C5", "GO:BP", msig_species))

# ------------------------- 6. 关键词提取 MSigDB 基因 ----------------------------
collect_keyword_genes <- function(tbl, pattern) {
  if (nrow(tbl) == 0) return(character())
  tbl %>% filter(str_detect(gs_name, regex(pattern, ignore_case = TRUE))) %>% pull(gene_symbol) %>% unique()
}
# 关键词覆盖 Hallmark、Reactome 和 GO；GO 名称可能带 GO_ 或 GOBP_ 前缀。
msig_pyroptosis <- c(collect_keyword_genes(msig_h, "INFLAMMATORY_RESPONSE|TNFA_SIGNALING_VIA_NFKB"), collect_keyword_genes(msig_reactome, "PYROPTOSIS|INFLAMMASOME|NLRP3|CASPASE|INTERLEUKIN_1|INTERLEUKIN_18"), collect_keyword_genes(msig_go, "PYROPTOSIS|INFLAMMASOME|REACTIVE_OXYGEN|REGULATION_OF_INFLAMMATORY_RESPONSE"))
msig_inflammasome <- c(collect_keyword_genes(msig_reactome, "INFLAMMASOME|NLRP3|AIM2|NLRC4|CASPASE-1|INTERLEUKIN-1|INTERLEUKIN-18"), collect_keyword_genes(msig_go, "INFLAMMASOME|NLRP3|CASPASE-1|INTERLEUKIN-1|INTERLEUKIN-18"))
msig_lipid <- c(collect_keyword_genes(msig_h, "FATTY_ACID_METABOLISM|CHOLESTEROL_HOMEOSTASIS|ADIPOGENESIS"), collect_keyword_genes(msig_reactome, "FATTY_ACID|LIPID|CHOLESTEROL|TRIGLYCERIDE|BETA_OXIDATION"), collect_keyword_genes(msig_go, "LIPID_METABOLIC|FATTY_ACID|CHOLESTEROL|TRIGLYCERIDE|LIPOPROTEIN"))
msig_ros <- c(collect_keyword_genes(msig_h, "REACTIVE_OXYGEN_SPECIES|HYPOXIA|P53_PATHWAY"), collect_keyword_genes(msig_reactome, "REACTIVE_OXYGEN|OXIDATIVE_STRESS|ANTIOXIDANT|NFE2L2|GLUTATHIONE"), collect_keyword_genes(msig_go, "REACTIVE_OXYGEN|OXIDATIVE_STRESS|CELLULAR_RESPONSE_TO_OXIDATIVE|ANTIOXIDANT|GLUTATHIONE"))

# 合并手工锚定基因和数据库基因，并去重排序。
gene_sets <- list(
  pyroptosis = unique(c(manual_sets$pyroptosis, msig_pyroptosis)),
  inflammasome = unique(c(manual_sets$inflammasome, msig_inflammasome)),
  lipid_metabolism = unique(c(manual_sets$lipid_metabolism, msig_lipid)),
  oxidative_stress_ROS = unique(c(manual_sets$oxidative_stress_ROS, msig_ros))
)
gene_sets <- map(gene_sets, ~ sort(unique(.x[!is.na(.x) & .x != ""])))

# ------------------------- 7. 覆盖率和核心轴质控 -------------------------------
check_gene_set_coverage <- function(gene_set, input_genes = core_axis) {
  input_genes <- unique(input_genes); present <- intersect(input_genes, gene_set)
  tibble(n_input = length(input_genes), n_present = length(present), coverage = n_present / length(input_genes), missing_core_axis = paste(setdiff(input_genes, present), collapse = ";"))
}
coverage_table <- imap_dfr(gene_sets, ~ check_gene_set_coverage(.x) %>% mutate(gene_set = .y, .before = 1))
print(coverage_table)
# 主轴必须完整保留在两个核心集合中；失败时停止，防止静默错误。
stopifnot(all(core_axis %in% gene_sets$pyroptosis), all(core_axis %in% gene_sets$inflammasome))

# ------------------------- 8. 审计表和来源记录 ---------------------------------
gene_sets_long <- imap_dfr(gene_sets, ~ tibble(gene_set = .y, gene_symbol = .x))
gene_set_summary <- gene_sets_long %>% count(gene_set, name = "n_genes")
pathway_membership_audit <- gene_sets_long %>% mutate(is_core_axis = gene_symbol %in% core_axis, is_manual_anchor = map2_lgl(gene_set, gene_symbol, ~ .y %in% manual_sets[[.x]]))
gene_set_provenance <- tibble(organism = organism, msig_species = msig_species, msigdb_version = tryCatch(as.character(msigdbr::msigdbr_version()), error = function(e) NA_character_), source_layers = "manual anchors; MSigDB H; C2 Reactome; C5 GO:BP", core_axis = paste(core_axis, collapse = "–"), created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# ------------------------- 9. 导出文件 ----------------------------------------
# RDS 供 R 分析；GMT 供 GSEA；CSV 便于人工审阅和跨软件使用。
saveRDS(list(organism = organism, core_axis = core_axis, gene_sets = gene_sets, coverage = coverage_table, provenance = gene_set_provenance), file.path(out_dir, "gene_sets.rds"))
gmt_lines <- imap_chr(gene_sets, ~ paste(c(.y, "MANUAL_AND_MSIGDB", .x), collapse = "\t"))
writeLines(enc2utf8(gmt_lines), file.path(out_dir, "gene_sets.gmt"), useBytes = TRUE)
write_csv(gene_sets_long, file.path(out_dir, "gene_sets_long.csv"))
write_csv(gene_set_summary, file.path(out_dir, "gene_set_summary.csv"))
write_csv(pathway_membership_audit, file.path(out_dir, "pathway_membership_audit.csv"))
write_csv(gene_set_provenance, file.path(out_dir, "gene_set_provenance.csv"))
write_csv(coverage_table, file.path(out_dir, "gene_set_coverage.csv"))
