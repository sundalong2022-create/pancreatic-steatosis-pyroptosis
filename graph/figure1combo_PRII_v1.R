# ============================================================
#  Figure 1  (1a + 1b + 1c)   —— 170 mm 宽，合计 110 mm 高
#  直接整段运行，输出 Fig1a/1b/1c 单图 + 拼装版 Fig1
# ============================================================
library(ggplot2)

## ---------- Part 1：Figure 1a  五分组合饲料设计（170 x 34 mm） ----------
dat <- data.frame(
  group        = c("CD", "Safflower", "HO-Safflower", "Sunflower", "HO-Sunflower"),
  x            = 1:5,
  fat_kcal_pct = c(18, 45, 45, 45, 45),
  oil_source   = c("Control", "Safflower", "Safflower", "Sunflower", "Sunflower"),
  oleic_level  = c("", "Conventional", "High-oleic", "Conventional", "High-oleic"),
  n            = c(5, 5, 4, 4, 3),
  role         = c("Fat-quantity reference\n(descriptive only)",
                   "Primary contrast", "Primary contrast",
                   "Parallel comparison", "Parallel comparison"),
  fill         = c("#BDBDBD", "#D6604D", "#B2182B", "#92C5DE", "#2166AC"),
  stringsAsFactors = FALSE
)
dat$lab_col <- ifelse(dat$group == "CD", "grey20", "white")

brk <- function(x1, x2, y, drop = 0.09) {
  data.frame(x    = c(x1, x1, x2),
             xend = c(x1, x2, x2),
             y    = c(y + drop, y, y),
             yend = c(y, y, y + drop))
}

p1a <- ggplot(dat) +
  geom_tile(aes(x = x, y = 1, fill = fill),
            width = 0.82, height = 0.62, colour = "grey25", linewidth = 0.3) +
  scale_fill_identity() +
  geom_text(aes(x = x, y = 1.13, label = group),
            size = 2.1, fontface = "bold", colour = "grey10", lineheight = 0.9) +
  geom_text(aes(x = x, y = 0.90, label = paste0("n = ", n), colour = lab_col),
            size = 2.0) +
  scale_colour_identity() +
  geom_text(aes(x = x, y = 0.60, label = role),
            size = 1.8, colour = "grey35", vjust = 1, lineheight = 0.9) +
  ## 高脂四臂标注
  annotate("segment", x = 2, xend = 5, y = 1.60, yend = 1.60,
           linewidth = 0.4, colour = "grey25") +
  annotate("text", x = 3.5, y = 1.68, label = "High-fat arms: 45% kcal from fat",
           size = 2.1, colour = "grey25", vjust = 0) +
  ## CD 标注
  annotate("segment", x = 0.85, xend = 1.15, y = 1.60, yend = 1.60,
           linewidth = 0.4, colour = "grey55") +
  annotate("text", x = 1, y = 1.68, label = "CD: 18% kcal from fat",
           size = 2.1, colour = "grey45", vjust = 0) +
  ## 分隔虚线：脂肪量轴 ≠ 脂肪酸组成轴
  annotate("segment", x = 1.5, xend = 1.5, y = 0.45, yend = 1.78,
           linetype = "dashed", linewidth = 0.3, colour = "grey60") +
  ## 主对比 bracket
  geom_segment(data = brk(1.62, 3.38, 0.34),
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.4, colour = "#B2182B") +
  annotate("text", x = 2.50, y = 0.22,
           label = "Primary contrast\nisoenergetic, 45% kcal;\nfatty acid composition only",
           size = 1.8, colour = "#B2182B", fontface = "bold",
           vjust = 1, lineheight = 0.9) +
  ## 平行比较 bracket
  geom_segment(data = brk(3.62, 5.38, 0.34),
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.4, colour = "#2166AC") +
  annotate("text", x = 4.50, y = 0.22,
           label = "Parallel comparison\n(second oil source)",
           size = 1.8, colour = "#2166AC", vjust = 1, lineheight = 0.9) +
  annotate("text", x = 0.85, y = 0.22, label = "Total n = 21",
           size = 2.0, colour = "grey25", hjust = 0, vjust = 1) +
  coord_cartesian(xlim = c(0.4, 5.6), ylim = c(-0.35, 2.02), clip = "off") +
  theme_void(base_size = 7) +
  theme(plot.margin = margin(2, 4, 2, 4),
        text = element_text(family = "Arial"))
p1a

## ---------- Part 2：Figure 1b  喂养时间轴（170 x 22 mm） ----------
x0 <- 6; x1 <- 18

## 遗传背景：所有字面量一律 plain()，避免 plotmath 解析报错
geno_expr <- 'Male~plain("C57BL/6J")~italic(Kras)^plain("LSL-G12D/+")*plain(";")~italic(Pdx1)*plain("-")*plain("Cre")'
## 若上标仍异常，可整行替换为：
## geno_expr <- "Male C57BL/6J KrasLSL-G12D/+; Pdx1-Cre"
## 并把下面那一层的 parse = TRUE 改为 parse = FALSE

p1b <- ggplot() +
  annotate("segment", x = x0, xend = x1, y = 1, yend = 1,
           linewidth = 8, colour = "#EFEFEF", lineend = "butt") +
  annotate("segment", x = x0, xend = x1 - 0.18, y = 1, yend = 1,
           linewidth = 0.65, colour = "grey20",
           arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  annotate("point", x = x0, y = 1, size = 1.6, colour = "grey10") +
  annotate("text", x = 6, y = 2.45, label = geno_expr, parse = TRUE,
           hjust = 0, size = 2.0, colour = "grey15") +
  annotate("text", x = 12, y = 1.78,
           label = "12 weeks on assigned diet (n = 21)",
           size = 2.3, colour = "grey25") +
  annotate("segment", x = c(6, 9, 12, 15, 18), xend = c(6, 9, 12, 15, 18),
           y = 1.00, yend = 1.06, linewidth = 0.25, colour = "grey55") +
  annotate("text", x = c(6, 9, 12, 15, 18), y = 1.22,
           label = c("6", "9", "12", "15", "18"),
           size = 1.9, colour = "grey45") +
  annotate("text", x = 18, y = 1.78, label = "Age (weeks)",
           hjust = 1, size = 1.9, colour = "grey45") +
  annotate("segment", x = x0, xend = x0, y = 0.82, yend = 0.52,
           linewidth = 0.3, colour = "grey55") +
  annotate("text", x = 6, y = 0.42,
           label = "6 wk of age: assigned to one of five diets",
           hjust = 0, vjust = 1, size = 2.1, colour = "grey20") +
  annotate("segment", x = x1, xend = x1, y = 0.82, yend = 0.52,
           linewidth = 0.3, colour = "grey55") +
  annotate("text", x = 18, y = 0.42,
           label = "18 wk of age: pancreatic tissue collected",
           hjust = 1, vjust = 1, size = 2.1, colour = "grey20") +
  coord_cartesian(xlim = c(5.2, 18.8), ylim = c(0.05, 2.90), clip = "off") +
  theme_void(base_size = 7) +
  theme(plot.margin = margin(2, 4, 2, 4),
        text = element_text(family = "Arial"))
p1b

## ---------- Part 3：Figure 1c  分析流程图（170 x 54 mm） ----------
L0 <- 8;  L1x <- 74
R0 <- 96; R1x <- 162
cx_l <- (L0 + L1x) / 2      # 41
cx_r <- (R0 + R1x) / 2      # 129

rows <- c(48.5, 41.5, 34.5, 27.5, 20.5, 13.5)   # 行距 7
bh   <- 6                                       # 盒高（原 5 → 加高到 6）

mk <- function(xmin, xmax, yc, label, fill = "white", border = "grey45") {
  data.frame(xmin = xmin, xmax = xmax,
             ymin = yc - bh/2, ymax = yc + bh/2,
             xc = (xmin + xmax)/2, yc = yc,
             label = label, fill = fill, border = border,
             stringsAsFactors = FALSE)
}

nodes <- rbind(
  ## 左列：小鼠 bulk
  mk(L0, L1x, rows[1], "Bulk pancreatic transcriptomes\n(normalisation set, n = 21)", "#F6F6F6"),
  mk(L0, L1x, rows[2], "Row-wise z-scores\ncomputed once and held fixed", "white"),
  mk(L0, L1x, rows[3], "PRII = mean of six gene-wise z-scores", "#F6F6F6"),
  mk(L0, L1x, rows[4], "Primary contrast (n = 9):\nHO-Safflower vs Safflower", "#FBE9E7", "#B2182B"),
  mk(L0, L1x, rows[5], "Bootstrap / permutation /\nleave-one-out / ssGSEA", "white"),
  mk(L0, L1x, rows[6], "Immune module associations", "#F6F6F6"),
  ## 右列：人类公共数据（两条互不连通的子链）
  mk(R0, R1x, rows[1], "TCGA-PAAD; GSE79668\n(human bulk tumours)", "#F6F6F6"),
  mk(R0, R1x, rows[2], "Cross-cohort replication of PRII", "white"),
  mk(R0, R1x, rows[3], "Deconvolution:\nxCell / EPIC / ESTIMATE", "#F6F6F6"),
  mk(R0, R1x, rows[4], "Purity-adjusted overall-survival\nanalysis", "#F6F6F6"),
  mk(R0, R1x, rows[5], "GSE154778: 14,926 cells\n(human scRNA-seq)", "white"),
  mk(R0, R1x, rows[6], "UCell / AUCell / AddModuleScore\nmyeloid localisation; pseudobulk", "#F6F6F6")
)

ed <- function(x, ya, yb) data.frame(x = x, y = ya, yend = yb)
arrows_df <- rbind(
  ed(cx_l, rows[1] - bh/2, rows[2] + bh/2),
  ed(cx_l, rows[2] - bh/2, rows[3] + bh/2),
  ed(cx_l, rows[3] - bh/2, rows[4] + bh/2),
  ed(cx_l, rows[4] - bh/2, rows[5] + bh/2),
  ed(cx_l, rows[5] - bh/2, rows[6] + bh/2),
  ed(cx_r, rows[1] - bh/2, rows[2] + bh/2),
  ed(cx_r, rows[2] - bh/2, rows[3] + bh/2),
  ed(cx_r, rows[3] - bh/2, rows[4] + bh/2),
  ed(cx_r, rows[5] - bh/2, rows[6] + bh/2)
)

y21_top <- rows[1] + bh/2
y21_bot <- rows[4] - bh/2
y9_top  <- rows[4] + bh/2
y9_bot  <- rows[4] - bh/2

p1c <- ggplot() +
  geom_rect(data = nodes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill, colour = border), linewidth = 0.3) +
  scale_fill_identity() + scale_colour_identity() +
  geom_text(data = nodes, aes(x = xc, y = yc, label = label),
            size = 2.3, colour = "grey15", lineheight = 0.92) +
  geom_segment(data = arrows_df, aes(x = x, xend = x, y = y, yend = yend),
               linewidth = 0.35, colour = "grey40",
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  ## 列标题
  geom_text(aes(x = cx_l, y = 53.6, label = "Mouse bulk transcriptomes"),
            size = 2.5, fontface = "bold", colour = "#B2182B") +
  geom_text(aes(x = cx_r, y = 53.6, label = "Human public datasets"),
            size = 2.5, fontface = "bold", colour = "#2166AC") +
  annotate("segment", x = L0, xend = L1x, y = 52.3, yend = 52.3, linewidth = 0.4, colour = "#B2182B") +
  annotate("segment", x = R0, xend = R1x, y = 52.3, yend = 52.3, linewidth = 0.4, colour = "#2166AC") +
  ## n = 21 brace
  annotate("segment", x = 77, xend = 77, y = y21_bot, yend = y21_top,
           linewidth = 0.3, colour = "grey55") +
  annotate("segment", x = 77, xend = 78.6, y = y21_top, yend = y21_top, linewidth = 0.3, colour = "grey55") +
  annotate("segment", x = 77, xend = 78.6, y = y21_bot, yend = y21_bot, linewidth = 0.3, colour = "grey55") +
  annotate("text", x = 80.6, y = (y21_top + y21_bot)/2,
           label = "normalisation set (n = 21)", angle = 90, size = 1.9, colour = "grey40") +
  ## n = 9 brace
  annotate("segment", x = 85.5, xend = 85.5, y = y9_bot, yend = y9_top,
           linewidth = 0.3, colour = "#B2182B") +
  annotate("segment", x = 85.5, xend = 87.1, y = y9_top, yend = y9_top, linewidth = 0.3, colour = "#B2182B") +
  annotate("segment", x = 85.5, xend = 87.1, y = y9_bot, yend = y9_bot, linewidth = 0.3, colour = "#B2182B") +
  annotate("text", x = 89.1, y = rows[4],
           label = "primary contrast (n = 9)", angle = 90, size = 1.9, colour = "#B2182B") +
  ## 单细胞支线末端注
  annotate("text", x = cx_r, y = 9,
           label = "cross-species, cross-dataset; absolute scores not comparable",
           size = 1.7, colour = "grey45", fontface = "italic") +
  ## 底部边界带
  annotate("rect", xmin = 8, xmax = 162, ymin = 0, ymax = 7,
           fill = "#F2F2F2", colour = "grey45", linetype = "dashed", linewidth = 0.35) +
  annotate("text", x = 85, y = 5.0, size = 2.2, colour = "grey20",
           label = "Pre-specified gene-set and scoring sensitivity analyses") +
  annotate("text", x = 85, y = 1.8, size = 2.1, colour = "grey40", fontface = "italic",
           label = "Context-dependent associations reported transparently") +
  coord_cartesian(xlim = c(0, 170), ylim = c(0, 55), clip = "off") +
  theme_void(base_size = 7) +
  theme(plot.margin = margin(2, 4, 2, 4),
        text = element_text(family = "Arial"))
p1c

## ---------- Part 4：导出单图 ----------
save_tiff <- function(plot_obj, filename, width_mm, height_mm) {
  tiff(filename, width = width_mm, height = height_mm, units = "mm",
       res = 600, compression = "lzw", type = "cairo")
  print(plot_obj)
  dev.off()
}
ggsave("Fig1_PRII_v1_a.pdf",  p1a, width = 170, height = 34, units = "mm", device = cairo_pdf)
save_tiff(p1a, "Fig1_PRII_v1_a.tiff", 170, 34)
ggsave("Fig1_PRII_v1_b.pdf",  p1b, width = 170, height = 22, units = "mm", device = cairo_pdf)
save_tiff(p1b, "Fig1_PRII_v1_b.tiff", 170, 22)
ggsave("Fig1_PRII_v1_c.pdf",  p1c, width = 170, height = 56, units = "mm", device = cairo_pdf)
save_tiff(p1c, "Fig1_PRII_v1_c.tiff", 170, 56)

## ---------- Part 5：拼装整张 Figure 1 ----------
library(patchwork)

p1a_t <- p1a + labs(tag = "a") + theme(plot.margin = margin(1, 1, 0, 4, unit = "mm"))
p1b_t <- p1b + labs(tag = "b") + theme(plot.margin = margin(1, 1, 0, 4, unit = "mm"))
p1c_t <- p1c + labs(tag = "c") + theme(plot.margin = margin(1, 1, 0, 4, unit = "mm"))

fig1 <- (p1a_t / p1b_t / p1c_t) +
  plot_layout(heights = c(34, 22, 54)) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"),
        plot.tag.position = "topleft")

# 总高 = 110 + 上下新增间距 ≈ 128 mm
ggsave("Fig1_PRII_v1.pdf",  fig1, width = 170, height = 128, units = "mm", device = cairo_pdf)
save_tiff(fig1, "Fig1_PRII_v1.tiff", 170, 128)
## library(patchwork)
## fig1 <- (p1a / p1b / p1c) + plot_layout(heights = c(34, 22, 54))
## ggsave("Fig1_PRII_v1.pdf",  fig1, width = 170, height = 110, units = "mm", device = cairo_pdf)
## ggsave("Fig1_PRII_v1.tiff", fig1, width = 170, height = 110, units = "mm", dpi = 600, compression = "lzw")



## ---------- 5 秒自检 ----------
range(nodes$ymin, nodes$ymax)              # 期望 10.5 – 51.5
invisible(parse(text = geno_expr))         # 不报错即安全
