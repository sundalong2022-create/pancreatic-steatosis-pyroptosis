# 基因集构建说明

本项目固定核心轴：

```text
TXNIP → NLRP3 → CASP1 → GSDMD → IL1B / IL18
```

脚本位置：`scripts/01_build_gene_sets.R`

运行方式：

```r
source("scripts/01_build_gene_sets.R")
```

默认按人类基因符号运行。若表达矩阵为小鼠官方 symbol，可在脚本开头将：

```r
organism <- "human"
```

改为：

```r
organism <- "mouse"
```

输出文件在 `results/gene_sets/`：

- `gene_sets.rds`：R list，供 GSVA/module score 直接使用
- `gene_sets.gmt`：GSEA 可用格式
- `gene_sets_long.csv`：长表，包含 gene set、gene symbol、是否核心轴
- `gene_set_summary.csv`：每个集合的基因数与核心轴数量
- `pathway_membership_audit.csv`：MSigDB 通路来源审计表
- `gene_set_provenance.csv`：物种、来源、msigdbr 版本、生成时间

四类主集合：

1. `pyroptosis`：焦亡通路，包含核心轴与 CASP4/CASP5/GSDME 等补充节点。
2. `inflammasome`：炎症小体，包含 NLRP3、PYCARD、CASP1、IL1B/IL18，以及 NLRP1/NLRC4/AIM2/MEFV 等平行传感器。
3. `lipid_metabolism`：脂质代谢，覆盖脂肪酸合成、脂肪酸氧化、脂滴、胆固醇转运等节点。
4. `oxidative_stress_ROS`：氧化应激/ROS，覆盖 NRF2 防御、抗氧化清除、NOX/DUOX 产 ROS 模块。

注意事项：

- 人类 symbol 通常全大写，如 `TXNIP/NLRP3/CASP1`。
- 小鼠 symbol 使用官方大小写，如 `Txnip/Nlrp3/Casp1`。
- 不要把人类和小鼠 symbol 混在同一次打分中。
- `CASP11` 主要用于小鼠；人类非经典炎症性 caspase 主要是 `CASP4/CASP5`。
- `TXNIP` 在本项目中作为 ROS 到 NLRP3 的耦联节点，不是所有炎症小体的通用结构组件。

与表达矩阵交集时建议报告覆盖率，不要静默删除核心轴：

```r
gene_sets <- readRDS("results/gene_sets/gene_sets.rds")
coverage_tbl <- check_gene_set_coverage(gene_sets, rownames(expr))
coverage_tbl
```
