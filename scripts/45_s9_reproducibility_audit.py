#!/usr/bin/env python3
"""Build the S9 reproducibility and parameter audit from verified project records.

This is a metadata-only compilation. It does not rerun or modify any analysis result.
Run from any working directory: python3 scripts/45_s9_reproducibility_audit.py
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "supplementary" / "S9_reproducibility_audit"
OUT.mkdir(parents=True, exist_ok=True)

PARAM_COLS = [
    "record_id", "module", "dataset", "source_file", "script_or_pipeline",
    "software", "software_version", "input_matrix", "preprocessing",
    "gene_id_rule", "score_method", "normalization_scope",
    "ortholog_mapping_rule", "gene_set_source", "gene_set_version",
    "fdr_family", "missing_value_rule", "sample_count_rule", "random_seed",
    "reported_parameter", "parameter_value", "status", "notes",
]
records: list[dict[str, str]] = []


def add(
    module: str,
    dataset: str,
    parameter: str,
    value: str,
    source: str,
    *,
    status: str = "reported",
    script: str = "",
    software: str = "",
    version: str = "",
    input_matrix: str = "",
    preprocessing: str = "",
    gene_id_rule: str = "",
    score_method: str = "",
    normalization_scope: str = "",
    ortholog: str = "",
    gene_set_source: str = "",
    gene_set_version: str = "",
    fdr: str = "",
    missing: str = "",
    sample_rule: str = "",
    seed: str = "",
    notes: str = "",
) -> None:
    records.append({
        "record_id": f"S9P{len(records)+1:03d}", "module": module,
        "dataset": dataset, "source_file": source, "script_or_pipeline": script,
        "software": software, "software_version": version,
        "input_matrix": input_matrix, "preprocessing": preprocessing,
        "gene_id_rule": gene_id_rule, "score_method": score_method,
        "normalization_scope": normalization_scope,
        "ortholog_mapping_rule": ortholog, "gene_set_source": gene_set_source,
        "gene_set_version": gene_set_version, "fdr_family": fdr,
        "missing_value_rule": missing, "sample_count_rule": sample_rule,
        "random_seed": seed, "reported_parameter": parameter,
        "parameter_value": value, "status": status, "notes": notes,
    })


PARAM_AUDIT = "results/manuscript_parameter_audit.md; results/prii_consistency_audit.md"
PRIMARY = "GSE292815_KC_12wk"
MOUSE_MATRIX = "results/qc/GSE292815_KC_12wk_matrix.rds"

# Study design and provenance
add("Study design", PRIMARY, "animal_sex", "All male", PARAM_AUDIT,
    sample_rule="12-week KC discovery cohort; n=21", notes="Source study/GEO metadata; not inferred from expression.")
add("Study design", PRIMARY, "strain_genotype", "C57BL/6J; KrasLSL-G12D/+; Pdx1-Cre (KC)", PARAM_AUDIT,
    sample_rule="12-week KC discovery cohort; n=21")
add("Study design", PRIMARY, "feeding_timeline", "Feeding began at 6 weeks; 12-week intervention; tissue collected at 18 weeks", PARAM_AUDIT)
add("Study design", PRIMARY, "group_counts", "CD n=5; Safflower n=5; HO-Safflower n=4; Sunflower n=4; HO-Sunflower n=3; total n=21", "data/metadata/sample_metadata.tsv; results/manuscript_parameter_audit.md",
    sample_rule="All 21 annotated 12-week KC samples; no samples excluded")
add("Study design", PRIMARY, "dietary_fat_energy", "CD 18% kcal fat; both HFD arms 45% kcal fat", PARAM_AUDIT,
    notes="The two planned primary comparisons are matched 45%-fat diets differing in fat source; CD is not described as HFD.")
add("Study design", PRIMARY, "allocation_randomization", "not reported", PARAM_AUDIT,
    status="not reported", notes="Available paper description and GEO metadata do not specify randomization, blocking, litter, or cage balancing; do not describe allocation as randomized.")

# Frozen PRII definition
add("PRII definition", PRIMARY, "core_genes", "Txnip; Nlrp3; Casp1; Gsdmd; Il1b; Il18", "scripts/05_psps_score.R; results/12wk_primary/README.md",
    script="scripts/05_psps_score.R", input_matrix=MOUSE_MATRIX, gene_id_rule="Mouse canonical gene symbols; duplicate symbols averaged when reading GEO matrix", score_method="Arithmetic mean of six gene-wise z-scores")
add("PRII definition", PRIMARY, "input_transform", "GEO normalized-count matrix; log2(x+1)", "scripts/02_bulk_de_analysis.R; scripts/05_psps_score.R",
    script="scripts/02_bulk_de_analysis.R; scripts/05_psps_score.R", input_matrix=MOUSE_MATRIX,
    preprocessing="Repeated gene symbols averaged; then log2(normalized expression + 1)", gene_id_rule="Gene_symbol; duplicate rows averaged")
add("PRII definition", PRIMARY, "z_score_formula", "Per gene across the 21 discovery samples, using sample SD (n-1); arithmetic mean across six genes", "scripts/05_psps_score.R; results/manuscript_parameter_audit.md",
    script="scripts/05_psps_score.R", score_method="row-wise z-score followed by column-wise arithmetic mean", normalization_scope="Within each cohort; absolute scores are not compared across cohorts", missing="Non-finite gene z-scores set to 0")
add("PRII definition", PRIMARY, "legacy_name_mapping", "PRII is the unchanged six-gene score historically stored as PSPS_core", "results/prii_rename_audit.md; results/prii_consistency_audit.md",
    status="reported", notes="Semantic rename only; no rescreening, gene replacement, or score recalculation.")
add("PRII definition", "Mouse-human mapping", "ortholog_strategy", "One-to-one canonical symbols: Txnip→TXNIP, Nlrp3→NLRP3, Casp1→CASP1, Gsdmd→GSDMD, Il1b→IL1B, Il18→IL18", "results/composition_lipid_integration/mouse_to_human_ortholog_map.csv; results/manuscript_parameter_audit.md",
    ortholog="One-to-one canonical symbols; species-specific scoring; no cross-species absolute score comparison")

# Discovery-cohort differential expression and enrichment
add("Primary DE", PRIMARY, "expression_filter", "Retain genes with rowSums(log2(expression+1) > 1) >= 2", "scripts/06_12wk_integrated_analysis.R; results/supplementary_S2/S2_data_manifest.txt",
    script="scripts/06_12wk_integrated_analysis.R", input_matrix=MOUSE_MATRIX, preprocessing="log2(normalized expression + 1)", sample_rule="All 21 samples in shared five-group design")
add("Primary DE", PRIMARY, "model_and_contrasts", "Shared five-group limma design; 9 prespecified contrasts", "scripts/06_12wk_integrated_analysis.R",
    script="scripts/06_12wk_integrated_analysis.R", software="R / limma", version="See results/12wk_primary/sessionInfo_12wk_integrated.txt", sample_rule="n=21; group sizes 5/5/4/4/3", seed="292815")
add("Primary DE", PRIMARY, "limma_fit", "eBayes(trend=TRUE, robust=TRUE)", "scripts/06_12wk_integrated_analysis.R",
    script="scripts/06_12wk_integrated_analysis.R", software="limma", version="3.62.2 in recorded primary session", seed="292815")
add("Primary DE", PRIMARY, "BH_family", "Benjamini-Hochberg across all tested genes separately within each contrast", "scripts/06_12wk_integrated_analysis.R; results/supplementary_S2/S2_data_manifest.txt",
    script="scripts/06_12wk_integrated_analysis.R", fdr="9 contrast-specific all-gene families", notes="FDR<0.05 defines statistical significance; displayed DEG lists additionally apply absolute logFC >= 1.")
add("Primary DE", PRIMARY, "display_DEG_threshold", "BH-adjusted P < 0.05 and |logFC| >= 1", "scripts/06_12wk_integrated_analysis.R; results/supplementary_S2/S2_data_manifest.txt")
add("Primary enrichment", PRIMARY, "FGSEA_parameters", "Moderated-t ranking; minSize=3; maxSize=2000; BH within each contrast pathway family", "scripts/06_12wk_integrated_analysis.R; results/supplementary_S2/S2_data_manifest.txt",
    script="scripts/06_12wk_integrated_analysis.R", software="fgsea", version="1.32.4 in recorded primary session", fdr="One pathway-result family per contrast", seed="292815")
add("Primary enrichment", PRIMARY, "ORA_universe_and_gate", "Expressed-gene universe; directional DEG set must contain >=5 genes; BH within contrast × direction × ontology/database", "scripts/06_12wk_integrated_analysis.R; results/supplementary_S2/S2_data_manifest.txt",
    script="scripts/06_12wk_integrated_analysis.R", software="clusterProfiler", version="4.14.6 in recorded primary session", fdr="Separate contrast × direction × ontology/database families", notes="KEGG output depends on available/cached KEGG mapping; see results/12wk_primary/enrichment/KEGG_status.txt.")

# Primary scores and robustness
add("Primary scores", PRIMARY, "score_contrast_family", "Five scores tested per planned contrast; BH within each contrast", "scripts/05_psps_score.R",
    script="scripts/05_psps_score.R", software="R / limma", version="See results/12wk_primary/sessionInfo_psps_12wk.txt", score_method="Five prespecified pathway/core scores", fdr="2 contrast-specific families, 5 scores each", sample_rule="21 complete discovery samples")
add("Primary scores", PRIMARY, "pairwise_score_correlations", "10 pairwise Spearman tests among five scores; BH across all 10", "scripts/05_psps_score.R", fdr="One family of 10 pairwise score correlations")
add("Primary scores", PRIMARY, "core_gene_tests", "Six core genes; BH within each planned contrast", "scripts/05_psps_score.R", fdr="2 contrast-specific families, 6 core genes each")
add("Robustness", PRIMARY, "random_gene_null", "5,000 six-gene draws from expressed non-core genes with nonzero cross-sample SD; environment override PSPS_PERM_N", "scripts/14_psps_robustness.R; results/12wk_primary/robustness/README.md",
    script="scripts/14_psps_robustness.R", software="R", version="See results/12wk_primary/robustness/sessionInfo_robustness.txt", seed="292815", sample_rule="Discovery n=21")
add("Robustness", PRIMARY, "label_permutation", "5,000 label permutations within each matched oil-background comparison; override PSPS_LABEL_PERM_N", "scripts/14_psps_robustness.R; results/12wk_primary/robustness/README.md",
    script="scripts/14_psps_robustness.R", seed="292815", notes="Empirical permutation P values are reported for the prespecified contrast; permutation does not establish randomized allocation.")
add("Robustness", PRIMARY, "bootstrap", "5,000 replicates by default; within-group resampling for score contrasts and sample-level resampling for correlations; override PSPS_BOOTSTRAP_B", "scripts/14_psps_robustness.R; results/12wk_primary/robustness/README.md",
    script="scripts/14_psps_robustness.R", seed="292815", notes="Bootstrap intervals are not multiplicity-adjusted.")
add("Robustness", PRIMARY, "alternative_ssGSEA", "GSVA ssGSEA; alpha=0.25; normalize=TRUE", "scripts/14_psps_robustness.R; results/12wk_primary/robustness/README.md",
    script="scripts/14_psps_robustness.R", software="GSVA", version="2.0.7 in recorded primary session", score_method="ssGSEA", normalization_scope="Within discovery cohort", seed="292815")
add("Robustness", PRIMARY, "formal_VIF_and_overlap_adjusted_partial_correlation", "not estimated", "results/prii_supplementary_tables_plan_v1_zh_en.md",
    status="not estimated", notes="Do not imply that PRII has formal multicollinearity-adjusted incremental value from these missing outputs.")

# Mouse immune modules and external human PAAD replication
add("Mouse immune association", PRIMARY, "module_score", "Mean gene-wise z-score for each marker module after excluding the six PRII genes", "scripts/09_12wk_immune_microenvironment.R; results/12wk_primary/immune/README.md",
    script="scripts/09_12wk_immune_microenvironment.R", input_matrix=MOUSE_MATRIX, score_method="Mean of present marker-gene z-scores", sample_rule="All 21 discovery samples")
add("Mouse immune association", PRIMARY, "Spearman_BH_family", "Nine prespecified non-overlapping modules in one BH family", "scripts/09_12wk_immune_microenvironment.R",
    script="scripts/09_12wk_immune_microenvironment.R", fdr="9 module-correlation tests")
add("Mouse immune association", PRIMARY, "group_contrast_BH", "Two planned contrasts; limma BH within the nine module outcomes for each contrast", "scripts/09_12wk_immune_microenvironment.R",
    script="scripts/09_12wk_immune_microenvironment.R", fdr="2 contrast-specific families, 9 modules each", seed="292815")
add("External human PAAD", "GSE79668", "expression_preprocessing", "Duplicate gene symbols summed; log2(CPM+1); cohort-wise gene z-score mean", "scripts/external_validation_gse79668.py; results/external_validation/README.md",
    script="scripts/external_validation_gse79668.py", software="Python; pandas/numpy/scipy/statsmodels", version="not recorded", input_matrix="data/external/GSE79668_counts.txt.gz", preprocessing="Raw counts to CPM and log2(CPM+1)", gene_id_rule="Uppercase human symbols; duplicate symbols summed", score_method="PRII mean-z; immune module mean-z", normalization_scope="Within GSE79668 cohort", fdr="BH across 8 immune modules", sample_rule="n=51; no usable OS time/event pair", status="partially reported")

# TCGA and deconvolution
add("TCGA-PAAD", "TCGA-PAAD", "expression_preprocessing", "Primary Tumor; exact sample barcode; tpm_unstrand; uppercase symbols; first duplicate retained; log2(TPM+1)", "scripts/08_tcga_paad_external_validation.R; results/tcga_paad_external/analysis_manifest.txt",
    script="scripts/08_tcga_paad_external_validation.R", input_matrix="TCGA_PAAD_RNAseq.rds", preprocessing="log2(TPM+1)", gene_id_rule="Uppercase symbols; first duplicate retained", score_method="Cohort-wise gene z-score mean", normalization_scope="Within TCGA-PAAD cohort", fdr="BH across 11 signature scores for signature associations", sample_rule="178 primary tumors")
add("TCGA-PAAD", "TCGA-PAAD", "complete_case_counts", "Primary tumor n=178; immune/purity shared complete n=149; survival n=160 (93 deaths); age+stage n=132; age+stage+purity n=111", "results/immune_validation/tcga_sample_count_audit.csv; results/tcga_paad_external/tcga_sample_count_audit.csv; results/prii_consistency_audit.md",
    sample_rule="Each count is an analysis-specific complete-case cohort; not a sequential attrition series")
add("TCGA immune", "TCGA-PAAD", "purity_adjusted_partial_correlation", "Rank-residual Spearman partial correlation; two-sided t approximation with df=n-3 for one covariate", "scripts/09_immune_microenvironment_validation.R; results/immune_validation/README.md",
    script="scripts/09_immune_microenvironment_validation.R", sample_rule="Exact-barcode matched purity; shared complete n=149", fdr="BH inherited separately by analysis family from source metric tables")
add("TCGA survival", "TCGA-PAAD", "survival_models", "Continuous PRII Cox; survival complete n=160, 93 events; age+stage n=132; age+stage+purity n=111; median-split KM descriptive", "results/tcga_paad_external/analysis_manifest.txt; results/prii_consistency_audit.md",
    script="scripts/08_tcga_paad_external_validation.R", fdr="No Cox multiplicity correction reported", sample_rule="Model-specific complete-case counts; KM is descriptive")
add("TCGA deconvolution", "TCGA-PAAD", "methods_and_convergence", "xCell, EPIC, ESTIMATE, local LM22-compatible SVR (permutations=0); EPIC converged in 167/178; official CIBERSORTx not run", "results/supplementary/S5_tcga_sensitivity/S5_run_manifest.txt; results/tcga_paad_external/deconvolution/CIBERSORTx_input_manifest.txt",
    script="scripts/16_tcga_formal_deconvolution.R; scripts/17_tcga_deconvolution_report.R; scripts/18_epic_convergence_sensitivity.R", sample_rule="Metric-specific complete cases; EPIC sensitivity excludes 11 non-converged samples", fdr="BH inherited from each full metric family; method comparisons BH across 35 pairs")

# Single-cell S4 and scoring-method boundaries
add("Single-cell context", "GSE154778", "dataset_and_statistical_unit", "14,926 cells; 16 samples (10 primary, 6 metastatic); inference at sample × cell type", "results/supplementary/S4_single_cell_context/S4_run_manifest.txt",
    input_matrix="data/raw/GSE154778/GSE154778_dgeMtx.csv.gz", sample_rule="16 samples; do not treat individual cells as independent biological replicates")
add("Single-cell context", "GSE154778", "UCell_and_pseudobulk", "UCell on complete gene-ranking background; maxRank=1500; ties averaged; sample_id × cell_type pseudobulk; minimum 10 cells per unit", "results/supplementary/S4_single_cell_context/S4_run_manifest.txt; scripts/19_single_cell_psps_upgrade.R",
    script="scripts/19_single_cell_psps_upgrade.R", software="UCell", version="Exact S4 runtime version not captured; UCell 2.10.1 appears in a separate earlier upgrade session", input_matrix="data/raw/GSE154778/GSE154778_dgeMtx.csv.gz", score_method="UCell", sample_rule="Pseudobulk units require >=10 cells; mixed model includes sample random intercept", fdr="BH within each six-comparison macrophage-reference model family", status="partially reported")
add("Single-cell context", "GSE154778", "AUCell_AddModuleScore_scope", "No S4 AUCell, AddModuleScore, or composition-adjusted PRII–immune association was performed", "results/supplementary/README_S3_zh_en.md; results/supplementary/S4_single_cell_context/S4_run_manifest.txt",
    status="not estimated", notes="Do not report these as S4 evidence.")
add("Scoring sensitivity", "GSE292815 and GSE154778", "AddModuleScore_compatible_parameters", "24 expression bins; up to 100 same-bin controls per target; candidate shortage sampled with replacement; empty bin falls back to all non-target genes; controls deduplicated; target minus control mean", "scripts/22_bulk_sc_composition_lipid_integration.R; results/manuscript_parameter_audit.md",
    script="scripts/22_bulk_sc_composition_lipid_integration.R", score_method="Seurat AddModuleScore-compatible implementation", seed="292815", sample_rule="Per cohort; bulk uses log2(expression+1), single-cell uses log1p(CPM×10^4)")
add("Scoring sensitivity", "GSE292815", "AUCell_bulk_cutoff", "AUCell cutoff grid 5%, 10%, 20%, 50%, 80%; conventional reference 5%; aucMaxRank=ceil(fraction × 17,830 genes)", "scripts/44_prii_aucell_bulk_cutoff_audit.R; results/supplementary/S3_prii_robustness/AUCell_run_manifest.txt",
    script="scripts/44_prii_aucell_bulk_cutoff_audit.R", software="AUCell", version="1.34.0", input_matrix=MOUSE_MATRIX, score_method="AUCell rank-based AUC", seed="20261011", sample_rule="17,830 genes × 21 samples", notes="At 5%, incomplete six-gene coverage prevents interpretation as full-signature agreement.")

# Independent mouse metabolic context cohorts
add("External mouse context", "GSE234956", "cohort_and_model", "45 WT mouse islet samples; 44% fat+fiber versus 10% fat+fiber; model adjusted for strain and sex; exposure duration unavailable", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; scripts/external_validation_gse234956.py",
    script="scripts/external_validation_gse234956.py", score_method="Within-cohort six-gene mean-z", fdr="Adjusted-table BH across 14 features; per-strain Mann–Whitney nominal", sample_rule="n=45; public metadata does not report duration", status="partially reported")
add("External mouse context", "GSE190812", "cohort_and_analysis", "Three arms; n=7 pooled-islet libraries per arm; each library combines islets from five mice; nominal Mann–Whitney; no meta-analysis", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; results/external_validation/",
    script="No tracked score-generation script located; S6 figure script reads precomputed outputs", score_method="Within-cohort six-gene mean-z; exact source preprocessing not recoverable", fdr="Nominal P values; no FDR reported", sample_rule="Library is the analysis unit; each library pools five mice", status="partially reported")
add("External mouse context", "GSE249790", "cohort_and_analysis", "Sorted beta cells; CD n=4 versus 2-day HFD n=3; Mann–Whitney; BH across 7 rows (PRII + six genes)", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; results/external_validation/",
    script="No tracked score-generation script located; S6 figure script reads precomputed outputs", score_method="Within-cohort six-gene mean-z; exact source preprocessing not recoverable", fdr="One BH family of 7 rows", sample_rule="n=4 vs n=3", status="partially reported", notes="Several genes at zero/low-count floor; negative direction is difficult to interpret.")
add("External mouse context", "GSE293790", "cohort_and_analysis", "Female KC/Kras whole-pancreas model; CD, DIO, and DIO→CD n=3 per arm; pairwise Mann–Whitney nominal; no gene-wise tests", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; scripts/external_validation_gse293790.py",
    script="scripts/external_validation_gse293790.py", score_method="Within-cohort six-gene mean-z", fdr="Nominal pairwise tests; no BH", sample_rule="n=3 per group", notes="Cancer-susceptible model, not normal fatty pancreas.")
add("External mouse context", "GSE234956; GSE190812; GSE249790; GSE293790", "cross_cohort_pooling", "No pooled effect estimate or cross-cohort meta-analysis", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt")

# Exploratory lipid programme and network audit
add("Exploratory lipid programme", PRIMARY, "gene_set_and_score", "Four external programmes: Hallmark fatty-acid metabolism, Hallmark adipogenesis, Reactome fatty-acid metabolism, KEGG PPAR signalling; GSVA ssGSEA alpha=0.25, normalize=TRUE", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt; results/12wk_primary/lipid_external/README.md",
    script="scripts/20_external_lipid_ssgsea.R", software="GSVA", version="S7 manifest says 2.x; exact S7 runtime version not captured", score_method="ssGSEA", normalization_scope="Within 21-sample GSE292815 cohort", ortholog="MSigDB human definitions mapped to mouse symbols via msigdbr", gene_set_source="MSigDB", gene_set_version="2026.1.Hs for the four audited external lipid sets only", sample_rule="n=21", status="partially reported")
add("Exploratory lipid programme", PRIMARY, "global_and_sensitivity_FDR", "Main global BH family=12 tests (4 pooled correlations + 8 diet contrasts); subscore pooled=16; group-adjusted=16; within-group diagnostic=80", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt; results/12wk_primary/lipid_external/README.md",
    fdr="Distinct prespecified families: 12, 16, 16, and 80 tests", sample_rule="Within-group n=3–5; diagnostic only, not plotted", notes="Do not call every S7 correction a single global FDR.")
add("Network audit", "GSE292815 / TCGA-PAAD / GSE154778", "evidence_layers_and_integrated_counts", "189 nodes; 253 integrated edges = 96 expression association + 144 database prior + 13 candidate ligand–receptor expression support", "results/supplementary/S8_network_evidence_audit/S8_run_manifest.txt; results/supplementary/S8_network_evidence_audit/integrated_edges_normalized.csv",
    sample_rule="TF bulk n=21; TCGA sensitivity n=149; ligand–receptor GSE154778 n=16", notes="Evidence labels are not causal/functional validation.")
add("Network audit", "Mus musculus STRING", "STRING_threshold", "Required combined score >=0.70; undirected edges", "results/supplementary/S8_network_evidence_audit/S8_run_manifest.txt; scripts/12_string_ppi_network.R",
    script="scripts/12_string_ppi_network.R", fdr="Not applicable to database-prior edge threshold")
add("Network audit", "DoRothEA / bulk and TCGA", "TF_family", "12 prespecified TF candidates; DoRothEA confidence A–C; observed-target expression and TF activity/expression association", "results/network_tf/analysis_manifest.txt; results/supplementary/S8_network_evidence_audit/S8_run_manifest.txt",
    script="scripts/10_tf_target_network.R", fdr="BH across 12 prespecified TF tests", sample_rule="Bulk n=21; TCGA sensitivity n=149")
add("Network audit", "OmniPath ligrecextra / GSE154778", "LR_prior_and_expression_support", "OmniPath prior; sample-level expression support uses mean log1p(CPM), detection fraction, and geometric-mean ligand/receptor expression; database release snapshot not recorded", "results/network_ligand_receptor/run_manifest.txt; results/supplementary/S8_network_evidence_audit/S8_run_manifest.txt",
    script="scripts/11_sc_ligand_receptor_network.py", software="Python packages", version="not recorded", sample_rule="16 samples", status="partially reported", notes="Prior-only and sample-supported edges are kept as separate evidence classes.")

# Explicit project-wide reproducibility gaps / status markers
add("Completeness", "GSE292815", "raw_FASTQ_QC", "not assessed; normalized GEO matrices only", "results/supplementary/S1_qc/S1_run_manifest.txt",
    status="not available", notes="No FASTQ-level QC or raw sequencing depth claim.")
add("Completeness", "GSE292815", "objective_outlier_rule", "No objective outlier-exclusion rule was applied; no samples removed", "results/supplementary/S1_qc/S1_run_manifest.txt")
add("Completeness", "GSE190812 and GSE249790", "score_generation_pipeline", "No tracked source scoring scripts or exact preprocessing pipeline recoverable from current project artifacts", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; scripts/",
    status="not recoverable", notes="S6 figure-generation script consumes precomputed score/test tables only.")
add("Completeness", "S1–S8 modules", "software_version_scope", "Recorded R versions differ by run: R 4.4.0 in primary/upstream analyses and R 4.6.1 in later supplementary audits; report by module", "results/**/sessionInfo*.txt; results/**/run_manifest.txt",
    status="reported with module-specific variation", notes="Do not present one project-wide R version.")
add("Completeness", "S4 GSE154778", "exact_UCell_version", "not captured in S4-specific sessionInfo; UCell 2.10.1 is in a separate earlier upgrade session only", "results/supplementary/S4_single_cell_context/S4_sessionInfo.txt; results/sc_location_GSE154778/upgrade/sessionInfo.txt",
    status="partially reported", notes="The earlier session does not prove the exact S4 runtime package version.")
add("Completeness", "S7 GSE292815", "exact_GSVA_runtime_version", "not captured in S7-specific sessionInfo; manifest says GSVA 2.x", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt; results/supplementary/S7_exploratory_lipid_context/S7_sessionInfo.txt",
    status="partially reported", notes="GSVA 2.0.7 exists in an earlier primary lipid session, but is not asserted as the exact S7 runtime.")
add("Completeness", "Python analyses", "Python_dependency_lock", "Python runtime/package versions and a requirements/lock file are not recorded for cohort scripts", "scripts/external_validation_gse79668.py; scripts/external_validation_gse234956.py; results/network_ligand_receptor/run_manifest.txt",
    status="not reported")
add("Completeness", "STRING / OmniPath / DoRothEA", "database_snapshot_versions", "Exact database release/snapshot versions are not captured in module-specific manifests", "results/network_ppi/README.md; results/network_tf/analysis_manifest.txt; results/network_ligand_receptor/run_manifest.txt",
    status="not reported")
add("Completeness", "TCGA-PAAD", "official_CIBERSORTx", "not run; available LM22-compatible output is local SVR with permutations=0", "results/tcga_paad_external/deconvolution/CIBERSORTx_input_manifest.txt; results/supplementary/S5_tcga_sensitivity/S5_run_manifest.txt")
add("Completeness", "GSE79668", "survival_validation", "not estimated; usable overall-survival time/event pair unavailable", "scripts/external_validation_gse79668.py; results/external_validation/README.md",
    status="not estimated")
add("Completeness", "S5 TCGA survival", "incremental_model_metrics", "delta_LRT, delta_AIC, and C_index not estimated in the current sensitivity package", "results/prii_supplementary_tables_plan_v1_zh_en.md",
    status="not estimated")
add("Completeness", "S10 incremental value", "incremental_signature_analysis", "analysis plan only; no complete nested-model / external replication results", "results/prii_incremental/statistical_analysis_plan_zh_en.md; results/prii_incremental_value_analysis_plan_zh_en.md",
    status="pending", notes="Do not claim independent incremental value until the planned analysis is completed.")

# Capture every explicit set.seed call and keep figure-generation seeds distinct
# from stochastic analysis seeds.
SEED_COLS = ["source_file", "line", "seed_expression", "source_context", "scope_note"]
seed_rows = []
seed_pattern = re.compile(r"set\.seed\s*\(\s*([^)]+)\s*\)")
for folder in (ROOT / "scripts", ROOT / "graph"):
    for path in sorted(folder.rglob("*.R")):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        rel = path.relative_to(ROOT).as_posix()
        for line_no, line in enumerate(lines, 1):
            hit = seed_pattern.search(line)
            if not hit:
                continue
            context = " ".join(x.strip() for x in lines[max(0, line_no-2):min(len(lines), line_no+1)])
            scope = "Script-level seed; intended use is not clear from the call site alone."
            if rel.startswith("graph/") or re.search(r"plot\(|layout_with|geom_jitter|position_jitter", context):
                scope = "Figure-generation script seed; generally stabilizes rendering/jitter/layout and is not a project-wide analysis seed."
            elif re.search(r"sample\(|runif\(|bootstrap|permutation|AUCell_buildRankings", context, flags=re.I):
                scope = "Stochastic analysis/scoring seed; applies to this script or function invocation only."
            seed_rows.append({"source_file": rel, "line": str(line_no), "seed_expression": hit.group(1).strip(),
                              "source_context": context, "scope_note": scope})
add("Random seeds", "Project scripts", "explicit_set_seed_calls",
    f"See S9_random_seed_inventory.csv ({len(seed_rows)} explicit set.seed call sites); values are analysis-specific or figure-specific, not one global seed",
    "scripts/*.R; graph/*.R", status="reported in source",
    notes="Inventory includes source location, line, expression, and context; figure seeds are not analysis-wide seeds.")


def write_csv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


write_csv(OUT / "Table_S9_reproducibility_parameter_audit.csv", PARAM_COLS, records)
write_csv(OUT / "S9_random_seed_inventory.csv", SEED_COLS, seed_rows)

# Analysis-specific denominators and biological/statistical units.
COUNT_COLS = ["dataset", "analysis_scope", "species_tissue_model", "n_total_or_input",
              "group_counts", "statistical_unit", "n_used_or_rule", "source_file", "notes"]
count_rows = [
    dict(dataset="GSE292815", analysis_scope="12-week discovery / primary", species_tissue_model="Male C57BL/6J KC mouse pancreas", n_total_or_input="21", group_counts="CD=5; Safflower=5; HO-Safflower=4; Sunflower=4; HO-Sunflower=3", statistical_unit="Mouse/sample", n_used_or_rule="All 21 annotated samples; no exclusions", source_file="data/metadata/sample_metadata.tsv; results/manuscript_parameter_audit.md", notes="6-week start; 12-week feeding; collection at 18 weeks."),
    dict(dataset="GSE292815", analysis_scope="QC-only other timepoints", species_tissue_model="3-week KC and 8-week KCO mouse pancreas", n_total_or_input="18", group_counts="3-week KC=8; 8-week KCO=10", statistical_unit="Mouse/sample", n_used_or_rule="All available samples in S1; not part of 12-week integration", source_file="results/supplementary/S1_qc/S1_run_manifest.txt", notes="Only normalized GEO matrices available."),
    dict(dataset="GSE79668", analysis_scope="Independent PAAD immune-module replication", species_tissue_model="Human pancreatic adenocarcinoma bulk expression", n_total_or_input="51", group_counts="Not applicable", statistical_unit="Tumor sample", n_used_or_rule="n=51 for PRII and 8 immune modules; no usable OS time/event pair", source_file="results/external_validation/README.md; scripts/external_validation_gse79668.py", notes="Human symbols; scores standardized within cohort."),
    dict(dataset="GSE154778", analysis_scope="S4 single-cell context", species_tissue_model="Human PDAC scRNA-seq, primary and metastatic", n_total_or_input="14,926 cells; 16 samples", group_counts="Primary=10 samples/8,000 cells; metastatic=6 samples/6,926 cells", statistical_unit="Sample × cell type pseudobulk", n_used_or_rule="Pseudobulk units require >=10 cells; macrophage-reference family has six comparisons", source_file="results/supplementary/S4_single_cell_context/S4_run_manifest.txt", notes="Cells are not independent biological replicates."),
    dict(dataset="TCGA-PAAD", analysis_scope="Primary tumor universe", species_tissue_model="Human PAAD primary tumor RNA-seq", n_total_or_input="178", group_counts="Primary Tumor only", statistical_unit="Patient/sample with exact barcode", n_used_or_rule="178 in unadjusted signature analyses", source_file="results/immune_validation/tcga_sample_count_audit.csv; results/tcga_paad_external/analysis_manifest.txt", notes="No fatty-pancreas status measured."),
    dict(dataset="TCGA-PAAD", analysis_scope="Purity-adjusted immune association", species_tissue_model="Human PAAD primary tumor RNA-seq + ABSOLUTE purity", n_total_or_input="149", group_counts="One shared complete-case cohort", statistical_unit="Exact-barcode tumor sample", n_used_or_rule="PRII, immune scores, and purity all complete", source_file="results/immune_validation/tcga_cohort_manifest.md", notes="Not a sequential attrition count."),
    dict(dataset="TCGA-PAAD", analysis_scope="Survival", species_tissue_model="Human PAAD clinical survival", n_total_or_input="160; 93 deaths", group_counts="One survival complete-case cohort", statistical_unit="Patient", n_used_or_rule="Continuous Cox model; model-specific n reported", source_file="results/tcga_paad_external/analysis_manifest.txt; results/prii_consistency_audit.md", notes="Median-split KM descriptive only."),
    dict(dataset="TCGA-PAAD", analysis_scope="Adjusted survival models", species_tissue_model="Human PAAD clinical + covariates", n_total_or_input="132 age+stage; 111 age+stage+purity", group_counts="Different complete-case models", statistical_unit="Patient", n_used_or_rule="Do not combine these n values with the n=160 survival cohort", source_file="results/tcga_paad_external/tcga_sample_count_audit.csv", notes="Counts are analysis-specific, not sequential attrition."),
    dict(dataset="TCGA-PAAD", analysis_scope="EPIC convergence sensitivity", species_tissue_model="Human PAAD deconvolution", n_total_or_input="167/178 converged", group_counts="11 non-converged excluded", statistical_unit="Tumor sample", n_used_or_rule="Converged samples only for EPIC sensitivity", source_file="results/supplementary/S5_tcga_sensitivity/S5_run_manifest.txt", notes="Other deconvolution methods have method-specific missingness."),
    dict(dataset="GSE234956", analysis_scope="External metabolic context", species_tissue_model="WT mouse islets", n_total_or_input="45", group_counts="44% fat+fiber vs 10% fat+fiber; detailed group n in source metadata", statistical_unit="Mouse/islet sample", n_used_or_rule="Adjusted model includes strain and sex", source_file="results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", notes="Exposure duration unavailable in public metadata."),
    dict(dataset="GSE190812", analysis_scope="External metabolic context", species_tissue_model="Mouse pooled islets", n_total_or_input="21 pooled libraries", group_counts="RD=7; HFD=7; HFD-O304=7", statistical_unit="Pooled library (5 mice per library)", n_used_or_rule="Do not interpret 7 libraries as 35 independent mice", source_file="results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", notes="Nominal Mann–Whitney tests."),
    dict(dataset="GSE249790", analysis_scope="External metabolic context", species_tissue_model="Sorted mouse beta cells", n_total_or_input="7", group_counts="CD=4; 2-day HFD=3", statistical_unit="Mouse/sample", n_used_or_rule="Small exploratory cohort", source_file="results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", notes="Several core genes have zero/low counts."),
    dict(dataset="GSE293790", analysis_scope="External metabolic context", species_tissue_model="Female KC/Kras mouse whole pancreas", n_total_or_input="9", group_counts="CD=3; DIO=3; DIO→CD=3", statistical_unit="Mouse/sample", n_used_or_rule="Cancer-susceptible model; no pooled estimate", source_file="results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", notes="Not normal fatty pancreas."),
]
write_csv(OUT / "S9_cohort_sample_count_audit.csv", COUNT_COLS, count_rows)

# FDR family dictionary; every row is an explicit, analysis-specific family.
FDR_COLS = ["family_id", "module", "dataset", "hypothesis_family", "family_size_or_rule", "correction", "source_file", "status", "notes"]
fdr_rows = [
    ["FDR01", "Primary DE", PRIMARY, "All tested genes within each contrast", "One family per each of 9 contrasts", "BH", "scripts/06_12wk_integrated_analysis.R", "reported", "DEG display also requires |logFC| >= 1."],
    ["FDR02", "Primary enrichment", PRIMARY, "Pathways from moderated-t FGSEA within each contrast", "One pathway family per contrast", "BH", "scripts/06_12wk_integrated_analysis.R", "reported", "minSize=3; maxSize=2000."],
    ["FDR03", "Primary enrichment", PRIMARY, "ORA terms", "Separate contrast × direction × ontology/database", "BH", "scripts/06_12wk_integrated_analysis.R", "reported", "Requires >=5 directional DEGs; expressed genes used as universe."],
    ["FDR04", "Primary score contrasts", PRIMARY, "Five pathway/core scores per planned contrast", "2 contrast-specific families × 5 scores", "BH", "scripts/05_psps_score.R", "reported", "Do not call one contrast's family the global project FDR."],
    ["FDR05", "Primary score correlations", PRIMARY, "Pairwise correlations among five scores", "10 tests", "BH", "scripts/05_psps_score.R", "reported", "One family."],
    ["FDR06", "Primary core-gene contrasts", PRIMARY, "Six core genes per planned contrast", "2 contrast-specific families × 6 genes", "BH", "scripts/05_psps_score.R", "reported", "Core gene tests."],
    ["FDR07", "Mouse immune association", PRIMARY, "PRII versus non-overlapping immune modules", "9 module tests", "BH", "scripts/09_12wk_immune_microenvironment.R", "reported", "Separate from non-core robustness module families."],
    ["FDR08", "Mouse immune group contrasts", PRIMARY, "Nine module scores per planned contrast", "2 contrast-specific families × 9 modules", "limma BH", "scripts/09_12wk_immune_microenvironment.R", "reported", "Two planned matched high-oleic comparisons."],
    ["FDR09", "External human PAAD", "GSE79668", "PRII association with immune modules", "8 modules", "BH", "scripts/external_validation_gse79668.py", "reported", "n=51."],
    ["FDR10", "TCGA score association", "TCGA-PAAD", "Signature-score associations", "11 signature scores", "BH", "scripts/08_tcga_paad_external_validation.R", "reported", "Unadjusted TCGA association family."],
    ["FDR11", "TCGA immune association", "TCGA-PAAD", "PRII associations with seven prespecified immune modules", "7 modules within each analysis (unadjusted and purity-adjusted)", "BH by analysis", "scripts/09_immune_microenvironment_validation.R", "reported", "Each analysis has its own seven-module family."],
    ["FDR12", "TCGA method comparison", "TCGA-PAAD", "Macrophage estimates across methods", "35 comparisons", "BH", "results/supplementary/S5_tcga_sensitivity/S5_run_manifest.txt", "reported", "Metric comparisons."],
    ["FDR13", "Single-cell context", "GSE154778", "Macrophage versus reference cell types", "Each six-comparison model family", "BH", "results/supplementary/S4_single_cell_context/S4_run_manifest.txt", "reported", "Pseudobulk/sample-level inference."],
    ["FDR14", "GSE234956", "GSE234956", "Adjusted feature model", "14 features", "BH", "scripts/external_validation_gse234956.py; results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", "reported", "Per-strain tests are nominal."],
    ["FDR15", "GSE190812", "GSE190812", "Dietary arm comparisons", "Per manifest; not adjusted", "None reported", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", "reported", "Nominal Mann–Whitney tests."],
    ["FDR16", "GSE249790", "GSE249790", "PRII + six core genes", "7 rows", "BH", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", "reported", "Two-day HFD sorted beta-cell cohort."],
    ["FDR17", "GSE293790", "GSE293790", "Pairwise group tests", "Per contrast; no BH reported", "None reported", "scripts/external_validation_gse293790.py; results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt", "reported", "n=3/group."],
    ["FDR18", "Exploratory lipid", PRIMARY, "Four pooled correlations + eight diet contrasts", "12 prespecified tests", "Global BH", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "The q values used for the main S7 figure."],
    ["FDR19", "Exploratory lipid subscores", PRIMARY, "Four PRII subscores × four lipid programs, pooled", "16 tests", "BH", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "Separate family."],
    ["FDR20", "Exploratory lipid subscores", PRIMARY, "Four PRII subscores × four lipid programs, group-adjusted", "16 tests", "BH", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "This is subscore analysis, not full six-gene PRII."],
    ["FDR21", "Exploratory lipid", PRIMARY, "Within-group diagnostic tests", "80 tests; group n=3–5", "BH", "results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "Diagnostic only; not plotted or used as stable inference."],
    ["FDR22", "TF network", "GSE292815; TCGA-PAAD", "12 prespecified TFs", "12 tests per declared analysis family", "BH", "results/network_tf/analysis_manifest.txt; results/supplementary/S8_network_evidence_audit/S8_run_manifest.txt", "reported", "Computational association, not TF binding."],
    ["FDR23", "Robustness contrasts", PRIMARY, "Leave-one-gene-out score contrasts", "7 score rows per contrast (full PRII + 6 leave-one-gene-out scores)", "BH within contrast", "scripts/14_psps_robustness.R", "reported", "The full six-gene score is included in the family."],
    ["FDR24", "Robustness contrasts", PRIMARY, "Non-core mean-z score contrasts", "3 non-core scores per contrast", "BH within contrast", "scripts/14_psps_robustness.R", "reported", "Separate from LOO and ssGSEA families."],
    ["FDR25", "Robustness contrasts", PRIMARY, "ssGSEA score contrasts", "3 ssGSEA scores per contrast", "BH within contrast", "scripts/14_psps_robustness.R", "reported", "Separate from LOO and non-core mean-z families."],
    ["FDR26", "Robustness correlations", PRIMARY, "LOO scores versus non-core inflammation", "7 score correlations", "BH", "scripts/14_psps_robustness.R", "reported", "One family."],
    ["FDR27", "Robustness correlations", PRIMARY, "Six mean-z/ssGSEA scores versus non-core inflammation", "6 correlations", "BH", "scripts/14_psps_robustness.R", "reported", "One family."],
    ["FDR28", "Robustness correlations", PRIMARY, "Full PRII versus three non-core modules", "3 correlations", "BH", "scripts/14_psps_robustness.R", "reported", "Separate family from LOO and scoring-method agreement."],
    ["FDR29", "PRII subscores", PRIMARY, "Pyro/Inf mean-z and ssGSEA contrasts", "4 score rows per planned contrast", "BH within contrast", "scripts/21_psps_subscores.R", "reported", "Does not replace the frozen six-gene PRII primary score."],
    ["FDR30", "PRII subscores", PRIMARY, "Subscore pooled lipid-program correlations", "4 subscores × 4 lipid sets = 16 tests", "BH", "scripts/21_psps_subscores.R; results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "Same family also summarized as S7 pooled subscore family."],
    ["FDR31", "PRII subscores", PRIMARY, "Subscore within-group lipid correlations", "4 subscores × 4 sets × 5 groups = 80 tests", "BH", "scripts/21_psps_subscores.R; results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "Diagnostic only; group n=3–5."],
    ["FDR32", "PRII subscores", PRIMARY, "Subscore group-adjusted partial lipid correlations", "4 subscores × 4 lipid sets = 16 tests", "BH", "scripts/21_psps_subscores.R; results/supplementary/S7_exploratory_lipid_context/S7_run_manifest.txt", "reported", "Separate family from pooled and within-group tests."],
    ["FDR33", "PRII subscores", PRIMARY, "LOO scores versus four lipid programs", "7 scores × 4 lipid sets = 28 tests", "BH", "scripts/21_psps_subscores.R", "reported", "One family."],
    ["FDR34", "PRII subscores", PRIMARY, "Each subscore/method versus immune modules", "One 9-module family per mean-z/ssGSEA × Pyro/Inf subscore (4 families)", "BH within subscore family", "scripts/21_psps_subscores.R", "reported", "The nine immune modules are core-excluded marker scores."],
    ["FDR35", "PRII subscores", PRIMARY, "LOO scores versus immune modules", "7 LOO/full scores × 9 modules = 63 tests", "BH", "scripts/21_psps_subscores.R", "reported", "One family."],
    ["FDR36", "PRII subscores", PRIMARY, "Conditional Pyro/Inf versus immune-module correlations", "18 tests per scoring method (9 modules × 2 axes); 2 separate method families", "BH within method family", "scripts/21_psps_subscores.R", "reported", "Rank-residual partial correlations controlling the opposing subscore."],
    ["FDR37", "GSE154778 sample-level", "GSE154778", "Cell fraction association with overall PRII and same-cell-type PRII", "Finite cell-type rows per outcome family", "BH by outcome", "scripts/09_immune_microenvironment_validation.R", "reported", "Macrophage focus uses >=5 macrophage cells and is separately reported."],
]
write_csv(OUT / "S9_fdr_family_dictionary.csv", FDR_COLS,
          [dict(zip(FDR_COLS, row)) for row in fdr_rows])

# Capture recorded environments without collapsing module-specific R versions.
SOFT_COLS = ["module", "source_file", "R_version", "platform", "package_versions_captured", "scope_note"]
def module_for(rel: str) -> str:
    if "/S1_qc/" in rel: return "S1"
    if "/S3_prii_robustness/" in rel: return "S3"
    if "/S4_single_cell_context/" in rel: return "S4"
    if "/S5_tcga_sensitivity/" in rel: return "S5"
    if "/S6_external_metabolic_context/" in rel: return "S6"
    if "/S7_exploratory_lipid_context/" in rel: return "S7"
    if "/S8_network_evidence_audit/" in rel: return "S8"
    if "/12wk_primary/" in rel: return "Primary 12-week"
    if "/tcga_paad_external/" in rel: return "TCGA-PAAD"
    if "/immune_validation/" in rel: return "TCGA immune validation"
    if "/network_tf/" in rel: return "TF network"
    if "/network_ppi/" in rel: return "STRING PPI"
    if "/network_ligand_receptor/" in rel: return "Ligand–receptor"
    if "/sc_location_GSE154778/" in rel: return "GSE154778 upgrade"
    if "/composition_lipid_integration/" in rel: return "Composition/lipid integration"
    if "/supplementary_S2/" in rel: return "S2 source analysis"
    return "Other"

soft_rows = []
for path in sorted((ROOT / "results").rglob("sessionInfo*.txt")):
    text = path.read_text(encoding="utf-8", errors="replace")
    rel = path.relative_to(ROOT).as_posix()
    first = text.splitlines()[0] if text.splitlines() else ""
    rver = first.strip() if first.startswith("R version") else "not captured"
    platform_match = re.search(r"^Platform:\s*(.+)$", text, flags=re.M)
    sections = text
    if "other attached packages:" in text:
        sections = text.split("other attached packages:", 1)[1]
    package_versions = sorted(set(re.findall(r"\b([A-Za-z][A-Za-z0-9.]*)_([0-9][A-Za-z0-9.+-]*)", sections)))
    compact = "; ".join(f"{name} {version}" for name, version in package_versions)
    note = "Recorded sessionInfo; version applies to this specific run only."
    if "S4_single_cell_context" in rel:
        note += " S4-specific session does not capture UCell version."
    if "S7_exploratory_lipid_context" in rel:
        note += " S7-specific session does not establish exact GSVA/msigdbr runtime versions."
    if "S6_external_metabolic_context" in rel or "S8_network_evidence_audit" in rel:
        note += " Figure/audit session may capture plotting packages rather than all upstream analysis packages."
    soft_rows.append({"module": module_for(rel), "source_file": rel, "R_version": rver,
                      "platform": platform_match.group(1).strip() if platform_match else "not captured",
                      "package_versions_captured": compact, "scope_note": note})
write_csv(OUT / "S9_software_environment_audit.csv", SOFT_COLS, soft_rows)

# Gaps are explicit statuses, not imputed parameters.
GAP_COLS = ["gap_id", "area", "item", "status", "evidence_or_source", "submission_handling"]
gaps = [
    ["G01", "Animal study", "Randomization / blocking / cage or litter balancing", "not reported", "results/manuscript_parameter_audit.md", "State that the source did not report allocation procedures; do not say randomized."],
    ["G02", "Raw sequencing QC", "FASTQ-level QC and raw sequencing-depth metrics", "not available", "results/supplementary/S1_qc/S1_run_manifest.txt", "Describe QC as normalized-matrix-level only."],
    ["G03", "GSE190812 / GSE249790", "Tracked score-generation scripts and exact preprocessing", "not recoverable", "results/supplementary/S6_external_metabolic_context/S6_run_manifest.txt; scripts/", "Use only documented manifests and label pipeline detail as unavailable."],
    ["G04", "S4 single cell", "Exact UCell version for the S4 runtime", "not captured", "results/supplementary/S4_single_cell_context/S4_sessionInfo.txt", "Do not substitute the version from an earlier upgrade session."],
    ["G05", "S7 lipid scoring", "Exact GSVA/msigdbr version in S7 runtime", "not captured", "results/supplementary/S7_exploratory_lipid_context/S7_sessionInfo.txt", "Report manifest as GSVA 2.x; do not claim the older runtime is exact."],
    ["G06", "Python analyses", "Python and dependency versions / lockfile", "not reported", "scripts/external_validation_gse79668.py; scripts/external_validation_gse234956.py; results/network_ligand_receptor/run_manifest.txt", "List package versions as unrecorded; preserve accession and source inputs."],
    ["G07", "Network databases", "Exact STRING, OmniPath, and DoRothEA snapshot/release", "not reported", "results/network_ppi/README.md; results/network_tf/analysis_manifest.txt; results/network_ligand_receptor/run_manifest.txt", "State database source and thresholds; avoid claiming a release snapshot."],
    ["G08", "Incremental value", "VIF and overlap-adjusted partial-correlation outputs", "not estimated", "results/prii_supplementary_tables_plan_v1_zh_en.md", "Do not infer independent incremental value from current analyses."],
    ["G09", "GSE79668", "Survival validation", "not estimated", "results/external_validation/README.md", "No usable OS time/event pair; do not present external survival replication."],
    ["G10", "S5 survival", "delta_LRT, delta_AIC, C-index", "not estimated", "results/prii_supplementary_tables_plan_v1_zh_en.md", "Do not fill with inferred values."],
    ["G11", "S10", "Incremental signature analysis", "pending", "results/prii_incremental/statistical_analysis_plan_zh_en.md", "Describe as planned only until model and external replication outputs exist."],
    ["G12", "Global software environment", "One project-wide R version", "not applicable", "results/**/sessionInfo*.txt", "R versions are module-specific (R 4.4.0 and R 4.6.1); report by module."],
]
write_csv(OUT / "S9_unreported_and_unestimated_items.csv", GAP_COLS,
          [dict(zip(GAP_COLS, g)) for g in gaps])

readme = """# Supplementary Table S9 — reproducibility and parameter audit

## 中文说明

本包汇总论文分析涉及的实验设计来源、队列样本口径、PRII 定义、统计参数、FDR 家族、随机种子、软件记录和未报告/未估计事项。它是对现有脚本、运行清单、sessionInfo 和结果审计文件的索引式核验，不重新运行分析，也不补造缺失参数。

机器可读主文件为 `Table_S9_reproducibility_parameter_audit.csv`，一行对应一个参数或状态；配套工作簿 `Table_S9_reproducibility_parameter_audit.xlsx` 便于人工筛选审阅。其余 CSV 分别列出分析队列分母、BH/FDR 家族、源码中显式 `set.seed` 调用、各运行记录的软件环境，以及明确的缺口。`S9_reproducibility_checklist_zh_en.md` 可作为投稿前核对单；CSV 是规范数据源，工作簿只呈现主表。

必须保留的解释边界：动物饮食分配是否随机未报告；TCGA 的 178/149/160/132/111 是不同分析的完整病例口径，不能作为串行筛除数字；R 4.4.0 与 R 4.6.1 按模块分别记录；MSigDB `2026.1.Hs` 只对应四个经审计的外部脂质基因集；GSE190812 和 GSE249790 的评分生成脚本/确切预处理不可从当前项目恢复；S10 增量价值分析仍是计划，不能写成已完成结果。

## English note

This package indexes the reported study design, cohort denominators, PRII definition, analysis parameters, FDR families, random seeds, recorded software environments, and explicitly unreported or unestimated items. It audits existing scripts, manifests, sessionInfo files, and consistency reports; it does not rerun analyses or impute missing parameters.

The primary machine-readable file is `Table_S9_reproducibility_parameter_audit.csv`, with one row per parameter/status. Companion CSVs provide analysis-specific sample counts, BH/FDR family definitions, explicit `set.seed` call sites in source code, captured software environments, and unresolved items. The checklist is suitable for submission review. Preserve these boundaries: allocation randomization was not reported; TCGA n=178/149/160/132/111 refer to different analysis-specific complete-case cohorts; R 4.4.0 and R 4.6.1 are module-specific; MSigDB `2026.1.Hs` applies only to the four audited external lipid sets; score-generation scripts/exact preprocessing for GSE190812 and GSE249790 are not recoverable; S10 remains planned rather than completed.

## Table title

中文：补充表 S9．可重复性、参数和完整性审计
English: Supplementary Table S9. Reproducibility, parameter, and completeness audit

## Column/status conventions

`reported` = explicitly present in a script or run record; `partially reported` = the broad method is documented but an exact version or preprocessing source is absent; `not reported` = searched source artifacts do not document the item; `not estimated` = no result was calculated; `not recoverable` = the producing pipeline is not present in the project; `pending` = planned work without completed result files. Blank fields mean not applicable or not individually stated; they are not zero.
"""
(OUT / "README_zh_en.md").write_text(readme, encoding="utf-8")

checklist = """# S9 reproducibility checklist / 可重复性核对单

| Check | Status / handling |
|---|---|
| Discovery animals and feeding design | Male C57BL/6J KC; start 6 weeks; 12-week feeding; harvest 18 weeks; n=21 and group counts verified. |
| Allocation procedure | Randomization/blocking/cage balancing not reported; do not claim random assignment. |
| Diet wording | CD 18% kcal fat; matched HFD arms 45% kcal fat; state fat-source contrasts precisely. |
| PRII definition | Fixed six genes; cohort-wise per-gene z using sample SD, then arithmetic mean; legacy PSPS_core is an alias. |
| Cross-cohort scores | Standardized within each cohort; do not compare absolute PRII values across cohorts. |
| Sample denominators | Keep cohort- and model-specific n; TCGA n values are separate complete-case analyses. |
| FDR | Cite the family dictionary; do not call module-specific BH corrections one global FDR. |
| Random seeds | Use per-analysis values and source context in `S9_random_seed_inventory.csv`; no single project-wide seed exists. |
| Single-cell inference | Sample/pseudobulk is the inferential unit; cell counts are not biological replicates. |
| Software | Cite module-specific sessionInfo; R 4.4.0 and R 4.6.1 both occur. |
| External metabolic cohorts | Keep tissues/models separate; no cross-cohort pooling/meta-analysis. |
| Database networks | Separate expression association, database prior, and candidate ligand–receptor expression support. |
| Missing items | Preserve `not reported`, `not estimated`, `not recoverable`, and `pending` as distinct statuses. |
| Incremental value | S10 is pending; do not claim independent added value before completing its planned analyses. |
"""
(OUT / "S9_reproducibility_checklist_zh_en.md").write_text(checklist, encoding="utf-8")

manifest = [
    "Supplementary Table S9: reproducibility, parameter, and completeness audit",
    f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
    "Build: scripts/45_s9_reproducibility_audit.py",
    f"Parameter/status rows: {len(records)}; cohort denominator rows: {len(count_rows)}; FDR family rows: {len(fdr_rows)}; sessionInfo records: {len(soft_rows)}; explicit gaps: {len(gaps)}; explicit set.seed call sites: {len(seed_rows)}",
    "This is a metadata compilation; no expression analysis or statistical result was recalculated.",
    "Primary output: Table_S9_reproducibility_parameter_audit.csv; companion outputs are machine-readable CSVs and bilingual Markdown.",
    "Known boundaries: random allocation not reported; module-specific software versions; S10 pending; no imputed missing parameters.",
]
(OUT / "S9_run_manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")

print("S9 audit files written to", OUT)
print(f"Parameter rows={len(records)}; FDR families={len(fdr_rows)}; sessionInfo records={len(soft_rows)}; gaps={len(gaps)}")
