# Pancreatic steatosis and pyroptosis: analysis code

This public repository contains source scripts used for the PRII transcriptomic analyses and the BMC Genomics figure layouts. PRII is the fixed six-gene transcriptomic proxy defined in the manuscript; the code does not measure functional pyroptosis.

## Data

Raw and processed expression matrices and TCGA sample-level clinical data are not redistributed. Derived figures, tables and numbered Additional files staged for the manuscript are available in the versioned v1.2.0 results release linked below. The analyses use public GEO datasets GSE292815, GSE154778, GSE79668, GSE293790, GSE234956, GSE190812, and GSE249790, plus TCGA-PAAD through the NCI Genomic Data Commons. Obtain data from the original repositories and follow their terms of use. `data/metadata/matrix_sample_map.tsv` records the GSE292815 matrix-to-sample mapping used by the discovery scripts.

## Code layout

- `scripts/`: analysis scripts for gene-set construction, differential expression, scoring, immune context, external cohorts, deconvolution, survival sensitivity, and network analyses.
- `graph/`: final figure-source scripts for the main and supplementary figures.
- `submission_bmc_genomics/scripts/`: journal-layout figure scripts.

Several scripts consume intermediate files produced by earlier analyses and expect the original project directory layout (for example, `data/raw/`, `data/external/`, and `results/`). This repository is a source-code release, not a self-contained one-command workflow; a complete run requires downloading/preparing the inputs and running the analysis stages in dependency order. No raw cohort-level data or private credentials are included.

## Software

The scripts use R and Python packages declared in their imports. Package requirements vary by analysis; the release does not yet include a single lockfile. Check each script and the associated published methods before running it. Statistical outputs should be regenerated from the public inputs and checked against the manuscript and its supplementary files.

## Reproducibility scope

The code is provided to document and support the reported analyses. Some external-cohort scripts use processed matrices or project-derived intermediate files, and the original preprocessing for those cohorts may not be recoverable. Do not interpret those scripts as establishing a uniform pipeline across datasets.

## Final supplementary analyses

- `scripts/47_s10_signature_audit.R`, `scripts/48_s10_incremental_models.R`, and `scripts/49_run_published_proxy_audit.R` reproduce the S10 comparator-overlap, conditional-association, nested-model, VIF, and published-proxy audit outputs staged with the manuscript.
- `scripts/50_mouse_immune_partial_spearman_by_diet.py` reproduces the exploratory diet-group-adjusted mouse immune-module sensitivity table (Table S11).

The S10 outputs are exploratory conditional bulk-transcriptomic associations. The Zuo score is a standardized directional proxy because the published expression scale was not fully reproducible; the Xie panel is descriptive only. The S10 and S11 scripts require the project input/results layout described above and do not download or redistribute raw data.

## Versioned submission results

The current source-code snapshot, including the clarification that the four-panel S8 redraw is an unselected alternative, is available in the [v1.2.1 code release](https://github.com/sundalong2022-create/pancreatic-steatosis-pyroptosis/releases/tag/v1.2.1). The derived figures, tables, Additional files, English legends, manifest and SHA-256 checksums staged for the BMC Genomics submission are published as the [v1.2.0 submission-results release](https://github.com/sundalong2022-create/pancreatic-steatosis-pyroptosis/releases/tag/v1.2.0). The downloadable package is [pancreatic-steatosis-pyroptosis-results-v1.2.0.zip](https://github.com/sundalong2022-create/pancreatic-steatosis-pyroptosis/releases/download/v1.2.0/pancreatic-steatosis-pyroptosis-results-v1.2.0.zip) with SHA-256 `828e25a48b39fe28253bc7911ba6024d251c0f601da19fdd49580b243be46ef9`. It contains derived submission outputs only; raw matrices and private credentials are not redistributed. Supplementary Figures S1–S8 are the original project TIFFs, S10 is the completed sensitivity figure, and S9 is a reproducibility supplementary table without a separate figure.
