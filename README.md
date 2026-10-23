# Pancreatic steatosis and pyroptosis: analysis code

This public repository contains source scripts used for the PRII transcriptomic analyses and the BMC Genomics figure layouts. PRII is the fixed six-gene transcriptomic proxy defined in the manuscript; the code does not measure functional pyroptosis.

## Data

Raw and processed expression matrices, TCGA sample-level data, large derived result matrices, rendered manuscript files, and figure exports are not redistributed here. The analyses use public GEO datasets GSE292815, GSE154778, GSE79668, GSE293790, GSE234956, GSE190812, and GSE249790, plus TCGA-PAAD through the NCI Genomic Data Commons. Obtain data from the original repositories and follow their terms of use. `data/metadata/matrix_sample_map.tsv` records the GSE292815 matrix-to-sample mapping used by the discovery scripts.

## Code layout

- `scripts/`: analysis scripts for gene-set construction, differential expression, scoring, immune context, external cohorts, deconvolution, survival sensitivity, and network analyses.
- `graph/`: final figure-source scripts for the main and supplementary figures.
- `submission_bmc_genomics/scripts/`: journal-layout figure scripts.

Several scripts consume intermediate files produced by earlier analyses and expect the original project directory layout (for example, `data/raw/`, `data/external/`, and `results/`). This repository is a source-code release, not a self-contained one-command workflow; a complete run requires downloading/preparing the inputs and running the analysis stages in dependency order. No raw cohort-level data or private credentials are included.

## Software

The scripts use R and Python packages declared in their imports. Package requirements vary by analysis; the release does not yet include a single lockfile. Check each script and the associated published methods before running it. Statistical outputs should be regenerated from the public inputs and checked against the manuscript and its supplementary files.

## Reproducibility scope

The code is provided to document and support the reported analyses. Some external-cohort scripts use processed matrices or project-derived intermediate files, and the original preprocessing for those cohorts may not be recoverable. Do not interpret those scripts as establishing a uniform pipeline across datasets.
