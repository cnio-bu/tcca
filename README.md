# The Therapeutic Cancer Cell Atlas (TCCA)

## Overview
This repository contains all code and workflows used to preprocess, analyze, and visualize the **Therapeutic Cancer Cell Atlas (TCCA)** — an integrated single-cell atlas with therapeutic predictions built from publicly available datasets. It includes a Snakemake-based pipeline for preprocessing GEO single-cell expression data, and modular scripts for genomic, therapeutic, functional and tumor microenvironment analyses.

## Workflow overview
This is a simplified diagram of the workflow steps:
<p align="center">
    <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./.img/general_workflow_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="./.img/general_workflow.png">
    <img alt="Fallback image description" src="./.img/general_workflow.png", width = 70%>
    </picture>
</p>

## Requirements

## Workflow steps
### 1. GEO data preparation
Before running the Snakemake workflow, raw GEO single-cell data must be converted to the standard 10x Genomics format (matrix.mtx, genes.tsv, barcodes.tsv). We provide two helper Python scripts for this step:
- `src/utils/sort_unzip_geo_v2.py` — for datasets using the new 10x format (features.tsv.gz).
- `src/utils/sort_unzip_geo_v3.py` — for datasets using the legacy format (genes.tsv.gz).
  
These scripts automatically:
- Detect sample identifiers (GSMxxxxx) in a directory.
- Create a subfolder for each sample.
- Unzip and rename the matrix, gene/feature, and barcode files.


### 2. Snakemake pipeline
The Snakemake workflow orchestrates the full analysis across all studies.

#### 2.1. Preprocessing
First, the pipeline runs the preprocessing rules (`preprocessing.smk`) independently  for each study, since each study requires a different script depending on its input files. These scripts apply study-specific Seurat filtering steps to remove low-quality cells. Each rule generates a standardized Seurat v4 object saved as `{results}/seurat/raw/{study_name}.rds`. Each `.rds` file is a list of Seurat objects, one for each sample in the study.

#### 2.2. Malignant cells subsetting and Beyondcell prediction
Next, the outputs from preprocessing are used by the Beyondcell rules (`beyondcell.smk`). For each study, study-specific Beyondcell scripts are applied to filter malignant cells from the Seurat objects and computes single-cell drug sensitivity scores, producing:

- List of Seurat v4 objects containing only malignant cells for each study.
- List of Beyondcell objects with predicted drug sensitivity profiles for each study.
  
For more details, see the [Beyondcell repository](https://github.com/cnio-bu/beyondcell).

#### 2.3. Functional enrichment
To compute functional enrichment of malignant cells across different gene expression signatures, the pipeline includes functional enrichment rules for each study. These rules take as input the malignant-cell objects and a reference gene set file (e.g., `reference/combined_gsets_functional.gmt`), and calculate enrichment scores for each cell using the Beyondcell Score (BCS) through the script `general_functional_enrichment.R`.

### 3. From Seurat

### 3. Inference of CNVs with SCEVAN
After completing the Snakemake pipeline, several standalone R scripts are used to process and integrate SCEVAN copy number alterations (CNVs) inference into interpretable, analysis-ready formats.
> These scripts are executed outside the Snakemake pipeline (`src/scevan/`)
> Many of them are designed for parallel computing on an HPC environment and are submitted as independent jobs using SLURM (sbatch), as indicated in the script headers.

#### 3.1. Running SCEVAN and preprocessing results
- `1_scevan_all_samples.R`: run SCEVAN for all studies.
- `2_scevan_metadata.R`: extract clonality metrics per study.
- `3_combine_clonality_metadata.R`: merge all clonality outputs into a single summary table.
- `4_scevan_cellxgene.R`: generate per-gene CNV matrices for each study.
- `5_filter_cellxgene_tables.R`: filter CNV matrices to level 1 (only high-quality malignant cells and genes).
- `6_tsv_to_bpcells_lvl*.R`: convert CNV matrices into on-disk BPCells objects for efficient storage and access.
- `7_final_cellxgene_v5_cnafill.R`: combine all BPCells matrices into a single, genome-wide CNV Seurat v5 object including clinical and clonality metadata.
- `8_sketch_rds.R`: generate downsampled Seurat “sketch” datasets (50k and 5k cells) for efficient visualization and exploratory analysis.
#### 3.2. Generating Interpretable CNV Matrices
These scripts summarize the raw SCEVAN output into biologically meaningful genomic units.
- `9_interpretable_matrix/1_regionxclone.R`: build a complete region × clone matrix across all samples.
- `9_interpretable_matrix/2_genomicranges.R`: create mapping files to associate SCEVAN regions with either custom genomic bins or cytogenetic bands.
- `9_interpretable_matrix/3_customregions.R`: aggregate CNV values into averaged custom genomic regions.
- `9_interpretable_matrix/3_cytobands.R`: aggregate CNV values into standard cytobands for more interpretable genome-level visualization.


### 4. Single-cell annotation and TME classification
This directory (`src/cell_annotation_tme`) contains the R notebook detailing the harmonization of cell type annotations. We first standardized all cell names into a unified nomenclature, integrating author-provided labels with results from automated tools (Azimuth and SingleR) used on un-annotated studies. Finally, we calculated sample-specific proportions of immune and stromal populations to define 12 distinct Tumor Microenvironment (TME) archetypes through clustering.



