# The Therapeutic Cancer Cell Atlas (TCCA)

## Overview
This repository contains all code and workflows used to preprocess, analyze, 
and visualize the **Therapeutic Cancer Cell Atlas (TCCA)** — an integrated 
single-cell atlas with therapeutic predictions built from publicly available datasets.  
It includes a Snakemake-based pipeline for preprocessing GEO single-cell expression data, 
and modular scripts for genomic, therapeutic, functional and tumor microenvironment analyses.


## Workflow overview
This is a simplified diagram of the workflow steps:
<p align="center">
    <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./.img/general_workflow_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="./.img/general_workflow.png">
    <img alt="Fallback image description" src="./.img/general_workflow.png", width = 50%>
    </picture>
</p>


## Requirements

## Workflow steps
### 1. GEO data preparation
Before running the Snakemake workflow, raw GEO single-cell data must be converted to 
the standard 10x Genomics format (matrix.mtx, genes.tsv, barcodes.tsv). We provide 
two helper Python scripts for this step:
- `sort_unzip_geo_v2.py` — for datasets using the new 10x format (features.tsv.gz).
- `sort_unzip_geo_v3.py` — for datasets using the legacy format (genes.tsv.gz).
  
These scripts automatically:
- Detect sample identifiers (GSMxxxxx) in a directory.
- Create a subfolder for each sample.
- Unzip and rename the matrix, gene/feature, and barcode files.

### 2. Snakemake pipeline
The Snakemake workflow orchestrates the full analysis across all studies.

##### Preprocessing
First, the pipeline runs the preprocessing rules (`preprocessing.smk`) independently 
for each study, since each study requires a different script depending on its input 
files. These scripts apply study-specific Seurat filtering steps to remove low-quality 
cells. Each rule generates a standardized Seurat v4 object saved as `{results}/seurat/raw/{study_name}.rds`. 
Each `.rds` file is a list of Seurat objects, one for each sample in the study.

##### Malignant cells subsetting and Beyondcell prediction
Next, the outputs from preprocessing are used by the Beyondcell rules (`beyondcell.smk`). 
For each study, study-specific Beyondcell scripts are applied to filter malignant cells 
from the Seurat objects and computes single-cell drug sensitivity scores, producing:

- List of Seurat v4 objects containing only malignant cells for each study.
- List of Beyondcell objects with predicted drug sensitivity profiles for each study.
  
For more details, see the [Beyondcell repository](https://github.com/cnio-bu/beyondcell).

##### Functional enrichment
To compute functional enrichment of malignant cells across different gene expression 
signatures, the pipeline includes functional enrichment rules for each study. 
These rules take as input the malignant-cell objects and a reference gene set file 
(e.g., `reference/combined_gsets_functional.gmt`), and calculate enrichment scores 
for each cell using the Beyondcell Score (BCS) through the script `general_functional_enrichment.R`.

##### 







