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
- Snakemake (>=7.0)
- Conda or Mamba (for environment management)
- For GPU-accelerated basecalling: CUDA-compatible GPU

## Setup the environment for code and data sharing
1. Clone this repository:
```bash
git clone https://github.com/cnio-bu/tcca.git
cd tcca
```
2. Download the preprocessed TCCA dataset, which contains 1,089,024 malignant cells and 750,218 non-malignant cells from 36 studies, along with metadata:
https://tcca.bioinfo.cnio.es/session/8609a42b7abfbebe6b39778186468dd4/download/download-dl_h5ad?w=

## Workflow steps
### 1. GEO data preparation
Before running the Snakemake workflow, raw GEO single-cell data must be converted to the standard 10x Genomics format (matrix.mtx, genes.tsv, barcodes.tsv). We provide two helper Python scripts for this step:
- `src/utils/sort_unzip_geo_v2.py` — for datasets using the new 10x format (features.tsv.gz).
- `src/utils/sort_unzip_geo_v3.py` — for datasets using the legacy format (genes.tsv.gz).
  
These scripts automatically:
- Detect sample identifiers (GSMxxxxx) in a directory.
- Create a subfolder for each sample.
- Unzip and rename the matrix, gene/feature, and barcode files.
<br>

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

<br>

### 3. Seurat object conversion and structuring
After preprocessing, Seurat v4 objects from each study are converted to optimized Seurat v5 formats using BPCells for efficient storage and analysis with scripts in `src/seurat`.
- `export_all_raw_to_v5.R`: converts Seurat v4 objects per study into Seurat v5, retains only protein-coding genes, and saves study-level Seurat v5 objects backed by BPCells matrices.
- `generate_all_levels_seurat.R`: builds merged Seurat v5 objects across all studies at three filtering levels (lvl3: all cells; lvl2: clinically annotated samples with ≥100 malignant cells; lvl1: malignant-only subset from lvl2).
> These scripts are run after the main Snakemake pipeline.

<br>

### 4. Single-cell annotation and TME classification
This directory `src/cell_annotation_tme` contains the R notebook detailing the harmonization of cell type annotations. It covers the standardization of cell names, integration of author-provided labels with automated annotations (Azimuth and SingleR), calculation of sample-specific immune and stromal proportions, and the definition of 12 Tumor Microenvironment (TME) archetypes through clustering.
**Scripts**
- `add_author_annotation.R`: adds missing or mismatched author-provided cell type labels to preprocessed Seurat objects, ensuring correct mapping of cell names and updating metadata.
- `rename_tme_archetype.R`: rename TME archetypes
  
<br>

### 5. Basic single-cell analysis with Scanpy
The following scripts (`src/scanpy`) are run in Python using Scanpy based on the Seurat v5 level 2 (lvl2) dataset exported as an .h5ad file.
- `envs/` directory provides environment specifications for running integration workflows and benchmarking analyses.
- `dimred_before_int.py`: normalization, HVG detection, PCA, and initial clustering before integration.
- `integration.py`: applies multiple state-of-the-art integration methods (e.g., scVI, SCANVI, Scanorama) o generate batch-corrected embeddings for integrated visualization across studies.
- `benchmark_integration.py`: evaluates integration performance using standard metrics to compare batch correction and biological conservation across the three integration methods.,
- `cluster.py`: visualizes integrated embeddings, performs clustering (all cells and malignant subset), and validates cell type annotation using known marker genes.

<br>

### 6. Beyondcell results analysis
After computing Beyondcell scores in the snakemake pipeline, several scripts in `src/beyondcell` integrate and analyze drug-response profiles across studies:
- `export_bcimmuno_to_BP.R`: converts individual Beyondcell results into Seurat v5 objects with BPCells matrices (cells × drugs).
- `merge_studywise_bps_immuno.R`: aligns and merges all Beyondcell matrices into a unified matrix with consistent drug identifiers and metadata.
- `therapeutic_clusters.R`: performs basic analysis on the global Beyondcell matrix, including dimensionality reduction, clustering of cells based on drugs scores using Seurat v5 and generates sketch datasets.
- `optimal_clustering.R`: evaluates clustering parameters (k, resolution) using internal validation metrics (Davies–Bouldin, Silhouette, Purity).
- `therapeutic_clusters_umap_heatmap.R`: visualizes final therapeutic clusters, drug activity patterns, and their distribution across cancer types.
- `get_sp_bc.R`: extract the switch point, a metric summarizing the balance between sensitive and resistant cell populations per drug and sample, allowing easier drug response comparison across samples.
- `bcgeneset_to_gmt.R`: converts Beyondcell gene set collection (SSc) to GMT format and merges it with the immunotherapy GMT file.
> All scripts are run outside the Snakemake pipeline and rely on Beyondcell outputs generated per study.

<br>

### 7. Inference of CNVs with SCEVAN
After completing the Snakemake pipeline, several standalone R scripts are used to process and integrate SCEVAN copy number alterations (CNVs) inference into interpretable, analysis-ready formats.
> These scripts are executed outside the Snakemake pipeline (`src/scevan/`).
> Many of them are designed for parallel computing on an HPC environment and are submitted as independent jobs using SLURM (sbatch), as indicated in the script headers.

**1. Running SCEVAN and preprocessing results**
- `1_scevan_all_samples.R`: run SCEVAN for all studies.
- `2_scevan_metadata.R`: extract clonality metrics per study.
- `3_combine_clonality_metadata.R`: merge all clonality outputs into a single summary table.
- `4_scevan_cellxgene.R`: generate per-gene CNV matrices for each study.
- `5_filter_cellxgene_tables.R`: filter CNV matrices to level 1 (only high-quality malignant cells and genes).
- `6_tsv_to_bpcells_lvl*.R`: convert CNV matrices into on-disk BPCells objects for efficient storage and access.
- `7_final_cellxgene_v5_cnafill.R`: combine all BPCells matrices into a single, genome-wide CNV Seurat v5 object including clinical and clonality metadata.
- `8_sketch_rds.R`: generate downsampled Seurat “sketch” datasets (50k and 5k cells) for efficient visualization and exploratory analysis.

**2. Generating interpretable CNV Matrices**  
These scripts summarize the raw SCEVAN output into biologically meaningful genomic units.
- `9_interpretable_matrix/1_regionxclone.R`: build a complete region × clone matrix across all samples.
- `9_interpretable_matrix/2_genomicranges.R`: create mapping files to associate SCEVAN regions with either custom genomic bins or cytogenetic bands.
- `9_interpretable_matrix/3_customregions.R`: aggregate CNV values into averaged custom genomic regions.
- `9_interpretable_matrix/3_cytobands.R`: aggregate CNV values into standard cytobands for more interpretable genome-level visualization.

<br>

### 8. Therapeutic analysis with scTherapy

The `src/sctherapy` folder contains the workflow for predicting subclone-level drug responses using scTherapy; characterizing therapeutic clusters (TCs); deriving transcriptomic signatures; and correlating these signatures with patient survival.

**1. Subclone-level predictions**  
- `sctherapy_fullcohort.R`: predicts drug responses for malignant subclones (only from samples containing normal or TME cells).

**2. Visualization**  
- `heatmap_subclones.R`: heatmap of predicted drug responses per subclone.  
- `jaccard_heatmap_clustering.R`: spectral clustering of subclones based on Jaccard similarity.  
- `plot_clustering_umap.R`: UMAP embedding of subclones colored by TC.

**3. Cluster characterization**  
- `cluster_characterization.R` and `clinical_chisq.R`: generate stacked barplots of TCs across clinical variables (sex, age, tumor type, TME composition, study).  
- Statistical tests (Chi-squared or Fisher’s exact) to identify drugs significantly associated with TCs.  
- Bubble and UpSet plots showing drugs grouped by mechanism of action.  
- **Significant drugs for each therapeutic cluster (TC):** FDR ≤ 0.05 and affecting ≥50% of subclones within the TC.
- `plot_bc_immuno.R`: visualizes Beyondcell-predicted immunotherapy sensitivity for malignant cells from scTherapy-defined subclones and summarizes mean scores per therapeutic cluster (TC).
- `generate_drug_table.R`: generates a summary table of predicted drugs, their mechanisms of action (MoAs), TC assignments, and associated clinical variables at the subclone level.

**4. TC transcriptomic signatures**  
- `marker_sctherapy_clusters.R`: subset malignant cells from the Seurat level 2 expression object with scTherapy predictions; normalize and identify highly variable and differentially expressed (DE) genes per TC.  
- Dotplots of top genes; export therapeutic cluster-specific signatures in GMT format.

**5. Survival analysis (TCGA)**  
- `survival_marker_genes.R`: sompute ssGSEA scores for TC signatures.  
- Fit Cox proportional hazards models adjusted for cancer type, sex, age, and tumor stage/grade.  
- Analyses for progression-free interval (PFI) and overall survival (OS); outputs include tables and forest plots.
- `envs/` directory provides environment specifications for running ssGSEA and survival analysis.

<br>

### 9. Functional analysis in the TCCA Cohort
This section describes the workflow for computing functional pathway enrichment in malignant cells from the TCCA cohort using complementary approaches: **GSVA, UCell, Beyondcell, PROGENy, NMF, and Hotspot**. The workflow is organized into four main stages:


**1. Pseudobulk generation and GSVA enrichment (`pseudobulk/`)**
- `bulk_from_expr.R`: aggregates malignant cell expression per sample to create a pseudobulk expression matrix.
- `bulk_functional.R`: computes gene set enrichment using GSVA with predefined gene sets (from MSigDB and published cancer transcriptional programs) and performs biclustering of pathways and samples to reveal functional patterns.
- `calculate_overlaps.R`: calculates pairwise Jaccard indices between gene sets to quantify redundancy and flag overlapping sets for potential filtering.
  
**2. Single-cell functional enrichment (UCell and Beyondcell) (`single-cell/`)**
- `functional_cell_level.R`: computes cell-level functional enrichment using AddModuleScore with metaprograms from Gavish et al., Nature 2023.
- `functional_ucell.R`: computes UCell-based functional enrichment per cell using signatures from MSigDB and cancer-related transcriptional programs (`/reference/combined_gsets_functional.gmt`).
Equivalent enrichment matrices are also generated using Beyondcell scores (instead of UCell) via Snakemake:
- `export_bcfunctional_to_BP.R`: converts Seurat v4 objects into Seurat v5 with BPCell matrices.
- `merge_studywise_functional_bps.R`: merges Beyondcell matrices from all studies into a single functional matrix across all malignant cells.
- `functional_heat_sketchs_all.R` and `functional_heat_sketchs_sctherapy.R`: generate sketches of 5,000 cells from the UCell enrichment matrix (including all malignant cell population or subsetting cells from scTherapy-defined subclones) to create ComplexHeatmap visualizations of functional activity with `functional_clusters_heatmap.R` and `functional_heatmap.R`.
  
**3. NMF-based metaprogram discovery (`nmf/`)**
- `compute_nmf.R`: runs sample-wise NMF decomposition on malignant cells to infer recurrent transcriptional programs. Executed in parallel across studies using SLURM (`run_nmf_study&subclone.sh`).
- `generate_mps_sample.R`and `robust_nmf_sample.R`: aggregate and cluster NMF-derived metaprograms across samples using a robust custom approach adapted from Gavish et al., Nature 2023 and Kinker et al., Nat Genet 2020.
- `mps_ucell.R`: calculates UCell enrichment for metaprograms at single cell level, generates heatmaps for 5k-cell sketches, and summarizes metaprogram activity with dot and violin plots (including CIN70 signature analysis) per therapeutic cluster.
- `assign_to_mps.py`: adds UCell scores of cancer-state metaprograms to malignant cells, assigns each cell to its dominant metaprogram, and visualizes metaprogram activity on scANVI-MDE UMAPs.
- `integration_from_mps.py`: adds UCell scores of cancer-state metaprograms to malignant cells and generates UMAP embeddings based on these signatures to integrate cells across samples and studies.
  
**4. Alternative pathway inference methods (`other_functional_infer/`)**
- `tcca_decoupler.py`: infers pathway activity scores using PROGENy gene signatures and decoupler’s multivariate linear model (MLM), producing per-cell activity matrices for downstream analysis.
- `hotspot_tcca.py`: identifies gene co-expression modules in malignant cells using the Hotspot algorithm. It computes local and global autocorrelations based on scANVI embeddings, clusters genes into modules, and outputs module definitions and per-cell activity scores.
- `annotate_hotspot.R`: performs functional enrichment of Hotspot-derived modules using MSigDB and custom gene sets via FGSEA overrepresentation analysis (ORA) to identify overrepresented pathways.

<br>

### 10. Intra-tumoral heterogeneity (ITH) analysis
The folder `/src/heterogeneity/`contains scripts to quantify, compare, and visualize different dimensions of intra-tumor heterogeneity (ITH) across samples and cancer types.
- `ith_score_subclone&sample.R`: quantifies ITH at both subclone and sample levels by computing transcriptomic, genomic, and therapeutic heterogeneity scores.
- `cor_genomic_vs_therapeutic_ith.R`: integrates genomic and therapeutic ITH metrics per sample, and explores their correlation across cancer types.
- `cor_genomic_ith_vs_sp.R`: computes correlation between genomic ITH and drug switch points (SP) per sample, both pan-cancer and by cancer type (primary and metastatic).

<br>

### 11. Additional analyses
This section contains supplementary scripts for extended analyses beyond the main workflow. It includes: enhanced cell and clinical annotations, correlation analyses linking scTherapy predictions, Beyondcell scores, and functional modules, as well as focused use cases for TC4 and TC10 subclones to explore CNV patterns, gene signatures, and drug sensitivity.

**1. Additional cell annotations (`cell_annotation/`)**
- `reharmonize_tcca_metadata.R`: harmonizes level-2 TCCA metadata by integrating subclone membership, scTherapy therapeutic cluster assignments, and a refined cancer-type classification (46 → 34 types).
- `predict_sex.R`: infers missing sex annotations from single-cell expression profiles with `cellXY` R package and updates the Seurat metadata accordingly.

**2. Clinical database (`clinical_annotation/`)**
- `annotate_cancer_cell_lines.R`: adds standardized metadata to profiled cancer cell lines (sex, age, disease, lineage, collection site).
- `filter_raw_clinical_data.R`: cleans and standardizes sample raw clinical metadata by removing duplicates/misassignments and harmonizing tumor types and sites.

**3. Correlations with Beyondcell??**
- `compare_sctherapy_vs_beyondcell.R`: correlates Beyondcell scores with scTherapy predictions (converted to numeric values) at the subclone level. 
- `cor_bc_vs_mps.R`: computes Pearson correlations between 43 functional MPs and 589 Beyondcell drugs.

**4. Use case of TC4 (`tc4_usecase/`)**
- `check_cnv_overlap_tc4.R`: identifies recurrent CNV co-amplifications across chr3q cytobands in TC4 subclones.
- `get_tc4_peak_gmt.R`: retrieves gene sets located in cytobands recurrently amplified in TC4 subclones (mainly chr3q peaks) and intersects them with top TC4 marker genes.
- `gdsc_tc4.R`: correlates TC4 CNV-derived gene signatures with GDSC2 drug sensitivity profiles to identify compounds associated with TC4 amplifications.
- `go_enrichment.R`: performs GO enrichment analysis on genes within CNV regions recurrently amplified in TC4 subclones to identify biological processes associated with these genomic alterations.
- `cor_bc_vs_tc4sigs.R`: correlates enrichment of TC4 amplification gene sets with Beyondcell drug sensitivity scores in TC4 cells.

**5. Use case of TC10 (`tc10_usecase/`)**
- `gdsc_tc10.R`: correlates TC10 gene signature enrichment with GDSC drug sensitivity across cell lines and tumor types.
- `brca_clinical.R`: visualizes BRCA subclone distributions across age groups, sample type, TME archetypes, and tumor subtypes (bars colored by TC).
- `brca_expr_clusters.R`: performs integration and clustering analysis of BRCA patient single-cell expression data using Seurat.
- `brca_bc_plots.R`: plots Beyondcell drug sensitivity scores for TC10-predicted drugs in BRCA patient cells. 

<br>


### 12. Figures
This section includes scripts used to generate the main and supplementary figures of the study, covering clinical summaries, TME archetype visualizations, and clonality analyses. Other figures are produced directly within the corresponding analysis scripts.
- `TCCA_palette.R`: defines color palettes for clinical variables, therapeutic clusters, and MPS groups.
- `fig1_sankey_cohort.R`: sankey plot showing sample distribution across key clinical variables.
- `tab1_summary_clinical_table.R`: generates a summary table of clinical features per cancer type.
- `fig1F_2B_S13_summary_plots.R`: summarizes single-cell and scTherapy data, including malignant vs. TME proportions, TME archetype distributions, top drugs and MoAs per TC, and TC composition.
- `fig1E_circular_tme.R`: circular plot displaying sample percentages by TME group and archetype.
- `sankey_three_layers.R`: three-layer Sankey plot linking subclones, TME archetypes, tumor types, and functional metaprograms.
- `figS5_S6_clonality_plots.R`: computes and visualizes clonality (subclones per 1,000 cells) across cancer types and clinical variables.