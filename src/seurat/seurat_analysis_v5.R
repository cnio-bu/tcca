library(BPCells)
library(Seurat)
library(tidyverse)

## set wd, required for bpcells relative path alloc.
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_cell_types")

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

seu <- readRDS(
    file = "merged_bc_metanalysis_full.Rds"
    )


raw_mat <- BPCells::open_matrix_dir(
    dir = "merged_bc_meta_counts"
    )

seu[["RNA"]]$counts <- raw_mat

print("LOADED")
## Import the database
clinical_metadata <- data.table::fread(
    "clinical_metadata_v2_clean.tsv"
    )

## Subset the full object to keep samples with malignant data
samples_to_keep <- intersect(unique(seu$sample), clinical_metadata$sample)

## Remove duplicated samples as identified by the annotation process
seu <- subset(seu, subset = sample %in% samples_to_keep)

## Normalize samples
#seu <- NormalizeData(seu)

## Fix a """duplicated""" sample called "t19" in two instances
meta.data <- seu@meta.data 
meta.data[c(1392:1762), "sample"] <- "T19_1"

clinical_metadata[
    clinical_metadata$sample == "T19" &
    clinical_metadata$study == "adrenalnb_rui_chong",
    "sample"
    ] <- "T19_1"

## Add clinical metadata
seu_cell_level_clinical <- seu@meta.data %>%
    rownames_to_column("cell") %>%
    left_join(
        y = clinical_metadata,
        by = "sample"
    )

rownames(seu_cell_level_clinical) <- seu_cell_level_clinical$cell
seu_cell_level_clinical$cell <- NULL

seu@meta.data <- seu_cell_level_clinical

## Split assay into n layers, each a sample
print("Splitting")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)

print("SPLITTED")
saveRDS(object=seu, file = "seu_split.Rds", destdir="/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_cell_types")
print("SAVED")
