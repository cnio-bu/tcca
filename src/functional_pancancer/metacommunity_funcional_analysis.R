library(BPCells)
library(Seurat)
library(tidyverse)

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

functional_mat <- open_matrix_dir(
    dir = "results/functional/full_mat_functional/"
    )

meta.data <- read_tsv(
    "results/annotation/functional_metadata_with_clinical.tsv"
    )

meta.data_full_clinical <- meta.data %>%
    as.data.frame()

rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

## remove empty rows
functional_mat <- functional_mat[1:130, ]

fc <- CreateSeuratObject(
    counts = functional_mat,
    assay = "RNA",
    project = "fc_pancancer",
    meta.data = meta.data_full_clinical
)

## load metacom mat
mcs <- readRDS("results/beyondcell_bp/beyondcell_pancancer.Rds")

## generate mats for limma
metacoms <- mcs@meta.data %>%
    select(metacom_untreated_1:metacom_treated_6) %>%
    rownames_to_column("cell_id")

functional_mat <- as.matrix(functional_mat)
