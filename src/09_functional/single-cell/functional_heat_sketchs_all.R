library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional/")

## Load Ucell score matrix
mat <- open_matrix_dir(dir = "./full_mat_ucell")

## Load full Seurat object with Ucell scores in metadata
seu <- readRDS("seurat_ucell.rds")
meta.data <- seu@meta.data[, setdiff(colnames(seu@meta.data), rownames(mat))]

colnames(mat) <- paste0("c", c(1:ncol(mat)))
meta.data$new_cell_id <-  paste0("c", c(1:nrow(meta.data)))
rownames(meta.data) <- meta.data$new_cell_id


seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "functional_pancancer",
    meta.data = meta.data
)

seu[["RNA"]]$data <- seu[["RNA"]]$counts

## Use the full dataset as variable feat.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu[["RNA"]]$counts)

## Go for sketch, it's faster
seu <- SketchData(
    object = seu,
    ncells = 50000,
    over.write = TRUE,
    sketched.assay = "sketch_50k",
    seed = 120394
)

seu <- Seurat::SketchData(
    object = seu,
    assay = "RNA",
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch_5k"
)

## export  sketches
sketched_mat <- seu[["sketch_50k"]]$data
sketched_mat_5k <- seu[["sketch_5k"]]$data
write_matrix_dir(mat = sketched_mat, dir = "results_ucell/sketch_mat_functional")
write_matrix_dir(mat = sketched_mat_5k, dir = "results_ucell/sketch_mat_functional_5k")

## export object
saveRDS(object = seu, file = "results_ucell/functional_pancancer.Rds")
