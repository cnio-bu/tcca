library(ComplexHeatmap)
library(uwot)
library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

## Load full beyondcell mat
mat <- open_matrix_dir(dir = "results/beyondcell_bp/full_mat_beyondcell/")
## Load full database
meta.data_full_clinical <- read_tsv(
    file = "results/annotation/beyondcell_metadata_with_clinical.tsv"
    ) %>%
    mutate(
        original_cell_id = gsub(
            pattern = "\\.\\.\\..*$", ## annoying ... by seurat
            replacement = "",
            x = cell
            ),
        new_cell_id = c(1:ncol(mat))
    )

meta.data_full_clinical <- meta.data_full_clinical %>%
    as.data.frame()

rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

## rename mat cells while preserving original in metadata
colnames(mat) <- c(1:ncol(mat))

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "beyondcell_pancancer",
    meta.data = meta.data_full_clinical
)

seu[["RNA"]]$data <- mat

## Use the full dataset as variable feat. We are talking 581. Let PCA summarise.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(mat)

## Go for sketch, it's faster
seu <- SketchData(
    object = seu,
    ncells = 50000,
    method = "LeverageScore",
    over.write = TRUE,
    sketched.assay = "sketch_50k",
    seed = 120394
)

DefaultAssay(seu) <- "sketch_50k"

seu <- FindVariableFeatures(
    object = seu,
    selection.method = "disp",
    nfeatures = 300
    )

seu <- ScaleData(seu, do.scale = TRUE, do.center = TRUE, scale.max = 20)
seu <- RunPCA(seu, npcs = 100, reduction.name = "pca")

seu <- FindNeighbors(
    object = seu,
    dims = 1:50,
    reduction = "pca"
    )

seu <- FindClusters(
    seu,
    resolution = 0.2,
    method = "igraph",
    algorithm = 2,
    random.seed = 120394,
    group.singletons = TRUE,
    cluster.name = "therapeutic_clusters_sketched_0.2"
    )

seu <- RunUMAP(seu, dims = 1:50, return.model = T, umap.method = "uwot")

seu <- Seurat::ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch_50k",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:50,
    refdata = list(therapeutic_clusters_0.2 = "therapeutic_clusters_sketched_0.2")
)

write_tsv(seu@meta.data, "results/annotation/tcs.tsv")

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
write_matrix_dir(mat = sketched_mat, dir = "results/beyondcell_bp/sketch_mat_beyondcell")
write_matrix_dir(mat = sketched_mat_5k, dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k")
