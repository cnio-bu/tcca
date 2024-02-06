library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

## Load full beyondcell mat
mat <- open_matrix_dir(dir = "results/functional/full_mat_functional/")
## Load full database
meta.data_full_clinical <- read_tsv(
    file = "results/annotation/functional_metadata_with_clinical.tsv"
)

meta.data_full_clinical <- meta.data_full_clinical %>%
    as.data.frame()

rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id
mat <- mat[1:130, ]

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "functional_pancancer",
    meta.data = meta.data_full_clinical
)

seu[["RNA"]]$data <- mat

## Use the full dataset as variable feat. We are talking 130 Let PCA summarise.
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
    cluster.name = "functional_clusters_sketched_0.2"
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
    refdata = list(functional_clusters_0.2 = "functional_clusters_sketched_0.2")
)

write_tsv(seu@meta.data, "results/annotation/fcs.tsv")

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
write_matrix_dir(mat = sketched_mat, dir = "results/functional/sketch_mat_functional")
write_matrix_dir(mat = sketched_mat_5k, dir = "results/functional/sketch_mat_functional_5k")

## export object
saveRDS(object = seu, file = "results/functional/beyondcell_pancancer.Rds")
