library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional/")

## Load Ucell score matrix
mat <- open_matrix_dir(dir = "full_mat_ucell")

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

# Remove cells with therapeutic cluster information from scTherapy
tcca_metadata <- read.table(
    "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t", header = TRUE
)
cells_tcs <- tcca_metadata %>%
    filter(malignancy == "True" & !is.na(scTherapy_cluster)) %>%
    pull(cell)

seu <- subset(seu, subset = cell %in% cells_tcs)
seu[["RNA"]]$data <- seu[["RNA"]]$counts

## Use the full dataset as variable feat.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu[["RNA"]]$counts)

## Go for sketch, it's faster
seu <- SketchData(
    object = seu,
    ncells = 5000,
    over.write = TRUE,
    sketched.assay = "sketch_5k",
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
    reduction = "pca",
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

write_tsv(seu@meta.data, "results/fcs.tsv")

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
write_matrix_dir(mat = sketched_mat_5k, dir = "results_ucell/sketch_mat_functional_5k_only_sctherapycells")

## export object
saveRDS(object = seu, file = "results_ucell/ucell_pancancer_only_sctherapy_cells.Rds")
