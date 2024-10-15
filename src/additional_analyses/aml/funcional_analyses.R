library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)

seu <- readRDS("results/aml/malignants_seu.rds")

## Manually manipulate pointer
## I SET UP A soft link. TODO: FIX CODE ON WS
    
## Project back clusters
seu <- ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch_500",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:50,
    refdata = list(cluster_full = "integrated.clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(seu) <- "RNA"

seu_clusters_full <- DimPlot(
    object = seu,
    reduction = "umap.full",
    group.by = "Patient_ID"
    )

malignant_markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    logfc.threshold = 1,
    only.pos = TRUE
    )
