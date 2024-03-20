library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

seu <- readRDS("results/lvl1/seu_lvl1_merged.Rds")
mat <- open_matrix_dir(dir = "results/lvl1/pancancer_merged_mat")

seu@assays$RNA$counts <- mat

## Remove preneoplastic brca
seu <- subset(seu, subset = study != "brca_bhupinder_pal" | tumor_subtype != "predicted_tumour")

seu <- NormalizeData(seu)

seu_bulk <- AggregateExpression(
    object = seu,
    slot = "counts",
    return.seurat = T,
    group.by = c("sample", "study")
    )

mat <- seu_bulk@assays$RNA@counts
mat <- as.matrix(mat)

write.table(x = mat, file = "results/functional/pancancer_pseudobulk.tsv")
