library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5")
seu <- readRDS("lvl1/seu_lvl1.rds")
seu <- JoinLayers(seu)

seu <- NormalizeData(seu)
seu$study_sample <- paste(seu$study, seu$sample, sep = "_")
seu_bulk <- AggregateExpression(
    object = seu,
    slot = "counts",
    return.seurat = T,
    group.by = c("study_sample")
    )

mat <- seu_bulk[["RNA"]]$counts
mat <- as.matrix(mat)

write.table(x = mat, file = "lvl1/pancancer_pseudobulk.tsv")