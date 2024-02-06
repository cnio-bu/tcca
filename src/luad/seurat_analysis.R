library(BPCells)
library(Seurat)
library(tidyverse)

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

mat <- readRDS("raw/luad_kim_nayoung/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
write_matrix_dir(mat = as.sparse(mat), dir = "results/luad/luad_nayoung_bp")
