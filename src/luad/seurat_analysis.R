library(BPCells)
library(Seurat)
library(tidyverse)

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

mat <- readRDS("raw/luad_kim_nayoung/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
write_matrix_dir(mat = as.sparse(mat), dir = "results/luad/luad_nayoung_bp")

mat <- open_matrix_dir(dir = "results/luad/luad_nayoung_bp")
meta.data <- data.table::fread(
    "raw/luad_kim_nayoung/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    )

meta.data <- as.data.frame(meta.data)
rownames(meta.data) <- meta.data$Index

seu <- CreateSeuratObject(
    counts = mat,
    assay = "luad_pancancer",
    meta.data = meta.data
    )

seu <- NormalizeData(seu)

## test for matches
is_tumor <- c("mBrain", "tLung", "PE", "mLN", "tL/B")
meta.data_tumor <- meta.data[meta.data$Sample_Origin %in% is_tumor, ]
meta.data_tumor <- meta.data_tumor[meta.data_tumor$Cell_subtype %in% c("Malignant cells", "tS1", "tS2", "tS3"), ]
