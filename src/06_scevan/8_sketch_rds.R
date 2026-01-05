#!/usr/bin/env Rscript

#sbatch -c 20 --job-name=sketch -o log_sketch.txt -e error_sketch.txt --mem=30G -t 200 --wrap "Rscript 8_sketch_rds.R"
#sbatch -c 20 --job-name=sketch -o log_sketch.txt -e error_sketch.txt --mem=30G -t 200 --dependency=afterok:4658433 --wrap "Rscript 8_sketch_rds.R"

library(BPCells)
library(Matrix)
library(tidyverse)
library(BPCells)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(SeuratObject)
#library(SeuratDisk)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata")
seu <- readRDS("full_genes_copynumber_lvl2.rds")

seu <- NormalizeData(seu)
seu[["RNA"]]$data <- seu[["RNA"]]$counts
seu <- FindVariableFeatures(seu, assay = "RNA", nfeatures = 2000, layer="data")

#seu <- SketchData(
#  object = seu,
#  ncells = 50000,
#  method = "LeverageScore",
#  sketched.assay = "sketch"
#)

#sketch_mat <- as.matrix(seu[["sketch"]]$counts)
#sketch_mat <- Matrix(sketch_mat, sparse = T)

#write_matrix_dir(
#  mat = sketch_mat,
#  dir = "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_cnv",
#  overwrite = TRUE)

#sketch_metadata <- seu@meta.data
#sketch_metadata <- sketch_metadata[rownames(sketch_metadata) %in% colnames(seu[["sketch"]]$counts),]

#write.table(sketch_metadata, 
#  "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_metadata.tsv",
#  row.names = T, 
#  sep = "\t")


DefaultAssay(seu) <- "RNA"

seu <- SketchData(
  object = seu,
  ncells = 5000, #Smaller sketch for plots
  method = "LeverageScore",
  sketched.assay = "sketch_5k"
)

sketch_metadata <- seu@meta.data
sketch_metadata <- sketch_metadata[rownames(sketch_metadata) %in% colnames(seu[["sketch_5k"]]$counts),]

write.table(sketch_metadata, 
  "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_metadata_5k_lvl2.tsv",
  row.names = T, 
  sep = "\t")



setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/")

saveRDS(
  object = seu,
  file = "full_genes_copynumber_lvl2.rds"
#  destdir = "full_genes_copynumber_sketch",
#  relative = TRUE
)

## Save sketched matrix already scaled and normalized for plotting
sketched_mat <- as.matrix(seu[["sketch_5k"]]$data)

saveRDS(sketched_mat, file = "full_genes_copynumber_sketch_5k_scaled_normcounts_mtx_lvl2.rds")
