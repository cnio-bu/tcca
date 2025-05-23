#!/usr/bin/env Rscript

#sbatch -c 20 --job-name=sketch -o log_sketch.txt -e error_sketch.txt --mem=200G -t 1140 --wrap "Rscript 8_sketch_rds.R"
#sbatch -c 20 --job-name=sketch -o log_sketch.txt -e error_sketch.txt --mem=200G -t 1140 --dependency=afterok:4381416 --wrap "Rscript 8_sketch_rds.R"

library(BPCells)
library(Matrix)
library(tidyverse)
library(BPCells)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata")
seu <- readRDS("full_genes_copynumber_1layer/full_genes_copynumber.rds")

seu <- NormalizeData(seu)
seu[["RNA"]]$data <- seu[["RNA"]]$counts
seu <- FindVariableFeatures(seu)
seu <- SketchData(
  object = seu,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

#sketch_mat <- as.matrix(seu[["sketch"]]$counts)
#sketch_mat <- Matrix(sketch_mat, sparse = T)

#write_matrix_dir(
#  mat = sketch_mat,
#  dir = "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_cnv",
#  overwrite = TRUE)

sketch_metadata <- seu@meta.data
sketch_metadata <- sketch_metadata[rownames(sketch_metadata) %in% colnames(seu[["sketch"]]$counts),]

write.table(sketch_metadata, 
  "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_metadata.tsv",
  row.names = T, 
  sep = "\t")


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
  "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/sketch_mat_metadata_5k.tsv",
  row.names = T, 
  sep = "\t")



setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/")

saveRDS(
  object = seu,
  file = "full_genes_copynumber.rds",
  destdir = "full_genes_copynumber_sketch",
  relative = TRUE
)







