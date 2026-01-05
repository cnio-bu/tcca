#!/usr/bin/env Rscript

library(BPCells)
library(tidyr)
library(tibble)
library(Matrix)
library(gtools)

## Get input
args <- commandArgs(trailingOnly = TRUE)
filename_complete <- as.character(args[1])
parts <- unlist(strsplit(filename_complete, "/"))
filename <- tools::file_path_sans_ext(parts[length(parts)])

#Read table
tsv <- read.table(filename_complete,
                  header = TRUE, sep = "\t")

#Combine metadata columns, replace NA and convert table to matrix
tsv <- tsv %>%
  unite(col = "barcode_study_sample", c("barcode", "study", "sample"), sep  = "__") %>%
  column_to_rownames(var = "barcode_study_sample") 

tsv <- tsv %>% na.replace(0)

tsv_mat <- as.matrix(tsv)
tsv_mat <- Matrix(tsv_mat, sparse = T)


# Write the matrix to a directory
write_matrix_dir(
  mat = tsv_mat,
  dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl2_bpcells/", filename),
  overwrite = TRUE)
