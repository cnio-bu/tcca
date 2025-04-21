#!/usr/bin/env Rscript

#sbatch -c 8 --job-name=clonexgene -o log.txt -e error.txt --mem=100G -t 60 --wrap "Rscript 7_final_clonexgene_v5_cnafill.R"

library(BPCells)
library(Matrix)
library(tidyverse)
library(BPCells)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(SeuratObject)
library(matrixStats)

#Function cbind.fill
cbind.fill<-function(mat.list, genes){
    
    fill_0 <- function(m) {
        m <- as.matrix(m)
        missing_rows <- setdiff(genes, rownames(m))
        if (length(missing_rows) > 0) {
            extra.mat <- matrix(0, length(missing_rows), ncol(m))
            rownames(extra.mat) <- missing_rows
            m <- rbind(m, extra.mat)
        }
        m <- m[order(rownames(m)),]
        m <- Matrix(m, sparse = T)
        
        return(m)
    }
    
    filled_mats <- map(mat.list, fill_0)
    result_m <- reduce(filled_mats, cbind)
    return(result_m)
}

## Open all bpcells mats
file.dir <- "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl3_bpcells/"
files.set <- list.dirs(file.dir, full.names = FALSE, recursive = FALSE)

# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()

for (i in 1:length(files.set)) {
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(file.dir, files.set[i]))
  # Add matrices to list (transposed to be genes x cells)
  data.list[[i]] <- t(mat)
}

names(data.list) <- files.set

#Generate base metadata table
clonality_table <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl3.tsv", sep = "\t", header = T)
clonality_table <- clonality_table %>%
  mutate(barcode_study_sample = paste(scevan_barcode, study, sample, sep  = "__"),
         study_sample = paste(study, sample, sep  = "__"))

clinical_metadata <- read.table("/home/lserranor/clinical_metadata_v4_clean.tsv", sep = "\t", header = T)
clinical_metadata <- clinical_metadata %>%
  mutate(study_sample = paste(study, sample, sep  = "__"))

full_metadata_table <- merge(clonality_table, clinical_metadata, by = "study_sample")
full_metadata_table <- full_metadata_table %>% 
  column_to_rownames(var = "barcode_study_sample") %>%
  subset(select = -c(study_sample, sample.y, study.y)) %>%
  rename(
    sample = sample.x,
    study = study.x
    ) %>%
  mutate(subclone_name = paste0(study, "__", sample, "__subclone", scevan_subclone)) ## Add subclone name for grouping

# Based on the annotation, generate a list of clones x genes matrixes
subclone_means.list <- list()

for (i in seq_along(data.list)) {
  # Select grouping variable as subclones
  bpcells <- t(data.list[[i]])
  metadata_filtered <- full_metadata_table[rownames(bpcells), , drop = FALSE]
  subclones <- metadata_filtered$subclone_name
  
  # Get subclone average per gene
  bpcells_by_subclone <- rowsum(as.matrix(bpcells), group = subclones)
  cell_counts_by_subclone <- table(subclones)
  subclone_gene_means <- sweep(bpcells_by_subclone, 1, cell_counts_by_subclone, FUN = "/")
  subclone_gene_means <- t(subclone_gene_means)
  
  # Save
  subclone_means.list[[i]] <- subclone_gene_means
}

names(subclone_means.list) <- files.set

# Get all gene names
all_genes <- NULL

for (i in subclone_means.list){
  all_genes <- c(all_genes, as.character(rownames(i)))
}

all_genes <- unique(all_genes)

# Compute cbind.filled matrix
full_mat <- cbind.fill(subclone_means.list, all_genes)
full_mat <- as.matrix(full_mat)

# Save
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/")

saveRDS(object = full_mat, file = "full_clone_gene_copynumber.rds")