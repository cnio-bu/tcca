#!/usr/bin/env Rscript

#sbatch -c 20 --job-name=Sv5 -o log.txt -e error.txt --mem=60G -t 80 --wrap "Rscript 7_final_cellxgene_v5_cnafill.R"

library(BPCells)
library(Matrix)
library(tidyverse)
library(BPCells)
library(tidyr)
library(dplyr)
library(tibble)
library(Seurat)
library(SeuratObject)

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
        
        name <- unlist(strsplit(colnames(m)[1], "__"))[2]

        write_matrix_dir(
            mat = m,
            dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl2_bpcells_merged/", name),
            overwrite = TRUE)
        
        m <- open_matrix_dir(dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl2_bpcells_merged/", name))

        return(m)
    }
    
    filled_mats <- map(mat.list, fill_0)
    result_m <- reduce(filled_mats, cbind)
    return(result_m)
}

## Open all bpcells mats
file.dir <- "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl2_bpcells/"
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


# Get all gene names
all_genes <- NULL

for (i in data.list){
  all_genes <- c(all_genes, as.character(rownames(i)))
}

all_genes <- unique(all_genes)

#Compute cbind.filled matrix
full_mat <- cbind.fill(data.list, all_genes)

# Write the matrix to a directory
write_matrix_dir(
  mat = full_mat,
  dir = "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl2_fullbpcellsmatrix",
  overwrite = TRUE)

#Generate base metadata table
clonality_table <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl2.tsv", sep = "\t", header = T)
clonality_table <- clonality_table %>%
  mutate(barcode_study_sample = paste(scevan_barcode, study, sample, sep  = "__"),
         study_sample = paste(study, sample, sep  = "__"))

clinical_metadata <- read.table("/home/lserranor/clinical_metadata_v4_clean.tsv", sep = "\t", header = T)
clinical_metadata <- clinical_metadata %>%
  mutate(study_sample = paste(study, sample, sep  = "__"))

full_metadata_table <- merge(clonality_table, clinical_metadata, by = "study_sample")
full_metadata_table <- full_metadata_table[!duplicated(full_metadata_table$barcode_study_sample), ]
rownames(full_metadata_table) <- NULL

full_metadata_table <- full_metadata_table %>% 
  column_to_rownames(var = "barcode_study_sample") %>%
  subset(select = -c(study_sample, sample.y, study.y)) %>%
  rename(
    sample = sample.x,
    study = study.x
    )

#Save matrix as Seurat RDS
seu <- CreateSeuratObject(counts = full_mat, meta.data = full_metadata_table)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/")

saveRDS(
  object = seu,
  file = "full_genes_copynumber_lvl2.rds",
#  destdir = "full_genes_copynumber_1layer",
#  relative = TRUE
)
