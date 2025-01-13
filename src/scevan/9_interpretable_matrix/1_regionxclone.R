#!/usr/bin/env Rscript

#sbatch -c 10 --job-name=segmCN -o log.txt -e error.txt --mem=40G -t 300 --wrap "Rscript 1_regionxclone.R"

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
library(stringr)

#Function cbind.fill
cbind.fill<-function(mat.list, genes){
    fill_na <- function(m) {
        m <- as.matrix(m)
        missing_rows <- setdiff(genes, rownames(m))
        if (length(missing_rows) > 0) {
            extra.mat <- matrix(NA, length(missing_rows), ncol(m))
            rownames(extra.mat) <- missing_rows
            m <- rbind(m, extra.mat)
        }
        
        name <- colnames(m)
        m <- as.matrix(m[order(rownames(m)),])
        colnames(m) <- name
        m <- Matrix(m, sparse = T)
        
        write_matrix_dir(
            mat = m,
            dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_bpcells_merged/", name),
            overwrite = TRUE)
        
        m <- open_matrix_dir(dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_bpcells_merged/", name))

        return(m)
    }
    
    filled_mats <- map(mat.list, fill_na)
    result_m <- reduce(filled_mats, cbind)
    return(result_m)
}

extract_clone_name <- function(path) {
  study <- str_extract(path, "(?<=single_cell/cna/)[^/]+")
  sample.clone <- str_extract(path, "([^/]+)_([^/]+)_CN\\.seg")
  sample <- str_extract(sample.clone, "(.*)(?=(_subclone))")
  clone <- str_extract(sample.clone, paste0("(?<=", sample, "_)(.*)")) %>%
    str_remove("_CN.seg")  
  name <- paste(study, sample, clone, sep = "__")
  return(name)
}

## Set paths and get clones names
cna_mtxs <- list.files(path = "/storage/scratch01/shared/projects/bc-meta/single_cell/cna", pattern = "_CN.seg$", full.names = T, recursive = T)
names <- lapply(cna_mtxs, extract_clone_name)

for (i in 1:length(cna_mtxs)){
  sample <- str_remove(names[[i]], "__subclone*.")
  
  loaded_data <- read.table(cna_mtxs[i])
  loaded_data <- mutate(loaded_data,
                        region = paste(Chr, Pos, End, sep = "_"))
  loaded_data <- loaded_data[,c("region", "CN")]
  colnames(loaded_data)[2] <- names[[i]]
  rownames(loaded_data) <- NULL
  loaded_data <- column_to_rownames(loaded_data, var = "region")
  loaded_data <- as.matrix(loaded_data)
  loaded_data <- Matrix(loaded_data, sparse = T)
  write_matrix_dir(
    mat = loaded_data,
    dir = paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_bpcells/", names[[i]]),
    overwrite = TRUE)
}

## Open all bpcells mats
file.dir <- "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_bpcells/"
files.set <- list.dirs(file.dir, full.names = FALSE, recursive = FALSE)

# Loop through directories and output BPCells matrices on-disk
data.list <- c()

for (i in 1:length(files.set)) {
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(file.dir, files.set[i]))
  # Add matrices to list (transposed to be genes x cells)
  data.list[[i]] <- mat
}

names(data.list) <- files.set

# Get all regions
all_regions <- NULL

for (i in data.list){
  all_regions <- c(all_regions, as.character(rownames(i)))
}

all_regions <- unique(all_regions)

#Compute cbind.filled matrix
full_mat <- cbind.fill(data.list, all_regions)

# Write the matrix to a directory
write_matrix_dir(
  mat = full_mat,
  dir = "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_fullbpcellsmatrix",
  overwrite = TRUE)

#Save matrix as Seurat RDS
full_mat <- open_matrix_dir("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_fullbpcellsmatrix/")
seu <- CreateSeuratObject(counts = full_mat)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/")

saveRDS(
  object = seu,
  file = "full_segments_copynumber.rds",
  destdir = "full_segments_copynumber_1layer",
  relative = TRUE
)