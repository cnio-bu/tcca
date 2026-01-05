#!/usr/bin/env Rscript

library("Seurat")
library("tidyverse")

## Get input
args <- commandArgs(trailingOnly = TRUE)
filename_complete <- as.character(args[1])
filename <- gsub("_v5", "", filename_complete)

## Set paths
annotated_seurat_list <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/", filename, ".rds")

## Set saving directories
where_to_save <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/clonality_metadata/", filename, ".tsv")

## LOAD ALL DATA
## Load all lvl2 data
seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
seu@meta.data <- mutate(seu@meta.data,
                        sample_study = paste0(sample, "__", study))

## Get study and sample names
original_seu <- subset(seu, study == filename)

## Read annotated data
annotated_seu <- readRDS(file = annotated_seurat_list)


## FUNCTIONS

#Get both cell names
get_cellnames <- function(sc){
  if (!is.null(sc)){
    newcells <- gsub("\\.", "", Cells(sc))
    newcells <- gsub("-", "", newcells)
    newcells <- gsub("_", "", newcells)
    
    cells <- Cells(sc)
    all <- data.frame(original_barcode = cells, scevan_barcode = newcells)
    return(all)
  }
  else{
    return(NULL)
  } 
}

#Get clonality information
get_clonality <- function(sc){
  if (!is.null(sc)){
    clonality <- select(sc@meta.data, sample, malignancy, scevan_prediction, confidentNormal, scevan_subclone)
    clonality <- tibble::rownames_to_column(clonality, "scevan_barcode")
    return(clonality)
  }
  else{
    return(NULL)
  } 
}

#Merge dataframes
merge_function <- function(df) {
  if (!is.null(df)){
    merged_df <- df %>% left_join(cellnames, by = "scevan_barcode")
    return(merged_df)
  }
  else{
    return(NULL)
  } 
}

## Get all cell names and clonality dataframes
cellnames <- get_cellnames(original_seu)
clonality_list <- lapply(annotated_seu, get_clonality)

## Use Map to merge corresponding elements from both lists
full_list <- lapply(clonality_list, merge_function)
full_table <- do.call(rbind, full_list)
full_table <- transform(full_table,
                        study = filename)
rownames(full_table) <- NULL

## Save
write.table(full_table, file = where_to_save, sep = "\t")