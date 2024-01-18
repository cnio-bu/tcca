#!/usr/bin/env Rscript

library("Seurat")
library("tidyverse")

## Get input
args <- commandArgs(trailingOnly = TRUE)
filename_complete <- as.character(args[1])
filename <- tools::file_path_sans_ext(filename_complete)

## Set paths
original_seurat_list <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/", filename_complete)
annotated_seurat_list <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/", filename_complete)

## Set saving directories
where_to_save <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/clonality_metadata/", filename, ".tsv")

## params
threads_to_use <- 19

## Read data
original_seu <- readRDS(file = original_seurat_list)
annotated_seu <- readRDS(file = annotated_seurat_list)

## FUNCTIONS

#Get both cell names
get_cellnames <- function(sc){
  if (!is.null(sc)){
    newcells <- gsub("\\.", "-", Cells(sc))
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
merge_function <- function(df_A, df_B) {
  if (!is.null(df_A) & !is.null(df_B)){
    merged_df <- merge(df_A, df_B, by = "scevan_barcode")
    return(merged_df)
  }
  else{
    return(NULL)
  } 
}

## Get all cell names and clonality dataframes
cellnames_list <- lapply(original_seu, get_cellnames)
clonality_list <- lapply(annotated_seu, get_clonality)

## Use Map to merge corresponding elements from both lists
full_list <- Map(merge_function, cellnames_list, clonality_list)
full_table <- do.call(rbind, full_list)
full_table <- transform(full_table,
                        study = filename)
rownames(full_table) <- NULL

## Save
write.table(full_table, file = where_to_save, sep = "\t")