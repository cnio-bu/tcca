#!/usr/bin/env Rscript

library("tidyverse")

## Get input
args <- commandArgs(trailingOnly = TRUE)
filename_complete <- as.character(args[1])
filename <- basename(filename_complete)

## Set paths
cna_mtxs <- list.files(path = paste0(filename_complete, "/output"), pattern = "CNAmtx.RData$", full.names = T)
cna_annotations <- list.files(path = paste0(filename_complete, "/output"), pattern = "count_mtx_annot.RData$", full.names = T)

## Set saving directories
where_to_save <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes/", filename, ".tsv")

## Read data
list_cna <- lapply(cna_mtxs, function(file) {
  loaded_data <- get(load(file))
  loaded_data <- t(loaded_data)
  name <- str_extract(file, "([^/]+)_CNAmtx.RData$") %>%
    str_remove("_CNAmtx.RData$")
  return(list(name = name, data = loaded_data))
})

names_list <- lapply(list_cna, function(x) x$name)
list_cna <- lapply(list_cna, function(x) x$data)
names(list_cna) <- names_list


list_annotation <- lapply(cna_annotations, function(file) {
  loaded_data <- get(load(file))
  loaded_data <- unite(loaded_data, gene_information, seqnames, start, end, gene_id, gene_name, sep = "_") #combine all information into a single column
  rownames(loaded_data) <- NULL
  name <- str_extract(file, "([^/]+)_count_mtx_annot.RData$") %>%
    str_remove("_count_mtx_annot.RData$")
  return(list(name = name, data = loaded_data))
})

names_list <- lapply(list_annotation, function(x) x$name)
list_annotation <- lapply(list_annotation, function(x) x$data)
names(list_annotation) <- names_list


# Get the common names (shall be all of them)
common_names <- intersect(names(list_cna), names(list_annotation))

# Initialize an empty list to store the combined data frames
annotated_list <- list()

# Loop through the common names and annotate each element of the list
for (name in common_names) {
  study <- filename
  sample <- name
  annotated_list[[name]] <- as.data.frame(list_cna[[name]])
  colnames(annotated_list[[name]]) <- list_annotation[[name]]$gene_information #Annotate gene names
  annotated_list[[name]] <- cbind(study, sample, annotated_list[[name]]) #Add sample and study
  annotated_list[[name]] <- rownames_to_column(annotated_list[[name]], var = "barcode")
}


##FUNCTION TO MERGE DATAFRAMES IN A LIST BY COLUMN NAME

merge_by_colname <- function(colnames, dataframes) {
  # Create a template data frame with NA for all columns
  template_df <- tibble::tibble(!!!setNames(rep(list(NA), length(colnames)), colnames))
  
  # Function to fill in missing columns with NA
  fill_na <- function(df) {
    missing_cols <- setdiff(colnames, colnames(df))
    df[, missing_cols] <- NA
    return(df)
  }
  
  # Fill in missing columns with NA for each dataframe
  filled_dataframes <- map(dataframes, fill_na)
  
  # Reduce (merge) the data frames into a single data frame
  result_df <- reduce(filled_dataframes, full_join, by = colnames)
  
  return(result_df)
}


# Get common gene names
common_genes <- NULL

for (i in annotated_list){
  common_genes <- c(common_genes, as.character(colnames(i)))
}

common_genes <- unique(common_genes)

#Get full table
full_table <- merge_by_colname(common_genes, annotated_list)

## Save
write.table(full_table, file = where_to_save, sep = "\t", row.names = FALSE)