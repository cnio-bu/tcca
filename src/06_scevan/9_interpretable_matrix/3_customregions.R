#!/usr/bin/env Rscript

#sbatch -c 10 --job-name=custommap -o log_regmapping.txt -e error_regmapping.txt --mem=50G -t 10 --wrap "Rscript 3_customregions.R"

library(BPCells)
library(tidyverse)
library(Matrix)

## Functions

average_duplicate_rows <- function(initial_mat) {
  # Initialize the final matrix
  unique_rows <- unique(rownames(initial_mat))
  
  final_mat <- matrix(NA, nrow = length(unique_rows), ncol = ncol(initial_mat),
                      dimnames = list(unique_rows, colnames(initial_mat)))
  
  # Loop through unique row names and calculate averages
  for (row_name in unique_rows) {
    
    # Calculate row-wise averages excluding NA values
    rows_with_name <- initial_mat[rownames(initial_mat) == row_name, , drop = FALSE]
    
    col_average <- colMeans(rows_with_name, na.rm = TRUE)
    
    # Replace NaN with NA
    col_average[!is.finite(col_average)] <- NA
    
    # Assign the average values to the final matrix
    final_mat[row_name, ] <- col_average
  }
  
  return(final_mat)
}

## Load data

mat <- open_matrix_dir("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_fullbpcellsmatrix")
mat <- as.matrix(mat)

renaming <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/scevan_to_custom_segments.tsv", sep = "\t")

##This has to be performed in the cluster
## Add rownames to original matrix
new_mat <- mat[match(renaming$scevan_name, rownames(mat)), , drop = FALSE]
rownames(new_mat) <- renaming$custom_name

## Merge duplicated rownames
final_mat <- average_duplicate_rows(new_mat)
write.table(final_mat, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_segments_clones_lvl2_customregions.tsv", sep = "\t")