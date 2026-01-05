#!/usr/bin/env Rscript

library(SCEVAN)
library(BPCells)
library(Seurat)
library(tidyverse)
library(Matrix)

## INPUT
## Get study path
args <- commandArgs(trailingOnly = TRUE)
full_seurat_list <- as.character(args[1])

## Set threads_to_use from argument or default
if (length(args) >= 2) {
  threads_to_use <- as.integer(args[2])
} else {
  threads_to_use <- 19
}

## Get file name from path
## /storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/file_name.rds
parts <- unlist(strsplit(full_seurat_list, "/"))
filename <- gsub("_v5", "", parts[length(parts)])

## Set saving directories
where_to_save <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/", filename, ".rds")
cna_dir <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna/", filename)

### FUNCTION
SCEVAN_pred <- function(sc){
    
    sample_name <- unique(sc@meta.data$sample_study) #Get sample for the error
    
    tryCatch(
        {
            #Correction of cellnames
            newcells <- gsub("\\.", "", Cells(sc)) #SCEVAN has problems with rownames containing both "." and "-" and it crashes with some samples ("undefined columns selected")
            newcells <- gsub("-", "", newcells)
            newcells <- gsub("_", "", newcells)  

            sc <- RenameCells(sc, new.names = newcells)
            
            #Get vector of normal cells
            normalcells <- rownames(sc@meta.data[sc@meta.data$malignancy == FALSE,])
            
            #Generate SCEVAN annotation
            sample <- unique(sc@meta.data$sample) #Get sample for naming the output results
            count_mtx <- Seurat::GetAssayData(object = sc, slot = "counts")
            count_mtx <- as.matrix(count_mtx)
            count_mtx <- Matrix(count_mtx, sparse = TRUE)
            results <- SCEVAN::pipelineCNA(count_mtx = count_mtx,
                                 sample = sample,
                                 par_cores = threads_to_use,
                                 SUBCLONES = TRUE,
                                 ClonalCN = FALSE, #otherwise it crashes with some samples ("replacement has 2 rows, data has 1")
                                 plotTree = FALSE,
                                 SCEVANsignatures = TRUE,
                                 organism = "human",
                                 norm_cell = normalcells) 
                                 
            #Change NA values
            results["confidentNormal"][is.na(results["confidentNormal"])] <- "no"
            results["subclone"][is.na(results["subclone"])] <- "non_tumor"
            
            #Rename columns to standarized names: scevan_prediction, scevan_subclone
            colnames(results)[colnames(results) == "class"] <- "scevan_prediction"
            colnames(results)[colnames(results) == "subclone"] <- "scevan_subclone"
            
            #Fill in SCEVAN predictions
            sc_scevan <- Seurat::AddMetaData(sc, metadata = results)
            return(sc_scevan)

        },

        error = function(e){
            #Print error to the log file
            cat("The sample", as.character(sample_name), "from", as.character(filename), "is screwing you up with the error: \n", conditionMessage(e), "\n", file=log_file, append = TRUE)
            return(NULL)

        }
    )

}

## CODE
# Create the dir for CNA
dir.create(cna_dir, showWarnings = TRUE)
setwd(cna_dir)

## Load all lvl2 data
seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
seu@meta.data <- mutate(seu@meta.data,
                        sample_study = paste0(sample, "__", study))

## Get study and sample names
sc <- subset(seu, study == filename)
all_samples <- unique(sc@meta.data$sample_study)

count_mats <- list()

for (sam in all_samples){
  sample_seu <- subset(sc, sample_study == sam)
  count_mats <- c(count_mats, list(sample_seu))
}

names(count_mats) <- all_samples

## Open log files for appending
log_file <- file("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/log_error.txt", "a")

## Fill in SCEVAN predictions with standarized format: scevan_prediction, scevan_subclone
full_annotated_list_scevan <- lapply(count_mats, SCEVAN_pred)

## Close log files
close(log_file)

## Save
saveRDS(object = full_annotated_list_scevan, file = where_to_save)

