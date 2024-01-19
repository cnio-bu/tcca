#!/usr/bin/env Rscript

library("SCEVAN")
library("Seurat")
library("tidyverse")

## Get input
args <- commandArgs(trailingOnly = TRUE)
full_seurat_list <- as.character(args[1])

## Get file name from path
## /storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/file_name.rds
parts <- unlist(strsplit(full_seurat_list, "/"))
filename <- tools::file_path_sans_ext(parts[length(parts)])

## Set saving directories
where_to_save <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/", filename, ".rds")
cna_dir <- paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna/", filename)

## params
threads_to_use <- 19

## Function definition
SCEVAN_pred <- function(sc){
    
    sample_name <- unique(sc@meta.data$sample) #Get sample for the error
    
    tryCatch(
        {
            #Correction of cellnames
            newcells <- gsub("\\.", "-", Cells(sc)) #SCEVAN has problems with rownames containing both "." and "-" and it crashes with some samples ("undefined columns selected")
            sc <- RenameCells(sc, new.names = newcells)
            
            #Get vector of normal cells
            normalcells <- rownames(sc@meta.data[sc@meta.data$malignancy == FALSE,])
            
            #Generate SCEVAN annotation
            sample <- unique(sc@meta.data$sample) #Get sample for naming the output results
            count_mtx <- Seurat::GetAssayData(object = sc, slot = "counts")
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


# Create the dir for CNA
dir.create(cna_dir, showWarnings = FALSE)
setwd(cna_dir)

## Read data
seu <- readRDS(file = full_seurat_list)

## Open log files for appending
log_file <- file("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/annotated/log_error.txt", "a")

## Fill in SCEVAN predictions with standarized format: scevan_prediction, scevan_subclone
full_annotated_list_scevan <- lapply(seu, SCEVAN_pred)

## Close log files
close(log_file)

## Save
saveRDS(object = full_annotated_list_scevan, file = where_to_save)

