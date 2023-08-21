library("SCEVAN")
library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
full_seurat_list <- snakemake@input[["seurat_list"]]
where_to_save <- snakemake@output[["annotated_list"]]
cna_dir <- snakemake@params[["cna_res"]]

## SNAKEMAKE params
threads_to_use <- snakemake@threads

SCEVAN_pred <- function(sc){
  #Correction of cellnames
  newcells <- gsub("\\.", "-", Cells(sc)) #SCEVAN has problems with rownames containing both "." and "-" and it crashes with some samples ("undefined columns selected")
  sc <- RenameCells(sc, new.names = newcells)
  
  #Get vector of normal cells
  normalcells <- rownames(sc@meta.data[sc@meta.data$malignancy == FALSE,])
  
  #Generate SCEVAN annotation
  sample <- levels(sc@meta.data$orig.ident) #Get sample for naming the output results
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
}

# Create the dir for CNA
dir.create(cna_dir, showWarnings = FALSE)
setwd(cna_dir)

## Read data
seu <- readRDS(file = full_seurat_list)

## Fill in SCEVAN predictions with standarized format: scevan_prediction, scevan_subclone
full_annotated_list_scevan <- lapply(full_seurat_list, SCEVAN_pred)

## Save
saveRDS(object = full_annotated_list_scevan, file = where_to_save)