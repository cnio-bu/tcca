library("SCEVAN")
library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
full_seurat_list <- snakemake@input[["seurat_list"]]
where_to_save <- snakemake@output[["annotated_list"]]
cna_dir <- snakemake@params[["cna_res"]]

## SNAKEMAKE params
threads_to_use <- snakemake@threads

## Function definitions
SCEVAN_pred <- function(sc){
  #In this case, the SCEVAN prediction is performed previously, in
  #skcm_chao_zhang_seurat, as there is no malignancy annotation from the authors. 
  
  #Change NA values
  sc@meta.data["confidentNormal"][is.na(sc@meta.data["confidentNormal"])] <- "no"
  sc@meta.data["subclone"][is.na(sc@meta.data["subclone"])] <- "non_tumor"
  
  #Rename columns to standarized names: scevan_prediction, scevan_subclone
  colnames(sc@meta.data)[colnames(sc@meta.data) == "class"] <- "scevan_prediction"
  colnames(sc@meta.data)[colnames(sc@meta.data) == "subclone"] <- "scevan_subclone"

  return(sc)
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