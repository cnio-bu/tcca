library("tidyverse")
library("Seurat")

## Load common funs
source("src/sc_functions.R")
cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys
data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")

## Load data
load(file = paste0(data_folder,"/raw/adrenalnb_rui_chong/human_NB_subset_tumor.rda"))

## Reset idents
human_NB_subset_tumor$old_ident <- human_NB_subset_tumor$orig.ident
human_NB_subset_tumor$orig.ident <- human_NB_subset_tumor$Sample

## Split the merged obj
samples_list <- Seurat::SplitObject(object = human_NB_subset_tumor, 
                                    split.by = "Sample")

## Filter cells
samples_filtered <- lapply(
    samples_list,
    filter_sc,
    res_dir = paste0(data_folder, "/qc/adrenalnb_rui_chong"))

## Normalize and scale data
malignant_cells_result <- lapply(samples_filtered, normalize_and_scale)

## Seurat object
saveRDS(
  object = malignant_cells_result,
  file = paste0(data_folder, "/obj/adrenalnb_rui_chong/all_malignant.rds"))