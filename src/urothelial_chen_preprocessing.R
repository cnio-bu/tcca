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
load(file = paste0(data_folder,"/raw/urothelial_chen/sc_Chen_BUC.rda"))

## Reset idents
sc_Chen_BUC$old_ident <- sc_Chen_BUC$orig.ident
sc_Chen_BUC$orig.ident <- sc_Chen_BUC$Sample

## New column TRUE-FALSE Carcinoma-Normal sample
sc_Chen_BUC@meta.data$Malignant <- ifelse(!grepl("NBM",sc_Chen_BUC$Sample),
                                          TRUE,FALSE)

## New column TRUE-FALSE Epithelial-Non epithelial cell
sc_Chen_BUC@meta.data$Epithelial <- ifelse(grepl(
  "Epithelial",sc_Chen_BUC$cluster_sample),TRUE,FALSE)

## Get epithelial and malignant cells
sc_Chen_BUC <- subset(sc_Chen_BUC,(
  subset = (Epithelial == TRUE & Malignant == TRUE)))

## Split the merged obj
samples_list <- Seurat::SplitObject(object = sc_Chen_BUC, 
                                    split.by = "Sample")

## Filter cells
samples_filtered <- lapply(
    samples_list,
    filter_sc,
    res_dir = paste0(data_folder, "/qc/urothelial_chen"))

## Normalize and scale data
malignant_cells_result <- lapply(samples_filtered, normalize_and_scale)

## Seurat object
saveRDS(
  object = malignant_cells_result,
  file = paste0(data_folder, "/obj/urothelial_chen/all_malignant.rds"))

## Number of malignant cells
sum(sapply(malignant_cells_result, ncol))
