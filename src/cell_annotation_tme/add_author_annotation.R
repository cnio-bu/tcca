# This script aims to add the cell annotation from the authors that was
# not properly added with the preprocessing pipeline, which includes some NAs
# to some cells due to a mismatch in rownames.
library(Seurat)
library(BPCells)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/")

## Add cell annotation to study "adrenalnb_rui_chong" ##
ad <- readRDS("./seurat/all_cell_types/adrenalnb_rui_chong_v5.Rds")

tumor_data <- readr::read_csv("./raw/adrenalnb_rui_chong/tumor_dataset_annotation.csv")

tumor_data <- tumor_data %>%
  select(cellname, celltype) %>%
  as.data.frame()

rownames(tumor_data) <- tumor_data$cellname
tumor_data$cellname <- NULL
rownames(tumor_data) <- gsub("\\.1$|-1$", "", rownames(tumor_data))
selected_cells <- gsub("\\.1$|-1$", "", rownames(ad@meta.data))
ad$cell_type <- tumor_data[selected_cells, "celltype"]
saveRDS(ad, "./seurat/all_cell_types/adrenalnb_rui_chong_v5.Rds")


## Add cell annotation to study "aml_audrey_lasry" ##
aml <- readRDS("./seurat/all_cell_types/aml_audrey_lasry_v5.Rds")
raw_metadata <- readr::read_csv("./raw/aml_audrey_lasry/metadata_clustering_w_header_upd.csv") %>%
  as.data.frame() %>% column_to_rownames("NAME")

aml$cell_type <- raw_metadata[colnames(aml), "Cell_type_identity"]

saveRDS(aml, "./seurat/all_cell_types/aml_audrey_lasry_v5.Rds")

## Add cell annotation to study "brmets_hugo_gonzalez" (only 1650 cells) ##
brmets <- readRDS("./seurat/all_cell_types/brmets_hugo_gonzalez_v5.Rds")
annot_files <- list.files("./raw/brmets_hugo_gonzalez",
                          pattern = "_Cell_Types_Annotations.csv",
                          full.names = TRUE)

cell_annot.list <- list()
for (sample in unique(brmets$sample)) {
  file <- grep("GSM5645888_Melan_1", annot_files, value = TRUE)
  cell_annot <- read.csv(file, row.names = 1, stringsAsFactors = FALSE)
  selected_cells <- rownames(brmets@meta.data)[brmets$sample == sample]
  selected_cells <- gsub("_[0-9]*", "", selected_cells)
  rownames_annot <- gsub("^.*_", "", rownames(cell_annot))
  matching_rows <- rownames(cell_annot)[rownames_annot %in% selected_cells]
  cell_annot <- cell_annot[grep(selected_cells, rownames(cell_annot))]
  
  # Check all cells in the preprocessed object have annotation in the raw annotation files.
  annot_n <- unlist(lapply(cell_annot.list, nrow))
  brmets_n <- as.vector(table(brmets$sample))
  annot_n - brmets_n # The numbers do not match (redo annotations --> Oscar)
}