#!/usr/bin/env Rscript

library("tidyverse")

## Get input
args <- commandArgs(trailingOnly = TRUE)
filename_complete <- as.character(args[1])
parts <- unlist(strsplit(filename_complete, "/"))
filename <- tools::file_path_sans_ext(parts[length(parts)])

##FILTER CELL X GENE TABLES TO GO FROM LVL 2 TO LVL 3

## Get input
clonality_table <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl3.tsv", sep = "\t", header = T)
cellxgene_table <- read.table(filename_complete, sep = "\t", header = T)

## Add study__sample column to filter
clonality_table <- clonality_table %>%
  mutate(barcode__study__sample = paste(scevan_barcode, study, sample, sep = "__"))

cellxgene_table <- cellxgene_table %>%
  mutate(barcode__study__sample = paste(barcode, study, sample, sep = "__"))

#Filter table and save
cellxgene_table <- cellxgene_table %>%
  filter(barcode__study__sample %in% clonality_table$barcode__study__sample) #Samples in lvl3
cellxgene_table <- subset(cellxgene_table, select = -barcode__study__sample) 

write.table(cellxgene_table, paste0("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl3/", filename, ".tsv"), sep = "\t", row.names = FALSE)
