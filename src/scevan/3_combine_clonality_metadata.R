## COMBINE MATRIXES OF CLONALITY METADATA
library(tidyverse)
library(Seurat)
library(beyondcell)

#Load all samples
all_samples <- list.files("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/clonality_metadata",
                         recursive = FALSE,
                         full.names = FALSE
)

full_paths <- list.files("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/clonality_metadata",
                          recursive = FALSE,
                          full.names = TRUE
)

all_mats <- lapply(X = full_paths, FUN = function(x){read.table(x, header = T)})
names(all_mats) <- all_samples

## Generate full matrix
full_table <- do.call(rbind, all_mats)
rownames(full_table) <- NULL

## Save lvl 2 full table
write.table(full_table, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl2.tsv", sep = "\t")

##Filter table to contain only malignant cells (from lvl 2 to lvl 1)
full_table_1 <- full_table %>%
  filter(malignancy == TRUE)

write.table(full_table_1, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl1.tsv", sep = "\t")