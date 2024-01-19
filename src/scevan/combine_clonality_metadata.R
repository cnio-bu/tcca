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

##FILTER CLONALITY TABLE TO CONTAIN ONLY MALIGNANT CELLS AND REMOVE DUPLICATED SAMPLES
#(FROM LVL 1 TO LVL 3)

## Get input
clinical_metadata <- read.table("/home/lserranor/clinical_metadata_v4_clean.tsv", sep = "\t", header = T)

#Correct clonality table for M-P1 sample and save (LVL 1)
full_table$sample[full_table$sample == " M-P1"] <- "M-P1"
write.table(full_table, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl1.tsv", sep = "\t")

#Select only lvl 2 samples (including at least 100 malignant cells, excluiding duplicated samples) (LVL 2)
full_table <- full_table %>%
  mutate(study__sample = paste0(study, "__", sample))

clinical_metadata <- clinical_metadata %>%
  mutate(study__sample = paste0(study, "__", sample))

full_table_2 <- full_table %>%
  filter(study__sample %in% clinical_metadata$study__sample)

full_table_2 <- subset(full_table_2, select = -study__sample)

write.table(full_table_2, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl2.tsv", sep = "\t")

#Select only malignant cells (LVL 3)
full_table_3 <- full_table_2 %>%
  filter(malignancy == TRUE)

write.table(full_table_3, "/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl3.tsv", sep = "\t")