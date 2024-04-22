library(BPCells)
library(Seurat)
library(tidyverse)

setwd("Documentos/bc-meta")

## load bc mat
bc <- readRDS("results/beyondcell_bp/beyondcell_pancancer.Rds")

paste0("BC object size:", print(dim(bc)))

## subset bc to remove bhupinder pal samples
bc_subset <- subset(bc, 
                    subset = (tumor_type == "BRCA" & tumor_subtype != "predicted_tumour" | tumor_type != "BRCA")
                    )

paste0("Subset bc object size", print(dim(bc_subset)))

# load func from old mat
functional_mat <- readRDS("results/functional/metaprograms_sc_mat.rds")
paste0("func mat object size:", print(dim(functional_mat)))

## load seu mat
seu <- readRDS("results/tcca/tcca_seurat_raw.rds")

seu <- subset(seu, subset = malignancy == TRUE)
paste0("seu  object size:", print(dim(seu)))

paste0("Cell diff ", as.character(ncol(seu) - ncol(bc_subset)))

## missing dalia?
samples_diff <- setdiff(unique(seu@meta.data$sample), unique(bc_subset@meta.data$sample))

print(nrow(seu@meta.data[seu@meta.data$sample %in% samples_diff, ]))

b <- seu@meta.data[seu@meta.data$sample %in% samples_diff, ]      

bc_set_data <- table(bc_subset@meta.data$sample)
seu_set_data <- table(seu@meta.data$sample)

seu_sample_cells <- as.data.frame(seu_set_data)
bc_sample_cells <- as.data.frame(bc_set_data)

## bc data names
colnames(bc_sample_cells) <- c("sample", "n.cells")
colnames(seu_sample_cells) <- c("sample", "n.cells")

## merge data
merged_set <- seu_sample_cells %>%
    full_join(
        y = bc_sample_cells,
        by = c("sample")
    )

colnames(merged_set) <- c("sample", "n.cells_seu", "n.cells_bc")
merged_set[is.na(merged_set$n.cells_bc), "n.cells_bc"] <- 0

merged_set$set_diff <- merged_set$n.cells_seu - merged_set$n.cells_bc