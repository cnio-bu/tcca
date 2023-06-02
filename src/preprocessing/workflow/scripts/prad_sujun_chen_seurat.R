library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
raw_matrix <- snakemake@input[["raw_matrix"]]
cell_annot <- snakemake@input[["cell_annot"]]
where_to_save <- snakemake@output[["seurat_list"]]

## Function definitions
filter_sc <- function(sc) {
    sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
    sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "percent.ribo")
    
    sc_filtered <- subset(x = sc, subset = (percent.mt <= 10) &
                              (nFeature_RNA >= 1000 & nFeature_RNA <= 7000) &
                              (nCount_RNA > 500) & (percent.ribo <= 40)
    )
    
    this_counts <- GetAssayData(sc_filtered, slot = "counts")
    nonzero_genes <- this_counts > 0
    
    # Keep genes whose expression is found in at least 5% of the sample
    sample_cell_cutoff <- round(ncol(sc_filtered) / 100 * 5, digits = 0)
    genes_to_keep <- Matrix::rowSums(nonzero_genes) >= sample_cell_cutoff
    
    new_filtered_sc <- CreateSeuratObject(
        counts = this_counts[genes_to_keep, ],
        meta.data = sc_filtered@meta.data
    )

    return(new_filtered_sc)
}

normalize_and_scale <- function(sc) {
    sc <- Seurat::NormalizeData(sc,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000
                                )
    sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst")
    sc <- Seurat::ScaleData(sc, features = rownames(sc))
    return(sc)
}

## Read data
mat <- data.table::fread(input = raw_matrix) %>%
  as.data.frame()

rownames(mat) <- mat$V1
mat$V1 <- NULL
mat <- as.matrix(mat)

## Read metadata
metadata <- data.table::fread(cell_annot) %>%
    as.data.frame()

rownames(metadata) <- metadata$Cell
metadata$Cell <- NULL 

## Merge data + metadata
seu <- Seurat::CreateSeuratObject(counts = mat,
                                  project = "prad_sujun_chen",
                                  meta.data = metadata
                                  )

## Split the merged obj
sample_list <- Seurat::SplitObject(object = seu, split.by = "Sample")

names(sample_list) <- sapply(sample_list, function(sc){unique(sc$Sample)})

## Filter cells
filtered_sc <- lapply(sample_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )
