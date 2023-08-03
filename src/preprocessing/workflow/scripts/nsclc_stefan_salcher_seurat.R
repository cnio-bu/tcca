library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file      <- snakemake@input[["matrix"]] 
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

seu <- readRDS(mat_file) 
DefaultAssay(seu) <- "RNA"

seu_list <- Seurat::SplitObject(object = seu, split.by = "sample")
names(seu_list) <- unique(seu$sample)
samples_to_filter <- c("Leader_Merad_2021_336", "Leader_Merad_2021_298")

seu_list <- seu_list[!(names(seu_list) %in% samples_to_filter)]
seu_list <- lapply(seu_list, filter_sc)

## Get rid of NULLs samples w/o malignants left
seu_list[sapply(seu_list, is.null)] <- NULL
## In this datasets, some samples had very few cells, less than the minimum
## for calculating variable features
seu_list <- seu_list[sapply(seu_list, ncol) > 10]

seu_list <- lapply(seu_list, normalize_and_scale)

saveRDS(seu_list, where_to_save)