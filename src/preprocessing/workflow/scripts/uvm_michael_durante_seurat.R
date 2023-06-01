library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file      <- snakemake@input[["matrix"]]
metadata      <- snakemake@input[["metadata"]]
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

scale_data_find_variables <- function(sc) {
    sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst")
    sc <- Seurat::ScaleData(sc, features = rownames(sc))
    return(sc)
}

mat <- Seurat::Read10X_h5(filename = mat_file)
metadata <- data.table::fread(input = metadata)  %>%
    as.data.frame()

rownames(metadata) <- metadata$Cell
metadata$Cell <- NULL

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    meta.data = metadata,
    project = "uvm"
    )

seu_list <- Seurat::SplitObject(object = seu, split.by = "Sample")

seu_list <- lapply(seu_list, filter_sc)
seu_list <- lapply(seu_list, scale_data_find_variables)

saveRDS(seu_list, where_to_save)