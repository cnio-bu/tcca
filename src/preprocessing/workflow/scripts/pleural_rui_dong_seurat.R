library("Seurat")
library("tidyverse")
#hdf5r

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
metadata <- read.csv(metadata, sep = "\t")

rownames(metadata) <- metadata$Cell

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    project = "Pleuropulmonary blastoma",
    meta.data = metadata
    )

## This is dumb but this way 
## we'll standarize this study with the others with >1 sample
seu <- list(seu)
seu <- lapply(seu, filter_sc)
seu <- lapply(seu, scale_data_find_variables)

saveRDS(seu, where_to_save)