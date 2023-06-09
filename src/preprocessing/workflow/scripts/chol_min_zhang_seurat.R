library("Seurat")
library("tidyverse")
library("stringr")

## SNAKEMAKE I/O
metadata      <- snakemake@input[["metadata"]]
where_to_save <- snakemake@output[["seurat_list"]]
data_directory <- snakemake@params[["data_dir"]]

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

annotate_clinical_data <- function(sc){

    common_cells <- intersect(colnames(sc), metadata$Cell)
    this_meta <- metadata[metadata$Cell %in% common_cells, ]
    rownames(this_meta) <- this_meta$Cell
    sc <- sc[, common_cells]
    sc <- AddMetaData(sc, metadata = this_meta)
    return(sc)

}


all_samples <- list.files(data_directory, full.names = TRUE)

seu_list <- list()

for(sample in all_samples){

    sample_name <- basename(sample)
    sample_name <- gsub("GSM[0-9]+_", "", gsub("_UMI.csv", "", sample_name))
    mat <- read.csv(sample, sep = ",", header = T, row.names = 1)
    seu <- Seurat::CreateSeuratObject(
        counts = mat,
        assay = "RNA",
        project = sample_name
    )

    seu$orig.ident <- NULL
    seu <- AddMetaData(seu, sample_name, col.name = "orig.ident")

    seu_list <- c(seu_list, seu)
}

names(seu_list) <- sapply(basename(all_samples), function(x){
    substring(x, 1, (str_length(x)-4))
})

metadata <- data.table::fread(input = metadata)  %>%
    as.data.frame()
metadata$Cell <- gsub("^","X",metadata$Cell)
rownames(metadata) <- metadata$Cell


seu_list <- lapply(seu_list, filter_sc)
seu_list <- lapply(seu_list, normalize_and_scale)
seu_list <- lapply(seu_list, annotate_clinical_data)

saveRDS(seu_list, where_to_save)
