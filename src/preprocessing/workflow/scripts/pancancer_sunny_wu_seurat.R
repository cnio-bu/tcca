library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
metadata <- snakemake@input[["metadata"]]
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


all_cells <- Seurat::Read10X(
data.dir = data_directory,
    gene.column = 1
)

meta_data <- read.delim(
    file = metadata
    ) %>%
    slice(2:n()) %>%
    as.data.frame()

rownames(meta_data) <- meta_data$NAME
meta_data$NAME <- NULL

all_cells <- Seurat::CreateSeuratObject(counts = all_cells,
                                        project = "pancancer_sunny_wu",
                                        meta.data = meta_data
                                        )


all_cells$old_ident <- all_cells$orig.ident
all_cells$orig.ident <- all_cells$biosample_id
    
sample_list <- Seurat::SplitObject(object = all_cells, split.by = "biosample_id")
names(sample_list) <- sapply(sample_list, function(sc){unique(sc$"biosample_id")})

# Filter cells
samples_filtered <- lapply(
    sample_list,
    filter_sc,
)

## Normalize and scale data
samples_filtered <- lapply(samples_filtered, normalize_and_scale)

## Get rid of the NULL elements/samples_filtered with no malignant cells left
samples_filtered[sapply(malignant_cells, is.null)] <- NULL

saveRDS(
    object = samples_filtered,
    file = where_to_save
    )
