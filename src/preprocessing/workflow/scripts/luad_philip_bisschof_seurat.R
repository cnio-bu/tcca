library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
metadata <- snakemake@input[["metadata"]]
infercnv_scores <- snakemake@input[["infercnv_scores"]]
data_directory <- snakemake@params[["data_dir"]]

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

## Function: Annotate with clinical data
fill_clinical <- function(sc) {
    
    this_meta <- patient_data %>%
        filter(sample_id == unique(sc$orig.ident)) %>%
        full_join(y = data.frame(cell_id = colnames(sc)), by = character()) %>%
        as.data.frame()
    
    rownames(this_meta) <- this_meta$cell_id
    this_meta$cell_id <- NULL
    
    sc <- AddMetaData(object = sc, metadata = this_meta)
    return(sc)
    
}

## Get all samples
all_samples <- list.dirs(data_directory,
                         recursive = FALSE,
                         full.names = FALSE
                         )

full_paths <- paste0(data_directory,
                     "/",
                     all_samples
                     )

all_mats <- lapply(X = full_paths, FUN = Seurat::Read10X)

patient_data <- readxl::read_excel(
    path = metadata
    ) %>%
    filter(tissue_type == "Tumor")

infer_cnv_scores <- read.csv(
    file = infercnv_scores,
    sep = "\t"
    ) %>%
    filter(tissue_type == "Tumor", !is.na(cell_id), cna_clone == "CNA") %>%
    mutate(
        cell_id = gsub(pattern = "^p.*_", replacement = "", x = cell_id)
    )

names(all_mats) <- all_samples

generate_seurat_objects <- function(dgMat, sample) {
       sample_name <- sample
       seu <- CreateSeuratObject(counts = dgMat,
                                 project = sample,
                                 assay = "RNA",
                                 )
       
       return(seu)
}

all_seurat_objects <- lapply(seq_along(all_mats), function(id) {
    generate_seurat_objects(all_mats[[id]], sample = all_samples[id])
})

# QC
filtered_sc <- lapply(all_seurat_objects, filter_sc)
# Normalize
filtered_sc <- lapply(filtered_sc, normalize_and_scale)
# Annotation
filtered_sc <- lapply(filtered_sc, fill_clinical)

saveRDS(object = filtered_sc, file = where_to_save)