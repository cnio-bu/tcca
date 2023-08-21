library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
seurat_object_list <- snakemake@input[["object_list"]]
where_to_save <- snakemake@output[["seurat_list"]]

## Load the source file... big RData obj
seurat_object <- readRDS(seurat_object_list)

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

## Function definitions
rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE))
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  
  sc@meta.data <- mutate(sc@meta.data,
  patient = gsub("^(.*?)\\..*$", "\\1", sample))

  return(sc)
}

## This is a seurat obj. with an unknown v. Converting it to Seurat 4.
mat <- GetAssayData(seurat_object, slot = "counts")
meta.data <- seurat_object@meta.data

synovial_object <- Seurat::CreateSeuratObject(counts = mat, project = "synovial",
                                       meta.data = meta.data)

# Change default ident for later split and QC
synovial_object@meta.data$orig.ident <- synovial_object@meta.data$characteristics..sample

## Split the merged obj
samples_list <- Seurat::SplitObject(object = synovial_object,
                                               split.by = "characteristics..sample"
                                               )

names(samples_list) <- unique(synovial_object@meta.data$characteristics..sample)

## Filter cells
filtered_sc <- lapply(samples_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Get rid of "NULL" samples, there are no cells left in these
filtered_sc[sapply(filtered_sc, is.null)] <- NULL

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "characteristics..cell.type", 
                              malignant_names = c("Malignant"),
                              cell_type_colname = "characteristics..cell.type",
                              sample_colname = "characteristics..sample")

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )
        