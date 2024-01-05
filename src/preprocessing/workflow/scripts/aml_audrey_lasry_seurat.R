library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file      <- snakemake@input[["matrix"]]
features      <- snakemake@input[["features"]]
cells         <- snakemake@input[["cells"]]
metadata      <- snakemake@input[["metadata"]]
aditional_data <- snakemake@input[["additional_metadata"]]

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

## Function definitions
rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname, patient_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE))
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  colnames(sc@meta.data)[colnames(sc@meta.data) == patient_colname] <- "patient"
  
  return(sc)
}

mat <- Matrix::readMM(file = mat_file)

cells <- data.table::fread(
    input = cells,
    header = FALSE
    )
genes <- data.table::fread(
    input = features,
    header = FALSE
    )

rownames(mat) <- genes$V1
colnames(mat) <- cells$V1

meta.data <- data.table::fread(
    input = metadata
    ) %>%
    select(-X, -Y) %>%
    as.data.frame()

rownames(meta.data) <- meta.data$NAME
meta.data$NAME <- NULL

seu <- Seurat::CreateSeuratObject(counts = mat,
                                  project = "aml_audrey_lasry",
                                  meta.data = meta.data
                                  )

rm(feats, genes, mat, meta.data)
gc()

additional_meta.data <- data.table::fread(aditional_data) %>%
    filter(
        NAME %in% colnames(seu)
    ) %>%
    select(
        NAME,
        donor_id,
        biosample_id,
        sex,
        organ__ontology_label,
        disease__ontology_label
        ) %>%
    as.data.frame()

rownames(additional_meta.data) <- additional_meta.data$NAME
additional_meta.data$NAME <- NULL

seu <- AddMetaData(object = seu, metadata = additional_meta.data)
seurat_list <- Seurat::SplitObject(object = seu, split.by = "biosample_id")

filtered_sc <- lapply(seurat_list, filter_sc)
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "malignant", 
                              malignant_names = c("malignant"),
                              cell_type_colname = "Cell_type_identity",
                              sample_colname = "biosample_id", 
                              patient_colname = "donor_id"
                              )

saveRDS(object = filtered_sc, file = where_to_save)
