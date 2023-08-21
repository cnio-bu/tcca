library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file <- snakemake@input[["mat"]]
annotations <- snakemake@input[["annotations"]]
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


## Load matrix
mat <- data.table::fread(mat_file, sep = "\t")
rownames(mat) <- mat$gene
mat$gene <- NULL 

## Load metadata
annot <- data.table::fread(annotations, sep = "\t")
rownames(annot) <- annot$V1
annot$V1 <- NULL 


full_seu <- CreateSeuratObject(counts = mat, 
                               project = "rcell_r_li",
                               assay = "RNA",
                               meta.data = annot
                               )

rm(mat, annot)
gc()

## Split Seurat object by patient
seurat_list <- Seurat::SplitObject(
    object = full_seu,
    split.by = "patient"
    )

rm(full_seu)
gc()

# QC
filtered_sc <- lapply(seurat_list, filter_sc)
# Normalize
filtered_sc <- lapply(filtered_sc, normalize_and_scale)
## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "broad_type", 
                              malignant_names = c("RCC"),
                              cell_type_colname = "broad_type",
                              sample_colname = "orig.ident", 
                              patient_colname = "patient")

## save the resulting list
saveRDS(
    object = filtered_sc,
    file = where_to_save
)