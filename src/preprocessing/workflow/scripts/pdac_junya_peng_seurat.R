library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
seurat_object_list <- snakemake@input[["object_list"]]
celltype <- snakemake@input[["celltype"]]
where_to_save <- snakemake@output[["seurat_list"]]

## Load the source file
load(seurat_object_list)

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

#Add anotation to cells
anot <- read.table(celltype,
                   sep = "\t",
                   header = T,
                   row.names = 1
                   )

sc_Peng_PDAC@meta.data <- merge(sc_Peng_PDAC@meta.data,
                                anot,
                                by = "row.names",
                                all.x = T
                                )

row.names(sc_Peng_PDAC@meta.data) <- sc_Peng_PDAC@meta.data$Row.names
sc_Peng_PDAC@meta.data <- sc_Peng_PDAC@meta.data[,-1]

# Change default ident for later split and QC
sc_Peng_PDAC@meta.data$orig.ident <- sc_Peng_PDAC@meta.data$Sample
Idents(sc_Peng_PDAC) <- "orig.ident"

## Split the merged obj
samples_list <- Seurat::SplitObject(object = sc_Peng_PDAC,
                                    split.by = "orig.ident"
                                    )

## Filter cells
filtered_sc <- lapply(samples_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Get rid of "NULL" samples, there are no cells left in these
filtered_sc[sapply(filtered_sc, is.null)] <- NULL

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "cluster", 
                              malignant_names = c("Ductal cell type 2"),
                              cell_type_colname = "cluster",
                              sample_colname = "orig.ident", 
                              patient_colname = "Sample")

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )