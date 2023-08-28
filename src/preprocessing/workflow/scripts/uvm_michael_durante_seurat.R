library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file      <- snakemake@params[["data_dir"]]
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

mat <- Seurat::Read10X(data.dir = mat_file)

metadata <- data.table::fread(input = metadata) %>%
    mutate(
    cell_id = gsub(pattern = "^.*@", replacement = "", x = Cell)
  ) %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  as.data.frame()

metadata$Cell <- NULL
#Add patient column (only one sample per patient)
metadata <- mutate(metadata, sample = Patient) 

## Keep cells with annot only
cells_to_keep <- intersect(colnames(mat), metadata$cell_id)
mat <- mat[, cells_to_keep]
metadata <- metadata[metadata$cell_id %in% cells_to_keep, ]
rownames(metadata) <- metadata$cell_id


seu <- Seurat::CreateSeuratObject(
    counts = mat,
    meta.data = metadata,
    project = "uvm",
    assay = "RNA"
    )

seu_list <- Seurat::SplitObject(object = seu, split.by = "Patient")

seu_list <- lapply(seu_list, filter_sc)
seu_list <- lapply(seu_list, normalize_and_scale())
## Add and rename standarized columns: malignancy, cell_type, sample, patient
seu_list <- lapply(seu_list, rename_columns, 
                              malignancy_colname = "Celltype..major.lineage.",
                              malignant_names = c("Malignant"),
                              cell_type_colname = "Celltype..major.lineage.",
                              sample_colname = "sample",
                              patient_colname = "patient")

saveRDS(seu_list, where_to_save)