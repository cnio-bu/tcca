library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
data_directory <- snakemake@params[["data_dir"]]
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

annotate_clinical_data <- function(sc){
    this_meta <- metadata %>%
        filter(Sample == unique(sc$orig.ident)) %>%
        mutate(
            barcode = str_remove_all(string = Cell, pattern = "^GSM.*@")
        ) %>%
        as.data.frame()
    
    common_cells <- intersect(colnames(sc), this_meta$barcode)
    this_meta <- this_meta[this_meta$barcode %in% common_cells, ]
    rownames(this_meta) <- this_meta$barcode
    sc <- sc[, common_cells]
    sc <- AddMetaData(sc, metadata = this_meta)
    return(sc)
    
}

rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname, patient_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE))
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  colnames(sc@meta.data)[colnames(sc@meta.data) == patient_colname] <- "patient"
  
  return(sc)
}

all_samples <- list.files(data_directory, full.names = TRUE)
seu_list <- list()

for(sample in all_samples){
    sample_name <- basename(sample)
    mat <- Seurat::Read10X(data.dir = sample)
    seu <- Seurat::CreateSeuratObject(
        counts = mat,
        assay = "RNA",
        project = sample_name
    )
    seu_list <- c(seu_list, seu)
}

names(seu_list) <- basename(all_samples)

metadata <- data.table::fread(input = metadata)  %>%
    mutate(
        Patient = str_remove_all(string = Patient, pattern = "^patient: ")
    ) %>%
    as.data.frame()

seu_list <- lapply(seu_list, filter_sc)
seu_list <- lapply(seu_list, normalize_and_scale)
seu_list <- lapply(seu_list, annotate_clinical_data)
## Add and rename standarized columns: malignancy, cell_type, sample, patient
seu_list <- lapply(seu_list, rename_columns, 
                              malignancy_colname = "Celltype..major.lineage.", 
                              malignant_names = c("Malignant"),
                              cell_type_colname = "Celltype..major.lineage.",
                              sample_colname = "Sample", 
                              patient_colname = "Patient")


saveRDS(seu_list, where_to_save)