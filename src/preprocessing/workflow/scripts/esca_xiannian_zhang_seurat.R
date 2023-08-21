library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file1      <- snakemake@input[["matrix1"]] #CD45- cells
mat_file2      <- snakemake@input[["matrix2"]] #CD45+ cells
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
        mutate(
            barcode = gsub("\\-", "\\.", str_remove_all(string = Cell, pattern = ".*@"))
        ) %>%
        as.data.frame()
    
    common_cells <- intersect(colnames(sc), this_meta$barcode)
    this_meta <- this_meta[this_meta$barcode %in% common_cells, ]
    rownames(this_meta) <- this_meta$barcode
    sc <- sc[, common_cells]
    sc <- AddMetaData(sc, metadata = this_meta)
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


mat1 <- data.table::fread(mat_file1) %>%
    data.frame(row.names = 1)
mat2 <- data.table::fread(mat_file2) %>%
    data.frame(row.names = 1)

metadata <- data.table::fread(input = metadata)  %>%
    as.data.frame()

seu1 <- Seurat::CreateSeuratObject(
    counts = mat1,
    project = "esca"
)

seu2 <- Seurat::CreateSeuratObject(
    counts = mat2,
    project = "esca"
)

seu1_annotated <- annotate_clinical_data(seu1)
seu2_annotated <- annotate_clinical_data(seu2)
seu_list1 <- Seurat::SplitObject(object = seu1_annotated, split.by = "Sample")
seu_list2 <- Seurat::SplitObject(object = seu2_annotated, split.by = "Sample")
seu_list <- c(seu_list1, seu_list2) 

seu_list <- lapply(seu_list, filter_sc)
seu_list <- lapply(seu_list, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
seu_list <- lapply(seu_list, rename_columns, 
                              malignancy_colname = "Celltype..major.lineage.", 
                              malignant_names = c("Malignant"),
                              cell_type_colname = "Celltype..major.lineage.",
                              sample_colname = "Sample", 
                              patient_colname = "Patient")

saveRDS(seu_list, where_to_save)