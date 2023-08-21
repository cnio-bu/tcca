library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
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

rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname, patient_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE))
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  colnames(sc@meta.data)[colnames(sc@meta.data) == patient_colname] <- "patient"
  
  return(sc)
}


generate_seurat_objects <- function(dgMat, sample, md_path) {
  sample_name <- sample
  metadata <- read.delim(sep = ",",
                         row.names = 1,
                         file = md_path
                         )
  seu <- CreateSeuratObject(counts = dgMat,
                            project = sample,
                            assay = "RNA",
                            meta.data = metadata
  )
  return(seu)
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

metadata_paths <- list.files(path = data_directory,
                             recursive = TRUE,
                             pattern = "*metadata.csv",
                             full.names = TRUE)

all_mats <- lapply(X = full_paths, FUN = Seurat::Read10X, gene.column = 1)

## Get Seurat objects
all_seurat_objects <- lapply(seq_along(all_mats), function(id) {
    generate_seurat_objects(all_mats[[id]], sample = all_samples[id], md_path = metadata_paths[id])
})

## Split objects
sample_list <- lapply(seq_along(all_seurat_objects), function(id) {
  Seurat::SplitObject(object = all_seurat_objects[[id]], split.by = "orig.ident")
  }) %>% unlist()

## Filter cells
filtered_sc <- lapply(sample_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "CellType", 
                              malignant_names = c("Cancer"),
                              cell_type_colname = "CellType",
                              sample_colname = "orig.ident", 
                              patient_colname = "PatientNumber")

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )
