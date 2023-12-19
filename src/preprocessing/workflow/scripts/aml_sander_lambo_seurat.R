library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
data_dir       <- snakemake@params[["data_dir"]]
geo_to_samples <- snakemake@input[["geo_to_samples"]]
where_to_save  <- snakemake@output[["seurat_list"]]


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

sample_dirs <- list.dirs(
  path = data_dir,
  recursive = FALSE
)

samples <- basename(sample_dirs)

## load metadata
meta.files <- list.files(
  path = data_dir,
  pattern = "^GSM.*.tsv",
  recursive = FALSE,
  full.names = TRUE
  )


meta.data <- meta.files %>%
  map(read_tsv) %>%
  bind_rows()


## add geo ids
geos <- data.table::fread(geo_to_samples)

meta.data <- meta.data %>%
  left_join(y = geos, by = c("GEO_ID" = "sample_id"))

seu_list <- list()  

for (i in c(1:length(samples))){
  mat_dir <- sample_dirs[i]
  sample_name <- samples[i]
  
  mat <- Seurat::Read10X(data.dir = mat_dir, gene.column = 1)
  
  sample_meta <- meta.data %>%
    filter(gsm_id == sample_name) %>%
    as.data.frame()
  
  rownames(sample_meta) <- sample_meta$Cell_Barcode
  
  seu <- CreateSeuratObject(
    counts = mat,
    project = "aml_sander",
    assay = "RNA",
    meta.data = sample_meta
    )
  seu_list <- c(seu_list, seu)
}

## Filter low QC cells
seu_filtered <- lapply(seu_list, FUN = filter_sc)
## Normalize
seu_filtered <- lapply(seu_list, FUN = normalize_and_scale)
## Add and rename standarized columns: malignancy, cell_type, sample, patient
seu_filtered <- lapply(seu_filtered,
                       rename_columns, 
                       malignancy_colname = "Malignant",
                               malignant_names = c("Malignant cells"),
                               cell_type_colname = "Classified_Celltype",
                               sample_colname = "GEO_ID",
                               patient_colname = "Patient_ID"
                       )


saveRDS(object = seu_filtered, file = where_to_save)