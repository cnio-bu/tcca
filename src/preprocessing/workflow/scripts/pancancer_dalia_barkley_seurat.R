library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
seurat_object_list <- snakemake@input[["object_list"]]
where_to_save <- snakemake@output[["seurat_list"]]

## Load the source file... big RData obj
load(seurat_object_list)

## Change the default assay to RNA for all samples
change_to_rna <- function(sc) {
    DefaultAssay(sc) <- "RNA"
    return(sc)
}

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

rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(
        malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE),
        ) %>%
      select(-sample)
    
  if("cell_type" %in% colnames(sc@meta.data)){
      sc@meta.data <- sc@meta.data %>%
          select(-cell_type)
  }
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  sc@meta.data <- sc@meta.data %>%
    mutate(patient = sample)
  
  return(sc)
}

srt.list.primary.all <- lapply(srt.list.primary.all, change_to_rna)

filtered_sc <- lapply(srt.list.primary.all, filter_sc)
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Get rid of "NULL" samples, there are no cells left in these
filtered_sc[sapply(filtered_sc, is.null)] <- NULL


## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "type", 
                              malignant_names = c("malignant"),
                              cell_type_colname = "pop",
                              sample_colname = "orig.ident")

saveRDS(
    object = filtered_sc,
    file = where_to_save
    )
