library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
metadata_directory <- snakemake@input[["metadata"]]
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

## Function definitions
rename_columns <- function(sc, malignancy_colname, malignant_names, cell_type_colname, sample_colname, patient_colname){
  sc@meta.data <- sc@meta.data %>%
    mutate(malignancy = ifelse(sc@meta.data[, malignancy_colname] %in% malignant_names, TRUE, FALSE))
  
  colnames(sc@meta.data)[colnames(sc@meta.data) == cell_type_colname] <- "cell_type"
  colnames(sc@meta.data)[colnames(sc@meta.data) == sample_colname] <- "sample"
  colnames(sc@meta.data)[colnames(sc@meta.data) == patient_colname] <- "patient"
  
  return(sc)
}


## Read data
all_cells <- Seurat::Read10X(data.dir = data_directory,
                             gene.column = 1
                             )

## Read metadata
metadata <- read.delim(sep="\t",row.names = 1,
                       file = metadata_directory
                       )

## Merge data + metadata
seu <- Seurat::CreateSeuratObject(counts = all_cells,
                                  project = "rcell_kevin_bi",
                                  meta.data = metadata
                                  )

## Fix missing donor and biosample annotations
seu@meta.data <- seu@meta.data %>%
  mutate(donor_id = toupper(str_remove(rownames(seu@meta.data),
                            pattern = "[A-Z]+\\.")),
         biosample_id = paste0(donor_id, "_scRNA")
         )

## Switch ident to biosample
seu$old_ident <- seu$orig.ident
seu$orig.ident <- seu$biosample_id

## Split the merged obj
sample_list <- Seurat::SplitObject(object = seu, split.by = "biosample_id")
names(sample_list) <- sapply(sample_list,
                             function(sc){unique(sc$"biosample_id")}
                             )

## Filter cells
filtered_sc <- lapply(sample_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "FinalCellType", 
                              malignant_names = c("TP1", "TP2", "Cycling Tumor"),
                              cell_type_colname = "FinalCellType",
                              sample_colname = "biosample_id", 
                              patient_colname = "donor_id")

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )
