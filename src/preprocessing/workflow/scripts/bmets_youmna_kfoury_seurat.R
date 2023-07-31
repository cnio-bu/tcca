library(Seurat)
library(tidyverse)

## SNAKEMAKE I/O
metadata <- snakemake@input[["metadata"]]
data_directory <- snakemake@params[["data_dir"]]
where_to_save <- snakemake@output[["seurat_list"]]


## FUNCTIONS DEFINITION
#Generate seurat objects
generate_seurat_objects <- function(dgMat, sample) {
  dgMat <- dgMat[!duplicated(dgMat$X),] #Remove duplicate rows
  rownames(dgMat) <- dgMat$X
  dgMat <- dgMat[,-1]
  seu <- CreateSeuratObject(counts = dgMat,
                            project = sample,
                            assay = "RNA",
                            names.delim = "z" #Names delim is set to a random character in order not to set original ident by barcode id code
  )
  
  return(seu)
}

#Standard filtering, normalization and scaling
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

#Annotation with metadata from the authors (which includes copykat predictions)
fill_metadata <- function(sc) {
  meta_set <- meta %>%
    filter(sample == unique(sc$orig.ident)) %>%
    as.data.frame()
  meta_set <-  mutate(meta_set,
    barcode = gsub("\\-", "\\.", meta_set$barcode)
  )
  cells_to_keep <- intersect(colnames(sc), meta_set$barcode)
  #We perform our filters and then keep cells in the intersection 
  #between their annotated cells and our seurat object
  sc_subset <- subset(sc, cells = cells_to_keep)
  
  meta_set <- meta_set %>%
    filter(barcode %in% cells_to_keep) %>%
    as.data.frame()
  rownames(meta_set) <- meta_set$barcode
  meta_set$barcode <- NULL
  
  sc_subset <- AddMetaData(object = sc_subset, metadata = meta_set)
  return(sc_subset)
  
}


##CODE
#Load all samples
all_samples <- list.files(data_directory,
                         recursive = FALSE,
                         full.names = FALSE
)

full_paths <- paste0(data_directory,
                     "/",
                     all_samples
)

all_samples <- lapply(all_samples, function(x){
  gsub(".count.csv", "", x)
})

class(all_samples) <- "character"

all_mats <- lapply(X = full_paths, function(x){
  read.table(x, header = T, sep = ",")
  })
names(all_mats) <- all_samples

#Generate Seurat objects
all_seurat_objects <- lapply(seq_along(all_mats), function(id) {
  generate_seurat_objects(all_mats[[id]], sample = all_samples[id])
})

# QC
filtered_sc <- lapply(all_seurat_objects, filter_sc)

# Normalize
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Annotation
#Load metadata from the study and annotate seurat objects
meta <- read.table(metadata, header=TRUE, sep =",")
meta$sample <- gsub("^(.*?)_.*$", "\\1", meta$barcode)

filtered_sc_clinical <- lapply(filtered_sc, fill_metadata)

#Save
saveRDS(object = filtered_sc_clinical, file = where_to_save)
