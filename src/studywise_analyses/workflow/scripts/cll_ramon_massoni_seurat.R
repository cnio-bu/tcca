library("tidyverse")
library("Seurat")
library("Matrix")

## SNAKEMAKE I/O
data_directory <- snakemake@params[["data_dir"]]
where_to_save <- snakemake@output[["seurat_list"]]

## Function definitions
filter_sc <- function(sc) {
    sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
    sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "percent.ribo")
    
    sc_filtered <- subset(x = sc, subset = (percent.mt <= 15) &
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

## Get all samples
all_samples <- list.dirs(data_directory,
                         recursive = FALSE,
                         full.names = FALSE
                         )

## Get Seurat objects

### Demultiplexing process from https://github.com/massonix/sampling_artifacts/blob/master/2-CLL/01-demultiplexing.Rmd
# Demultiplex
# Load expression matrix, gene and cell metadata
cll_list <- list()
for (sample in all_samples) {
  expression_matrix <- readMM(str_c(data_directory, "/", sample, "/", "GSE132065_", "matrix_", sample, ".mtx.gz"))
  barcodes <- read_csv(str_c(data_directory, "/", sample, "/", "GSE132065_", "barcodes_", sample, ".tsv.gz"), col_names = FALSE)
  colnames(barcodes) <- "barcode"
  features <- read_tsv(str_c(data_directory, "/", sample, "/", "GSE132065_", "features_", sample, ".tsv.gz"), col_names = FALSE)
  colnames(features) <- c("ensembl", "symbol", "feature_type")
  rownames(expression_matrix) <- features$symbol
  colnames(expression_matrix) <- barcodes$barcode

  # Separate HTO and RNA matrices
  hto_ind <- which(str_detect(features$feature_type, "Antibody Capture"))
  rna_ind <- which(str_detect(features$feature_type, "Gene Expression"))
  cll_hto <- expression_matrix[hto_ind, ]
  cll_rna <- expression_matrix[rna_ind, ]

  # Setup Seurat object
  cll <- CreateSeuratObject(counts = cll_rna)

  # Normalize RNA data with log normalization
  cll <- NormalizeData(cll)

  # Find and scale variable features
  cll <- FindVariableFeatures(cll, selection.method = "vst")
  cll <- ScaleData(cll, features = VariableFeatures(cll))

  # Add HTO as an independent assay
  cll[["HTO"]] <- CreateAssayObject(counts = cll_hto)
  cll <- NormalizeData(cll, assay = "HTO", normalization.method = "CLR")

  # Demultiplex
  cll <- HTODemux(cll, assay = "HTO", positive.quantile = 0.99)

  # Append to list of Seurat objects
  cll_list[[sample]] <- cll
}

# Save demultiplexed Seurat object
cll_merged <- merge(
  cll_list$`1472_4C`,
  y = c(cll_list$`1892_4C`),
  add.cell.ids = all_samples,
  project = "cll_ramon_massoni"
)

# Recode and retain important variables
cll_merged$donor <- str_remove(colnames(cll_merged), "_.*$")
cll_merged$time <- str_remove(cll_merged$hash.ID, "....-..-")
cll_merged$temperature <- cll_merged$hash.ID %>%
  str_remove("^....-") %>%
  str_remove("-.+h$")
selection <- c("nCount_RNA", "nFeature_RNA", "HTO_classification", "hash.ID",
               "time", "donor", "temperature", "orig.ident")
cll_merged@meta.data <- cll_merged@meta.data[, selection]
###

## Subset to remove Negative and Doublet cells
cll_merged_subset <- subset(x = cll_merged, subset = hash.ID != "Negative" & hash.ID != "Doublet")

## Set donor_hour label as orig.ident
cll_merged_subset$orig.ident <- cll_merged$hash.ID

## Annotate cell type

### Replicate annotation from https://github.com/massonix/sampling_artifacts/blob/master/2-CLL/03-annotation.Rmd
#Find Variable Genes
cll_merged_subset <- FindVariableFeatures(cll_merged_subset)
#Scale data
cll_merged_subset <- ScaleData(cll_merged_subset)
#PCA
cll_merged_subset <- RunPCA(cll_merged_subset, features = VariableFeatures(cll_merged_subset))
#Determine statistically significant principal components with elbow plot:
ElbowPlot(cll_merged_subset)
#Cluster cells: same dimensions and resolutions than in the authors workflow
cll_merged_subset <- FindNeighbors(cll_merged_subset, dims = 1:14)
cll_merged_subset <- FindClusters(cll_merged_subset, resolution = 0.01)
#Dimensionality reduction
cll_merged_subset <- RunUMAP(cll_merged_subset, reduction = "pca", dims = 1:14)
###

## Assigning cell type identity to clusters
cll_merged_subset$cell_type <- ifelse(cll_merged_subset$seurat_clusters == 0 | cll_merged_subset$seurat_clusters == 1, "CLL", ifelse(cll_merged_subset$seurat_clusters == 2, "T and NK", ifelse(cll_merged_subset$seurat_clusters == 3, "Monocytes","")))

## Split objects
sample_list <- Seurat::SplitObject(object = cll_merged_subset, split.by = "orig.ident")

## Filter cells
filtered_sc <- lapply(sample_list, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "cell_type", 
                              malignant_names = c("CLL"),
                              cell_type_colname = "cell_type",
                              sample_colname = "orig.ident", 
                              patient_colname = "donor")

## Seurat object
saveRDS(object = filtered_sc,
        file = where_to_save
        )
