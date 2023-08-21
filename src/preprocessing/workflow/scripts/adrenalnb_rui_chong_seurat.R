library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
object_list <- snakemake@input[["object_list"]]
reference_gene_annotation <- snakemake@input[["reference_gene_annotation"]]
tumor_data <- snakemake@input[["tumor_data"]]
gland_data <- snakemake@input[["gland_data"]]

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

load(object_list)

gene_annot <- readr::read_tsv(reference_gene_annotation)

## Extract count mat
mat <- human_NB_subset_tumor@assays$SCT@counts
ensembl_genes <- rownames(mat)

## Annotate to HGNC
genes_to_keep <- gene_annot %>%
select(symbol, ensembl_gene_id) %>%
filter(ensembl_gene_id %in% ensembl_genes)

gene_dict <- genes_to_keep %>%
select(ensembl_gene_id, symbol) %>% 
deframe()

annotated_hugo <- gene_dict[ensembl_genes]

## Remove NAs
annotated_hugo <- annotated_hugo[!is.na(annotated_hugo)]

mat_annot <- mat[names(annotated_hugo), ]
rownames(mat_annot) <- annotated_hugo[rownames(mat_annot)]

## Extract the original meta.data
metadata <- human_NB_subset_tumor@meta.data

## Free mem
rm(human_NB_subset_tumor)
gc()

## Regenerate Seurat object
full_seu <- Seurat::CreateSeuratObject(counts = mat_annot, meta.data = metadata)

## Set idents to sample
full_seu$orig.ident <- full_seu$Sample

## Add missing metadata
tumor_data <- readr::read_csv(tumor_data) %>%
select(cellname, celltype) %>%
as.data.frame()

rownames(tumor_data) <- tumor_data$cellname 
tumor_data$cellname <- NULL

gland_data <- readr::read_csv(gland_data) %>% 
select(cell_id, annotation) %>% 
as.data.frame()

rownames(gland_data) <- gland_data$cell_id
gland_data$cell_id <- NULL 

full_seu <- AddMetaData(full_seu, metadata = tumor_data)
full_seu <- AddMetaData(full_seu, metadata = gland_data)


## Split the merged obj
samples_list <- Seurat::SplitObject(object = full_seu, 
                                    split.by = "Sample"
                                    )
filtered_sc <- lapply(samples_list, filter_sc)
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## Add and rename standarized columns: malignancy, cell_type, sample, patient
filtered_sc <- lapply(filtered_sc, rename_columns, 
                              malignancy_colname = "Cell_annotation", 
                              malignant_names = c("tumor"),
                              cell_type_colname = "celltype",
                              sample_colname = "orig.ident", 
                              patient_colname = "Sample"
                              )


saveRDS(
    object = filtered_sc,
    file = where_to_save
)