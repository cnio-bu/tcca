library("tidyverse")
library("Seurat")

## SNAKEMAKE I/O
metadata <- snakemake@input[["metadata"]]
reference_gene_annotation <- snakemake@input[["reference_gene_annotation"]]
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

generate_seurat_objects <- function(dgMat, sample) {
  sample_name <- sample

  ## Annotate to HGNC
  ensembl_genes <- rownames(dgMat)
  genes_to_keep <- gene_annot %>%
  select(symbol, ensembl_gene_id) %>%
  filter(ensembl_gene_id %in% ensembl_genes)

  gene_dict <- genes_to_keep %>%
  select(ensembl_gene_id, symbol) %>%
  deframe()

  annotated_hugo <- gene_dict[ensembl_genes]

  ## Remove NAs
  annotated_hugo <- annotated_hugo[!is.na(annotated_hugo)]

  mat_annot <- dgMat[names(annotated_hugo), ]
  rownames(mat_annot) <- annotated_hugo[rownames(mat_annot)]

  seu <- CreateSeuratObject(counts = mat_annot,
                            project = sample,
                            assay = "RNA",
  )
  return(seu)
}

fill_metadata <- function(sc) {
    
    meta_set <- meta %>%
        filter(Sample == unique(sc$orig.ident))
    rownames(meta_set) <- sapply(strsplit(rownames(meta_set),'@'),'[',2)
    cells_to_keep <- intersect(colnames(sc), rownames(meta_set))
    #We perform our filters and then keep cells in the intersection 
    #between their annotated cells and our seurat object
    sc_subset <- subset(sc, cells = cells_to_keep)
    
    meta_set <- meta_set %>%
        filter(rownames(meta_set) %in% cells_to_keep) %>%
        as.data.frame()
    
    sc_subset <- AddMetaData(object = sc_subset, metadata = meta_set)
    return(sc_subset)
    
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

all_mats <- lapply(X = full_paths, FUN = Seurat::Read10X, gene.column = 1)

##Get gene annotations
gene_annot <- readr::read_tsv(reference_gene_annotation)

## Get Seurat objects
all_seurat_objects <- lapply(seq_along(all_mats), function(id) {
    generate_seurat_objects(all_mats[[id]], sample = all_samples[id])
})

## Filter cells
filtered_sc <- lapply(all_seurat_objects, filter_sc)

## Normalize and scale data
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Annotation
#Load metadata from the study and annotate seurat objects
print(metadata)
meta <- read.table(metadata, header=TRUE, row.names = 1, sep ="\t")

filtered_sc_clinical <- lapply(filtered_sc, fill_metadata)

## Seurat object
saveRDS(object = filtered_sc_clinical,
        file = where_to_save
        )

