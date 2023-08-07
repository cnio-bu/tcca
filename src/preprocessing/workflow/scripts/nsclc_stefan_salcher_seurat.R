library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat_file      <- snakemake@input[["matrix"]]
reference_gene_annotation <- snakemake@input[["reference_gene_annotation"]]
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

seu <- readRDS(mat_file)
DefaultAssay(seu) <- "RNA"

gene_annot <- readr::read_tsv(reference_gene_annotation)

## Extract count mat
mat <- seu@assays$RNA@counts
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
metadata <- seu@meta.data

## Free mem
rm(seu)
gc()

## Regenerate Seurat object
full_seu <- Seurat::CreateSeuratObject(
    counts = mat_annot,
    meta.data = metadata,
    assay = "RNA"
    )

## Set idents to sample
full_seu$orig.ident <- full_seu$sample

seu_list <- Seurat::SplitObject(object = full_seu, split.by = "sample")
names(seu_list) <- unique(full_seu$sample)
samples_to_filter <- c("Leader_Merad_2021_336", "Leader_Merad_2021_298")

seu_list <- seu_list[!(names(seu_list) %in% samples_to_filter)]
seu_list <- lapply(seu_list, filter_sc)

## Get rid of NULLs samples w/o malignants left
seu_list[sapply(seu_list, is.null)] <- NULL
## In these datasets, some samples had very few cells, less than the minimum
## for calculating variable features
seu_list <- seu_list[sapply(seu_list, ncol) > 10]

seu_list <- lapply(seu_list, normalize_and_scale)

saveRDS(seu_list, where_to_save)