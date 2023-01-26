library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mat <- snakemake@input[["mat_object"]]
annotations <- snakemake@input[["annotations"]]
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

## load full mat
full_mat <- readRDS(file = mat)
cell_annot <- read.delim(annotations)
rownames(cell_annot) <- cell_annot$Index

all_seurat <- Seurat::CreateSeuratObject(
    counts = full_mat,
    project = "all_lung",
    assay = "RNA",
    meta.data = cell_annot
    )

# met. Brain, early tumor, pleural efussion, met. lympn node, late stage tumour"
is_tumor <- c("mBrain", "tLung", "PE", "mLN", "tL/B")

all_seurat <- subset(all_seurat, subset = Sample_Origin %in% is_tumor)

# Reset idents
all_seurat$old_ident <- all_seurat$orig.ident
all_seurat$orig.ident <- all_seurat$Sample

# Create a named list of samples.
samples <- unique(all_seurat$Sample)
samples_list <- Seurat::SplitObject(object = all_seurat, split.by = "Sample")

names(samples_list) <- samples

# Filter, normalize, scale and get most variable feats
samples_list <- lapply(samples_list, filter_sc)
samples_list <- lapply(samples_list, normalize_and_scale)

## Get rid of the NULL elements
samples_list[sapply(samples_list, is.null)] <- NULL

saveRDS(object = samples_list, file = where_to_save)
