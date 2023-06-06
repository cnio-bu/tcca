library("SCEVAN")
library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
data_directory <- snakemake@params[["data_dir"]]
cna_dir       <- snakemake@params[["cna_res"]]
where_to_save <- snakemake@output[["seurat_list"]]

## SNAKE params
threads_to_use <- snakemake@threads



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

annotate_cna_clones <- function(sc){
    
    this_mat <- Seurat::GetAssayData(sc, slot = "counts")
    cna_pred <- SCEVAN::pipelineCNA(
        count_mtx = this_mat,
        sample = unique(sc$orig.ident),
        par_cores = threads_to_use,
        SUBCLONES = TRUE,
        plotTree = FALSE,
        organism = "human",
        SCEVANsignatures = TRUE
        
    )
    sc <- AddMetaData(object = sc, metadata = cna_pred)
    return(sc)
}

all_samples <- list.files(data_directory, full.names = TRUE)
all_seurat_objects <- list()

for(sample in all_samples){
    sample_name <- basename(sample)
    mat <- Seurat::Read10X(data.dir = sample)
    seu <- Seurat::CreateSeuratObject(
        counts = mat,
        assay = "RNA",
        project = sample_name
        )
    
    all_seurat_objects <- c(all_seurat_objects, seu)
                                    
}

names(all_seurat_objects) <- basename(all_samples)
filtered_sc <- lapply(all_seurat_objects, filter_sc)
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

## sample 18 fails in SCEVAN, temporarily removed
filtered_sc <- filtered_sc[-17]

## Create the dir for CNA
dir.create(cna_dir, showWarnings = FALSE)
setwd(cna_dir)
annotated_sc <- lapply(filtered_sc, annotate_cna_clones)

saveRDS(object = annotated_sc, file = where_to_save)
