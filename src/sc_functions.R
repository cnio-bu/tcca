library("Seurat")
library("tidyverse")

load_10x_from_geo <- function(sample){
    
    counts <- Matrix::readMM(paste0(sample, "_matrix.mtx"))
    
    genes <- read_tsv(paste0(sample, "_features.tsv"), col_names = FALSE) %>%
        pull("X2")
    
    cids <- read_tsv(paste0(sample, "_barcodes.tsv"), col_names = FALSE) %>%
        pull("X1")
    
    cids <- gsub(pattern = "-1$", replacement = "", x = cids)
    rownames(counts) <- genes
    colnames(counts) <- cids
    
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts,
                                             project = sample
                                             )
    return(seurat_obj)
}

filter_sc <- function(sc, res_dir) {
    
    this_sample <- unique(sc@meta.data$orig.ident)
    where_to_save <- paste0(res_dir, "/", this_sample)
    sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")
    sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "percent.ribo")
    
    ## get vlnplot for this sample
    this_sc_qc <- VlnPlot(sc,
                          features = c("nFeature_RNA",
                                       "nCount_RNA",
                                       "percent.mt",
                                       "percent.ribo"),
                          ncol = 4)
    
    ggsave(
        plot = this_sc_qc,
        filename = paste0(where_to_save, "_pre_qc.png"),
        dpi = 100,
        height = 7,
        width = 28
        )
    
    sc_filtered <- subset(x = sc, subset = (percent.mt <= 10) &
                              (nFeature_RNA >= 2000 & nFeature_RNA <= 7000) &
                              (nCount_RNA > 500) 
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
    
    this_sc_post_qc <- VlnPlot(new_filtered_sc,
                               features = c("nFeature_RNA",
                                            "nCount_RNA",
                                            "percent.mt",
                                            "percent.ribo"),
                               ncol = 4)
    
    ggsave(
        plot = this_sc_post_qc,
        filename = filename = paste0(where_to_save, "_post_qc.png"), 
        dpi = 100,
        height = 7,
        width = 28
        )
    
    return(new_filtered_sc)
}

normalize_and_scale <- function(sc) {
    sc <- Seurat::NormalizeData(
        sc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
        )
    sc <- Seurat::FindVariableFeatures(sc, selection.method = "vst")
    sc <- Seurat::ScaleData(sc, features = rownames(sc))
    return(sc)
}