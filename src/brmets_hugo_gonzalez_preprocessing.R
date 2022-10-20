library("tidyverse")
library("Seurat")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

setwd(cfg["working_dir"])

## Code for the preprocessing of brmets_hugo_gonzalez ##

## Not sourcing it since some black magic is necessary to annotate Cell Types
load_10x_from_geo <- function(sample){
    
    counts <- Matrix::readMM(paste0(sample, "_matrix.mtx"))
    
    genes <- read_tsv(paste0(sample, "_features.tsv"), col_names = FALSE) %>%
        pull("X2")
    
    cids <- read_tsv(paste0(sample, "_barcodes.tsv"), col_names = FALSE) %>%
        pull("X1")
    
    cell_annot <- read.csv(paste0(sample, "_Cell_Types_Annotations.csv"),
                           row.names = 1,
                           stringsAsFactors = FALSE
    )
    
    ## Integrated re-annot. is marked as s1/s2. Keep first the s2 if conflicts
    ## arise. If not, keep S1. Then, keep the unintegrated.
    cell_annot_filtered <- cell_annot %>%
        mutate(
            reported_barcode = rownames(cell_annot),
            Cell_Type = stringr::str_replace_all(Cell_Type, pattern = "MTCs", replacement = "MTC"),
            Cell_Type = as_factor(Cell_Type),
            Cell_Type = fct_relevel(Cell_Type, "MTC", after = 0)
        ) %>%
        separate(reported_barcode,
                 sep = "_",
                 into = c("integration", "samp", "barcode"),
                 remove = TRUE,
                 fill = "left"
        ) %>%
        arrange(desc(Cell_Type)) %>%
        distinct(barcode, .keep_all = TRUE) %>%
        select(barcode, Cell_Type) %>%
        as.data.frame()
    
    rownames(cell_annot_filtered) <- cell_annot_filtered$barcode
    
    cids <- gsub(pattern = "-1$", replacement = "", x = cids)
    rownames(counts) <- genes
    colnames(counts) <- cids
    
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts,
                                             project = sample,
                                             meta.data = cell_annot_filtered)
    return(seurat_obj)
}

filter_sc <- function(sc) {
    this_sample <- unique(sc@meta.data$orig.ident)
    where_to_save <- paste0(getwd(), "/single_cell/qc/brmets_hugo_gonzalez/")
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
        filename = paste0(where_to_save, this_sample, "_pre_qc.png"),
        dpi = 100,
        height = 7,
        width = 24
        )
    
    sc_filtered <- subset(x = sc, subset = (percent.mt <= 10) &
                              (nFeature_RNA >= 2000 & nFeature_RNA <= 7000) &
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
    
    this_sc_post_qc <- VlnPlot(new_filtered_sc,
                               features = c("nFeature_RNA",
                                            "nCount_RNA",
                                            "percent.mt",
                                            "percent.ribo"),
                               ncol = 4)
    
    ggsave(plot = this_sc_post_qc,
           filename = paste0(where_to_save, this_sample, "_post_qc.png"),
           dpi = 100,
           height = 7,
           width = 24
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

keep_all_malignants <- function(sc) {
    sc_filtered <- subset(x = sc, subset = Cell_Type == "MTC")
    return(sc_filtered)
}

## get all mats
mats <- list.files(
    path = paste0(getwd(), "/single_cell/raw/brmets_hugo_gonzalez/"),
    pattern = "_matrix.mtx",
    full.names = FALSE
)

mats <- gsub(pattern = "_matrix.mtx", replacement = "", x = mats)
mats <- mats[!grepl(pattern = "*_Mouse_*", x = mats)]
mats <- mats[mats != "GSM5645906_Rhabdomyosarcoma"] ## missing cell types

# Load data
raw_seurat_list <- lapply(mats, load_10x_from_geo)
names(raw_seurat_list) <- mats

# QC
filtered_sc <- lapply(raw_seurat_list, filter_sc)
# Normalize
filtered_sc <- lapply(filtered_sc, normalize_and_scale)
# Keep malignant cells
malignant_sc <- lapply(filtered_sc, keep_all_malignants)

res_dir <- paste0(getwd(), "/single_cell/obj/brmets_hugo_gonzalez/")

saveRDS(
    object = filtered_sc,
    file = paste0(res_dir, "all_samples_filtered.rds")
    )
saveRDS(
    object = malignant_sc,
    file = paste0(res_dir, "all_malignant.rds")
    )
