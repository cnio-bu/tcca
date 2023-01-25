library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
mats <- snakemake@input[["mats"]]
barcodes <- snakemake@input[["barcodes"]]
features <- snakemake@input[["features"]]
cell_types <- snakemake@input[["annotations"]]

where_to_save <- snakemake@output[["seurat_list"]]

dataset <- data.frame(
    mats = mats,
    barcodes = barcodes,
    features = features,
    annotation = cell_types
    )

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

all_seurat_objects = list()

for(i in 1:nrow(dataset)) {
    sample_name <- gsub(
        pattern="_matrix.mtx",
        replacement="",
        x=basename(dataset[i, "mats"])
        )

    counts <- Matrix::readMM(dataset[i, "mats"])
    genes <- read_tsv(dataset[i, "features"], col_names = FALSE, show_col_types = FALSE) %>%
        pull("X2")
    
    cids <- read_tsv(dataset[i, "barcodes"], col_names = FALSE, show_col_types = FALSE) %>%
        pull("X1")
    
    cell_annot <- read.csv(dataset[i, "annotation"],
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
                                             project = sample_name,
                                             meta.data = cell_annot_filtered
                                             )

    all_seurat_objects <- append(all_seurat_objects, seurat_obj)

}
# QC
filtered_sc <- lapply(all_seurat_objects, filter_sc)
# Normalize
filtered_sc <- lapply(filtered_sc, normalize_and_scale)


saveRDS(object = filtered_sc, file = where_to_save)