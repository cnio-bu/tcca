##SCEVAN function in other general script. Pay attention to crc scripts w/ local paths
library("SCEVAN")
library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
metadata <- snakemake@input[["metadata"]]
data_directory <- snakemake@params[["data_dir"]]
cna_dir       <- snakemake@params[["cna_res"]]
where_to_save <- snakemake@output[["seurat_list"]]

## SNAKE params
threads_to_use <- snakemake@threads



## Function definitions
generate_seurat_objects <- function(dgMat, sample) {
    seu <- Seurat::CreateSeuratObject(counts = dgMat,
                              project = sample,
                              assay = "RNA",
    )
    
    return(seu)
}

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

fill_clinical <- function(sc) {
    this_meta <- meta %>%
        filter(SampleID == as.character(unique(sc$orig.ident))) %>%
        uncount(weights = ncol(sc)) %>%
        as.data.frame()
    rownames(this_meta) <- colnames(sc)
    sc <- AddMetaData(object = sc, metadata = this_meta)
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

##CODE
#Load all samples
all_samples <- list.dirs(data_directory,
                         recursive = FALSE,
                         full.names = FALSE
)

full_paths <- paste0(data_directory, "/", all_samples) %>% 
              list.files(full.names = T)


all_mats <- lapply(X = full_paths, FUN = Seurat::Read10X_h5)
names(all_mats) <- all_samples

#Create Seurat objects
all_seurat_objects <- lapply(seq_along(all_mats), function(id) {
    generate_seurat_objects(all_mats[[id]], sample = all_samples[id])
})

#Filter, normalize and scale
filtered_sc <- lapply(all_seurat_objects, filter_sc)
filtered_sc <- lapply(filtered_sc, normalize_and_scale)

#Fill clinical data
meta <- read.csv(
    metadata,
    header=TRUE,
    row.names=NULL,
    sep="\t"
    )

filtered_sc_clinical <- lapply(filtered_sc, fill_clinical)

# Create the dir for CNA
dir.create(cna_dir, showWarnings = FALSE)
setwd(cna_dir)

#Annotate CNA
annotated_sc <- lapply(filtered_sc_clinical, annotate_cna_clones)

#Save
saveRDS(object = annotated_sc, file = where_to_save)
