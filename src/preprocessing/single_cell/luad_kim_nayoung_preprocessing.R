library("tidyverse")
library("Seurat")

## load common funs
source("src/sc_functions.R")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")

## load full mat
full_mat <- Matrix::readHB(paste0(
    data_folder,
    "/raw/luad_kim_nayoung/raw_umi_matrix.mtx")
)

full_mat <- readRDS(file = paste0(
    data_folder,
    "/raw/luad_kim_nayoung/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
    )

cell_annot <- read.delim(paste0(data_folder,
                                "/raw/luad_kim_nayoung/GSE131907_Lung_Cancer_cell_annotation.txt")
                         )


all_seurat <- Seurat::CreateSeuratObject(
    counts = full_mat,
    project = "all_lung",
    assay = "RNA",
    meta.data = cell_annot
    )

rm(full_mat, cell_annot)
gc()

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

# Keep freeing mem.
rm(all_seurat)
gc()

# Filter, normalize, scale and get most variable feats
samples_list <- lapply(samples_list, 
                       filter_sc,
                       res_dir = paste0(data_folder, "/qc/luad_kim_nayoung/")
                       )


samples_list <- lapply(samples_list, normalize_and_scale)

# keep malignants cells only
filter_malignants <- function(sc) {
    subtypes_to_keep <- c("Malignant cells", "tS1", "tS2", "tS3")
    if (sum(sc@meta.data$Cell_subtype %in% subtypes_to_keep) > 0) {
            sc_filtered <- subset(x = sc, subset = Cell_subtype %in% subtypes_to_keep)
            return(sc_filtered)
        } else {
            return(NULL)
        }
    
}

all_malignant <- lapply(samples_list, filter_malignants)


## Get rid of the NULL elements
all_malignant[sapply(all_malignant, is.null)] <- NULL

saveRDS(object = all_malignant,
        file = paste0(data_folder, "/obj/luad_kim_nayoung/all_malignant.rds")
        )
