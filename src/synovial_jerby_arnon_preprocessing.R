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


## Load Jerby Arnon data, preprocessed by M.J.Jimenez
all_synovial <- readRDS(paste0(data_folder,
                               "/raw/synovial_jerby_arnon/seurat_pre-qc.rds")
                        )

## This is a seurat obj. with an unknown v. Converting it to Seurat 4...
mat <- GetAssayData(all_synovial, slot = "counts")
meta.data <- all_synovial@meta.data

synovial <- Seurat::CreateSeuratObject(counts = mat, project = "synovial",
                                       meta.data = meta.data)

rm(mat, meta.data, all_synovial)
gc()

# Change default ident for later split and QC
synovial@meta.data$orig.ident <- synovial@meta.data$characteristics..sample

## Split the merged obj
synovial.samples <- Seurat::SplitObject(
    object = synovial,
    split.by = "characteristics..sample"
    )

names(synovial.samples) <- unique(synovial@meta.data$characteristics..sample)

synovial_filtered <- lapply(X = synovial.samples,
                            FUN = filter_sc,    
                            res_dir = paste0(data_folder, "/qc/synovial_jerby_arnon")
                            )


synovial_filtered <- lapply(X = synovial_filtered,
                            FUN = normalize_and_scale
)

# Get rid of "NULL" samples, there are no cells left in these
synovial_filtered[sapply(synovial_filtered, is.null)] <- NULL

## Filter and subset malignant cells
keep_malignant <- function(sc) {
    sc_filtered <- subset(
        x = sc,
        subset = `characteristics..cell.type` == "Malignant"
        )
    return(sc_filtered)
}

synovial_malignant <- lapply(synovial_filtered, keep_malignant)    

# Get rid of "NULL" samples, if there were no malignants left
synovial_malignant[sapply(synovial_malignant, is.null)] <- NULL

saveRDS(
    object = synovial_malignant, 
    file = paste0(data_folder, "obj/synovial_jerby_arnon/all_malignant.rds")
)
    
