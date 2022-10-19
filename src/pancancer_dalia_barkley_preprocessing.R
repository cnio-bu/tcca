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

## Load the source file... big RData obj
load(paste0(data_folder, "/raw/pancancer_dalia_barkley/srt.list.primary.all.RData"))

## Change the default assay to RNA for all samples
change_to_rna <- function(sc) {
    DefaultAssay(sc) <- "RNA"
    return(sc)
}

srt.list.primary.all <- lapply(srt.list.primary.all, change_to_rna)

## Filter cells
filtered_sc <- lapply(
    srt.list.primary.all,
    filter_sc,
    res_dir = paste0(data_folder, "/qc/pancancer_dalia_barkley")
    )

rm(srt.list.primary.all)
gc()

filtered_sc <- lapply(filtered_sc, normalize_and_scale)

# Get rid of "NULL" samples, there are no cells left in these
filtered_sc[sapply(filtered_sc, is.null)] <- NULL

## Discard samples without cells and keep malignants only
keep_all_malignants <- function(sc) {
    
    # This dataset might be empty before even filtering
    if (!is.null(sc)) {
        if (sum(sc@meta.data$type == "malignant") > 0) {
            sc_filtered <- subset(x = sc, subset = type == "malignant")
            return(sc_filtered)
        } else {
            return(NULL)
        }
    }
    
}

filtered_malignant <- lapply(filtered_sc, keep_all_malignants)

## Get rid of the NULL elements
filtered_malignant[sapply(filtered_malignant, is.null)] <- NULL

saveRDS(
    object = filtered_malignant,
    file = paste0(data_folder, "/obj/pancancer_dalia_barkley/all_malignant.rds")
    )
