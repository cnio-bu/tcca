library("Seurat")
library("tidyverse")

source("src/sc_functions.R")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")


## get all samples
all_samples <- list.dirs(path = paste0(data_folder, "/raw/luad_philip_bisschof"),
                         recursive = FALSE,
                         full.names = FALSE
                         )

full_paths <- paste0(data_folder,
                     "/raw/luad_philip_bisschof/",
                     all_samples
                     )


all_mats <- lapply(X = full_paths, FUN = Seurat::Read10X)

patient_data <- readxl::read_excel(
    path = paste0(
        data_folder,
        "/raw/luad_philip_bisschof/patients_metadata.xlsx")
    ) %>%
    filter(tissue_type == "Tumor")


infer_cnv_scores <- read.csv(
    file = paste0(data_folder, "/raw/luad_philip_bisschof/infercnv_clone_scores_nsclc.tsv"
                  ),
    sep = "\t"
    ) %>%
    filter(tissue_type == "Tumor", !is.na(cell_id), cna_clone == "CNA") %>%
    mutate(
        cell_id = gsub(pattern = "^p.*_", replacement = "", x = cell_id)
    )

names(all_mats) <- all_samples

generate_seurat_objects <- function(dgMat, sample) {
       sample_name <- sample
       seu <- CreateSeuratObject(counts = dgMat,
                                 project = sample,
                                 assay = "RNA",
                                 )
       
       return(seu)
}

all_seurat_objs <- lapply(seq_along(all_mats), function(id) {
    generate_seurat_objects(all_mats[[id]], sample = all_samples[id])
})

rm(all_mats)
gc()

## Filter samples 
all_samples_filtered <- lapply(
    all_seurat_objs,
    filter_sc,
    res_dir = paste0(data_folder, "/qc/luad_philip_bisschof/")
    )

## Normalizing and centering 
all_samples_filtered <- lapply(all_samples_filtered, normalize_and_scale)

## Filter and keep malignants
keep_malignants <- function(sc) {

    sample_name <- unique(sc$orig.ident)
    ## Remove the trailing "t" from names
    sample_name <- gsub(pattern = "t", replacement = "", x = sample_name)
    cells_to_keep <- infer_cnv_scores %>%
        filter(patient_id == sample_name) %>%
        pull(cell_id)
    
    sc_filtered <- subset(sc, cells = cells_to_keep)
    return(sc_filtered)
}

all_malignant <- lapply(all_samples_filtered, keep_malignants)

## Get rid of NULL elements if any
all_malignant[sapply(all_malignant, is.null)] <- NULL

## Annotate with clinical data
fill_clinical <- function(sc) {
    
    this_meta <- patient_data %>%
        filter(sample_id == unique(sc$orig.ident)) %>%
        full_join(y = data.frame(cell_id = colnames(sc)), by = character()) %>%
        as.data.frame()
    
    rownames(this_meta) <- this_meta$cell_id
    this_meta$cell_id <- NULL
    
    sc <- AddMetaData(object = sc, metadata = this_meta)
    return(sc)
    
}

all_malignant_annotated <- lapply(all_malignant, fill_clinical)

saveRDS(object = all_malignant_annotated,
        file = paste0(data_folder,
                      "/obj/luad_philip_bisschof/all_malignant.rds"
                      )
        )


