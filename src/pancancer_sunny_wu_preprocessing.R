library("Seurat")
library("tidyverse")

## load common funs
source("src/sc_functions.R")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")

all_cells <- Seurat::Read10X(
    data.dir = paste0(data_folder, "/raw/pancancer_sunny_wu"),
    gene.column = 1
)

meta_data <- read.delim(
    file = paste0(data_folder,
                  "/raw/pancancer_sunny_wu/Wu_etal_2021_metadata.txt"
                  )
    ) %>%
    slice(2:n()) %>%
    as.data.frame()

rownames(meta_data) <- meta_data$NAME
meta_data$NAME <- NULL

all_cells <- Seurat::CreateSeuratObject(counts = all_cells,
                                        project = "pancancer_sunny_wu",
                                        meta.data = meta_data
                                        )


all_cells$old_ident <- all_cells$orig.ident
all_cells$orig.ident <- all_cells$biosample_id
    
sample_list <- Seurat::SplitObject(object = all_cells, split.by = "biosample_id")
names(sample_list) <- sapply(sample_list, function(sc){unique(sc$"biosample_id")})

# Filter cells
samples_filtered <- lapply(
    sample_list,
    filter_sc,
    res_dir = paste0(data_folder, "/qc/pancancer_sunny_wu")
)

## Normalize and scale data
samples_filtered <- lapply(samples_filtered, normalize_and_scale)

filter_malignant <- function(sc) {
    
    types_to_keep <- c("Cancer/Epithelial",
                       "Cancer/ Epithelial Cycling",
                       "Cancer"
                       )
    if (sum(sc@meta.data$CellType %in% types_to_keep) > 1) {
        sc_filtered <- subset(x = sc, subset = CellType %in% types_to_keep)
        return(sc_filtered)
        
    } else {
        return(NULL)
    }
    
}

## Keep malignant cells
malignant_cells <- lapply(samples_filtered, filter_malignant)

## Get rid of the NULL elements/samples with no malignant cells left
malignant_cells[sapply(malignant_cells, is.null)] <- NULL

saveRDS(
    object = malignant_cells,
    file = paste0(data_folder, "/obj/pancancer_sunny_wu/all_malignant.rds")
    )
