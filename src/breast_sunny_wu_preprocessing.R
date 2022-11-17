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
  data.dir = paste0(data_folder, "/raw/breast_sunny_wu"),
  gene.column = 1
)

meta_data <- read.delim(sep=",",row.names = 1,
  file = paste0(data_folder,
                "/raw/breast_sunny_wu/metadata.csv"
  )
)

seu <- Seurat::CreateSeuratObject(counts = all_cells,
                                  project = "breast_sunny_wu",
                                  meta.data = meta_data
)


## Split the merged obj
sample_list <- Seurat::SplitObject(object = seu, split.by = "orig.ident")
names(sample_list) <- sapply(sample_list, function(sc){unique(sc$"orig.ident")})

# Filter cells
samples_filtered <- lapply(
  sample_list,
  filter_sc,
  res_dir = paste0(data_folder, "/qc/breast_sunny_wu")
)

## Normalize and scale data
samples_filtered <- lapply(samples_filtered, normalize_and_scale)

filter_malignant <- function(sc) {
  
  types_to_keep <- "Cancer Epithelial"
  if (sum(sc@meta.data$celltype_major == types_to_keep) > 1) {
    sc_filtered <- subset(x = sc, subset = celltype_major == types_to_keep)
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
  file = paste0(data_folder, "/obj/breast_sunny_wu/all_malignant.rds")
)
