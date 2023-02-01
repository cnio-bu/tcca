library("Seurat")
library("tidyverse")
## Load common funcs
source("src/sc_functions.R")
cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys
data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")
## Read data
all_cells <- Seurat::Read10X(
  data.dir = paste0(data_folder, "/raw/rcell_kevin_bi"),
  gene.column = 1
)
## Read metadata
meta_data <- read.delim(sep="\t",row.names = 1,
          file = paste0(data_folder,
                        "/raw/rcell_kevin_bi/metadata.txt"
        )
)
## Merge data + metadata
seu <- Seurat::CreateSeuratObject(counts = all_cells,
                                  project = "rcell_kevin_bi",
                                  meta.data = meta_data
)
## Fix missing donor and biosample annotations
seu@meta.data <- seu@meta.data %>%
  mutate(donor_id = toupper(str_remove(rownames(seu@meta.data), pattern = "[A-Z]+\\.")),
         biosample_id = paste0(donor_id, "_scRNA"))
## Switch ident to biosample
seu$old_ident <- seu$orig.ident
seu$orig.ident <- seu$biosample_id
## Split the merged obj
sample_list <- Seurat::SplitObject(object = seu, split.by = "biosample_id")
names(sample_list) <- sapply(sample_list, function(sc){unique(sc$"biosample_id")})
# Filter cells
samples_filtered <- lapply(
  sample_list,
  filter_sc,
  res_dir = paste0(data_folder, "/qc/rcell_kevin_bi")
)
## Normalize and scale data
samples_filtered <- lapply(samples_filtered, normalize_and_scale)