library(BPCells)
library(GSEABase)
library(Seurat)
library(tidyverse)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

readGMT <- function(x) {
    # --- Checks ---
    # Check x.
    if (length(x) != 1 | !is.character(x[1])) {
        stop('x must be a single string.')
    }
    # Check that x exists.
    if(!file.exists(x)) stop(paste0(x, ' does not exist.'))
    # Check that x is a GMT file.
    is.gmt <- stringr::str_detect(basename(x), pattern = "\\.gmt$")
    if (!is.gmt) stop(paste0(x, ' must be a GMT file.'))
    # --- Code ---
    # Read GMT.
    vector.gmt <- readLines(x)
    # Create list.
    list.gmt <- lapply(vector.gmt, FUN = function(x) {
        unique(unlist(strsplit(x, split = "\t")))
    })
    # List of gene sets (without name and description).
    gmt <- lapply(list.gmt, FUN = function(x) {
        x <- x[-c(1:2)]
        unique(x[which(x != "")])
    })
    # Add names.
    names(gmt) <- lapply(list.gmt, FUN = function(x) x[1])
    return(gmt)
}


seu <- readRDS("results/tcca/tcca_seurat_raw.rds")
mat <- BPCells::open_matrix_dir(dir = "results/tcca/raw_matrix_tcca")

seu@assays$RNA$counts <- mat

## filter non malignant cells
seu <- subset(seu, subset = malignancy == TRUE)

## make sure to filter out non malignant brca
seu <- subset(
    seu,
    subset = (tumor_type == "BRCA" & tumor_subtype != "predicted_tumour" | tumor_type != "BRCA")
)

## remove T10 sample, missmatches with bc indexes and we need them to match
seu <- subset(seu, subset = sample != "T10")
seu$new_cell_id <- paste0("c", 1:ncol(seu))

seu <- NormalizeData(seu)
    
## module ES
c3hallmarks <- readGMT(x = "reference/c4.3ca.v2023.2.Hs.symbols.gmt")
    
seu <- AddModuleScore(
    object = seu,
    features = c3hallmarks,
    seed = 120394,
    name = "metaprograms"
)

meta.programs <- seu@meta.data %>%
    select(metaprograms1:metaprograms41) %>%
    as.matrix()


saveRDS(
    object = meta.programs,
    file = "results/functional/metaprograms_sc_mat.rds"
    )

## export the meta.data too for matching
sc_mat_meta.data <- seu@meta.data

sc_mat_meta.data <- sc_mat_meta.data %>%
    rownames_to_column("original_rownames")

write_tsv(
    sc_mat_meta.data,
    file = "results/functional/cell_level_mps_metadata.tsv"
    )
