library(BPCells)
library(Seurat)
library(tidyverse)

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

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

seu <- readRDS("results/lvl1/seu_lvl1_merged.Rds")

gsets <- readGMT(x = "reference/refined_gsets_functional.gmt")

seu <- AddModuleScore(object = seu, features = gsets)
