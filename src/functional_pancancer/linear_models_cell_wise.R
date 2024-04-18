library(limma)
library(tidyverse)

mp_mat <- readRDS("results/functional/metaprograms_sc_mat.rds")

mp_meta <- data.table::fread("results/functional/cell_level_mps_metadata.tsv")

cell_level_metacoms <- data.table::fread(
    "results/modules/annotated/malignant_cells_best_metacoms_all_cohort.tsv"
)

## deduplicate T19
mp_meta[mp_meta$sample == "T19_1", "sample"] <- "T19"

original_bcodes <- gsub(pattern = "_[^_]+$", replacement = "", x = mp_meta$cell)

## remove trailing .xxxx
test_bcodes_mp_meta <- gsub(
    pattern = "\\.[^.]*$",
    replacement = "",
    x = mp_meta$cell
    )

## remove trailing hashes
test_bcodes_mp_meta <- gsub(
    pattern  = "-[^-]+$", 
    replacement = "",
    x = test_bcodes_mp_meta
)

## remove trailing .xxxx
test_bcodes_metacom <- gsub(
    pattern = "\\.[^.]*$",
    replacement = "",
    x = cell_level_metacoms$cell
    )

## remove trailing hashes
test_bcodes_metacom <- gsub(
    pattern  = "-[^-]+$", 
    replacement = "",
    x = test_bcodes_metacom
)

## remove leading underscores with sample name
test_bcodes_mp_meta <- gsub(
    pattern = ".*_",
    replacement = "",
    x = test_bcodes_mp_meta
)

test_bcodes_metacom <- gsub(
    pattern = ".*_",
    replacement = "",
    x = test_bcodes_metacom
)

## final char removal
test_bcodes_metacom <- gsub(pattern = "\\.", replacement = "", x = test_bcodes_metacom)
test_bcodes_mp_meta <- gsub(pattern = "\\.", replacement = "", x = test_bcodes_mp_meta)

a <- length(intersect(test_bcodes_metacom, test_bcodes_mp_meta))
    
## test
mp_meta$test_cell_id <- test_bcodes_mp_meta
cell_level_metacoms$test_cell_id <- test_bcodes_metacom

test_join <- mp_meta %>%
    inner_join(cell_level_metacoms[, c("test_cell_id", "sample", "study", "new_cell_id")],
               by = c("test_cell_id", "sample", "study")
               )
