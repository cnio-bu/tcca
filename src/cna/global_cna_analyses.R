library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")

cna_mat <- readRDS("raw/cna_level_1_laura/full_genes_copynumber_1layer/full_genes_copynumber.rds")

## Load cell wise metacommunities
cell_wise_metacommunities <- read_tsv(
    "results/modules/annotated/malignant_cells_best_metacoms_all_cohort.tsv"
    )

cna_mat$original_barcode

## test
intersection <- intersect(cna_mat$original_barcode, cell_wise_metacommunities$cell)

cna_mat <- subset(cna_mat, subset = original_barcode %in% intersection)

cell_wise_metacommunities <- cell_wise_metacommunities %>%
    filter(cell %in% intersection)

new_metadata <- cna_mat@meta.data %>%
    left_join(
        y = cell_wise_metacommunities[c("cell", "metacommunity")],
        by = c("original_barcode" = "cell")
        ) %>% 
    as.data.frame()

rownames(new_metadata) <- rownames(cna_mat@meta.data)

cna_mat@meta.data <- new_metadata

Idents(cna_mat) <- cna_mat$metacommunity

## remove cell lines
cna_mat <- subset(cna_mat, subset = study != "cell_lines_gabriella_kinker")

## remove liquid tumours as CNA prediction does not work well there
is_blood <- c("ALL", "CLL", "LAML")

`%nin%` = Negate(`%in%`)

cna_mat <- subset(cna_mat, subset = tumor_type %nin% is_blood)

## TODO UPDATE SEURAT
amps_dels <- FindAllMarkers(
    verbose = TRUE,
    object = cna_mat,
    assay = "RNA",
    slot = "counts", 
    test.use = "t",
    random.seed = 120394,
    return.thresh = 1
    )
