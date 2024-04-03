library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")

cna_mat <- readRDS("raw/cna_level_1_laura/full_genes_copynumber_1layer/full_genes_copynumber.rds")
cna_mat <- UpdateSeuratObject(cna_mat)

## Load cell wise metacommunities
cell_wise_metacommunities <- read_tsv(
    "results/modules/annotated/malignant_cells_best_metacoms_all_cohort.tsv"
    )


intersection <- intersect(cna_mat$original_barcode, cell_wise_metacommunities$cell)

cell_wise_metacommunities <- cell_wise_metacommunities %>%
    filter(cell %in% intersection)

## remove liquid tumours as CNA prediction does not work well there
is_blood <- c("ALL", "CLL", "LAML")
`%nin%` = Negate(`%in%`)

new_metadata <- cna_mat@meta.data %>%
    rownames_to_column("row_names") %>%
    filter(
        original_barcode %in% intersection,
        study != "cell_lines_gabriella_kinker",
        tumor_type %nin% is_blood
    ) %>%
    left_join(
        y = cell_wise_metacommunities[c("cell", "metacommunity")],
        by = c("original_barcode" = "cell")
        ) %>% 
    as.data.frame()

rownames(new_metadata) <- new_metadata$row_names

new_mat <- BPCells::open_matrix_dir(dir = "results/cna_mat_lvl1_gene/")
new_mat <- new_mat[, new_metadata$row_names]

seu <- CreateSeuratObject(counts = new_mat, meta.data = new_metadata)

Idents(seu) <- seu$metacommunity

## Done in new Seurat 5 due to bugs in the dev. version of v5
seu@assays[["RNA"]]$data <- new_mat

amps_dels <- FindAllMarkers(
    verbose = TRUE,
    object = seu,
    assay = "RNA",
    slot = "counts", 
    test.use = "LR",
    random.seed = 120394,
    logfc.threshold = 0,
    return.thresh = 1
    )

## bp test
amps_dels <- BPCells::marker_features(mat = new_mat, groups = new_metadata$metacommunity)
