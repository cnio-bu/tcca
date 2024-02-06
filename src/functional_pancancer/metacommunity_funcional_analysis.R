library(BPCells)
library(limma)
library(Seurat)
library(tidyverse)

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

functional_mat <- open_matrix_dir(
    dir = "results/functional/full_mat_functional/"
    )

meta.data <- read_tsv(
    "results/annotation/functional_metadata_with_clinical.tsv"
    )

meta.data_full_clinical <- meta.data %>%
    as.data.frame()

rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

## remove empty rows
functional_mat <- functional_mat[1:130, ]

fc <- CreateSeuratObject(
    counts = functional_mat,
    assay = "RNA",
    project = "fc_pancancer",
    meta.data = meta.data_full_clinical
)

## load metacom mat
mcs <- readRDS("results/beyondcell_bp/beyondcell_pancancer.Rds")

## generate mats for limma
metacoms <- mcs@meta.data %>%
    select(metacom_untreated_1:metacom_treated_6) %>%
    rownames_to_column("cell_id")

functional_mat <- as.matrix(functional_mat[1:5, ])

met_to_tests <- metacoms[, c("cell_id", "metacom_untreated_1")]
rownames(met_to_tests) <- met_to_tests$cell_id
met_to_tests$cell_id <- NULL

design <- model.matrix(~metacom_untreated_1, data = met_to_tests)

fit <- lmFit(functional_mat, design)
fit <- eBayes(fit)
fit2 <- topTable(
    fit = fit,
    coef = "metacom_untreated_1",
    number = Inf,
    adjust.method = "fdr"
    )

functional_dt <- as.data.frame(t(functional_mat))
functional_dt$metacom1 <- metacoms$metacom_untreated_1

test <- ggplot(data = functional_dt, aes(x = scale(metacom1), y = scale(`HALLMARK-MITOTIC-SPINDLE`))) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "loess") + 
    ggpubr::stat_cor()
