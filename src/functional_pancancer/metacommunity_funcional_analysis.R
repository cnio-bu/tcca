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

functional_mat <- as.matrix(functional_mat)

met_to_tests <- metacoms[, c("cell_id", "metacom_untreated_6")]
rownames(met_to_tests) <- met_to_tests$cell_id
met_to_tests$cell_id <- NULL

design <- model.matrix(~metacom_untreated_6, data = met_to_tests)

fit <- lmFit(functional_mat, design)
fit <- eBayes(fit)
fit2 <- topTable(
    fit = fit,
    coef = "metacom_untreated_6",
    number = Inf,
    adjust.method = "fdr"
    )

fit3_mt6 <- fit2[abs(fit2$logFC) >= 1, ]

all_fits <- list(fit3_mt1, fit3_mt2, fit3_mt3, fit3_mt4, fit3_mt5, fit3_mt6)
all_fits <- all_fits %>%
    map(as.data.frame) %>%
    map(rownames_to_column) %>%
    bind_rows(.id = "metacommunity")

write.table(
    x = all_fits,
    file = "results/annotation/metacommunity_functional.tsv"
    )

bastard_mat <- t(bastard_mat)
test <- cor(bastard_mat, method = "pearson")
test3 <- test[c(1:130), c(131:142)]
library(ComplexHeatmap)
rownames(test3) <- tolower(rownames(test3))
test2 <- ComplexHeatmap::Heatmap(
    t(test3),
    column_split = 2,
    column_names_gp = gpar((fontsize = 4)),
    column_names_rot = 90,
    row_split = 2
    )

png(
    filename = "results/figures/metaprograms_functional_correlation.png",
    width = 19,
    height = 6,
    units = "in",
    res = 200
    )

draw(test2)
dev.off()
