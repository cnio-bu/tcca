library(Seurat)
library(BPCells)
library(qlcMatrix)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(qusage)
library(UCell)
setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load TC4 genomic and marker based signatures
tc4_sigs <- read.gmt(
    "single_cell/sctherapy/results/marker_genes/survival_results/tc4_cnvs_new_translated.gmt"
)
cluster_markers <- read.gmt(
    "single_cell/sctherapy/results/marker_genes/marker_sigs_filtered.gmt"
)
tc4_sigs$Cluster04_UP <- cluster_markers$Cluster04_UP

# Load seurat object to compute UCell enrichment scores
seu <- readRDS(
    "single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds"
)
metadata <- read.table(
    "single_cell/seurat/v5/tcca_metadata.tsv",
    header = TRUE,
    sep = "\t",
    row.names = NULL
) %>%
    column_to_rownames("cell") 

seu@meta.data <- metadata[colnames(seu), ]
seu_malignant <- subset(seu, subset = malignancy == "True")
malignant_tc4 <- subset(seu, subset = scTherapy_cluster == 4)

# Compute UCell enrichment scores for TC4 signatures
tc4_ucell <- UCell::AddModuleScore_UCell(malignant_tc4, features = tc4_sigs)

tc4_ucell_mat <- seu@meta.data[, grep("UCell", colnames(seu@meta.data), value = TRUE)]

saveRDS(tc4_ucell, "single_cell/sctherapy/results/esca_validation_bc/tc4_ucell.rds")

# Load beyondcell matrix
mat_bc <- open_matrix_dir(dir = "beyondcell_immuno/full_mat_beyondcell")
colnames(mat_bc) <- colnames(seu)

# Transform to sparse matrix
mat_bc <- as(mat_bc, "sparseMatrix")
mat_tc4 <- as(tc4_ucell_mat, "sparseMatrix")

# Select only ESCA samples from patients
esca_samples <- malignant_tc4@meta.data %>%
    filter(malignancy == "True" & patient != "ccl" & refined_tumor_type == "ESCA") %>%
    pull(old_barcode)

# Extract ESCA samples for correlation
