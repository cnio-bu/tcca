#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
study_name <- args[1]


library(remotes)
library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)
library(Matrix)

setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf/subclone_wise")

# Load Seurat object with subclone annotations
seu <- readRDS("seu_lvl1_subclone_annot.rds")

# Subset study
seu_study <- subset(seu, subset = study == study_name)

# Split to compute normalization per sample
seu_study[["RNA"]] <- split(seu_study[["RNA"]], f = seu_study$sample)
seu_study <- Seurat::NormalizeData(seu_study,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)
seu_study <- Seurat::FindVariableFeatures(seu_study,
    selection.method = "vst",
    nfeatures = 7000
)
hvg <- VariableFeatures(seu_study)

# Create a list of Seurat objects (one for each subclone)
seu_study <- JoinLayers(seu_study)
print("Layers joined")
seu_subclone.list <- SplitObject(seu_study, split.by = "subclone_name")
print("Object splitted into subclones")


# Compute NMF for each subclone in the study
# geneNMF.programs <- GeneNMF::multiNMF(
#     seu_subclone.list,
#     # Use consensus variable features identified across batches
#     hvg = VariableFeatures(seu_study),
#     assay = "RNA",
#     slot = "data",
#     k = 2:9,
#     nfeatures = 7000
# )

seu_subclone.list <- seu_subclone.list[unlist(lapply(seu_subclone.list, ncol) >= 10)]
nmf_programs <- lapply(seu_subclone.list, function(subclone) {
    mat <- subclone[["RNA"]]$data[hvg, ]
    mat <- as(mat, "dgCMatrix")
    if (ncol(subclone) >= 10 && ncol(subclone) < 100) {
        k_values <- 2:3
    } else if (ncol(subclone) >= 100 && ncol(subclone) < 300) {
        k_values <- 3:5
    } else if (ncol(subclone) >= 300) {
        k_values <- 4:9
    }
    nmf_subclone <- lapply(k_values, function(k) {
        nmf <- RcppML::nmf(mat, k = k, verbose = FALSE, seed = 123)
        print(k)
        rownames(nmf$h) <- paste0("program", 1:nrow(nmf$h))
        colnames(nmf$h) <- colnames(mat)
        rownames(nmf$w) <- rownames(mat)
        colnames(nmf$w) <- paste0("program", 1:ncol(nmf$w))
        return(nmf)
    })
    print("NMF computed for a subclone")
    names(nmf_subclone) <- paste0("k", k_values)
    return(nmf_subclone)
})

nmf_programs <- unlist(nmf_programs, recursive = FALSE)
print("NMF computed for each subclone")

saveRDS(
    nmf_programs,
    paste0(
        "nmf_study/",
        study_name,
        "_geneNMFprograms_allsubclones.rds")
)