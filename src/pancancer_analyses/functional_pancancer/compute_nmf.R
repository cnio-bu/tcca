#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
study_name <- args[1]


library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)
library(Matrix)
library(RcppML)

setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load Seurat object with malignant cells
seu <- readRDS("single_cell/seurat/v5/lvl1/seu_lvl1.rds")
seu$study_sample <- paste0(seu$study, ".", seu$sample)

# Subset study
seu_study <- subset(seu, subset = study == study_name)
seu_sample.list <- SplitObject(seu_study, split.by = "study_sample")

# Run NMF per sample in the study
get_nmf_programs <- function(seu_sample, rank, seed = 123) {
    counts <- GetAssayData(seu_sample, assay = "RNA", slot = "counts")
    counts <- as.matrix(counts)
    CP100K_log <- log2(t(t(counts) / colSums(counts)) * 100000 + 1)
    print(dim(CP100K_log))
    CP100K_log <- CP100K_log[apply(
        CP100K_log, 1,
        function(x) length(which(x > 3.5)) > ncol(CP100K_log) * 0.02
    ), ]
    print(dim(CP100K_log))
    CP100K_log <- CP100K_log - rowMeans(CP100K_log)
    CP100K_log[CP100K_log < 0] <- 0
    nmf_sample <- lapply(rank, function(k) {
        nmf <- RcppML::nmf(CP100K_log, k = k, verbose = FALSE, seed = seed)
        print(k)
        rownames(nmf$h) <- paste0("program", 1:nrow(nmf$h))
        colnames(nmf$h) <- colnames(CP100K_log)
        rownames(nmf$w) <- rownames(CP100K_log)
        colnames(nmf$w) <- paste0("program", 1:ncol(nmf$w))
        return(nmf)
    })
    print("NMF computed for a sample")
    names(nmf_sample) <- paste0("k", rank)
    return(nmf_sample)
}

nmf_programs <- lapply(seu_sample.list, function(sample) {
    nmf_sample <- get_nmf_programs(sample, rank = 4:9)
    return(nmf_sample)
})

nmf_programs <- unlist(nmf_programs, recursive = FALSE)
print(paste("NMF computed for each sample in study", study_name))

saveRDS(
    nmf_programs,
    paste0(
        "functional_nmf/sample_wise/nmf_study/",
        study_name,
        "_nmfprograms_samplewise.rds"
    )
)
