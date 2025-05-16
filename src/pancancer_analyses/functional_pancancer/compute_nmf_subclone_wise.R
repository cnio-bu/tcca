#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
study_name <- args[1]


library(remotes)
library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)
library(Matrix)
library(RcppML)

setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf/subclone_wise")

# Load Seurat object with subclone annotations
seu <- readRDS("seu_lvl1_subclone_annot.rds")

# Subset study
seu_study <- subset(seu, subset = study == study_name)
seu_subclone.list <- SplitObject(seu_study, split.by = "subclone_name")
print("Object splitted into subclones")

# Select subclones with at least 10 cells
seu_subclone.list <- seu_subclone.list[unlist(lapply(seu_subclone.list, ncol) >= 10)]

# Run NMF per subclone in the study
get_nmf_programs <- function(seu_subclone, rank, seed = 123) {
    counts <- seu_subclone[["RNA"]]$counts
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
    if (ncol(CP100K_log) >= 10 && ncol(CP100K_log) < 100) {
        rank <- 2:3
    } else if (ncol(CP100K_log) >= 100 && ncol(CP100K_log) < 300) {
        rank <- 3:5
    } else if (ncol(CP100K_log) >= 300) {
        rank <- 4:9
    }
    nmf_subclone <- lapply(rank, function(k) {
        nmf <- RcppML::nmf(CP100K_log, k = k, verbose = FALSE, seed = seed)
        print(k)
        rownames(nmf$h) <- paste0("program", 1:nrow(nmf$h))
        colnames(nmf$h) <- colnames(CP100K_log)
        rownames(nmf$w) <- rownames(CP100K_log)
        colnames(nmf$w) <- paste0("program", 1:ncol(nmf$w))
        return(nmf)
    })
    print("NMF computed for a subclone")
    names(nmf_subclone) <- paste0("k", rank)
    return(nmf_subclone)
}

nmf_programs <- lapply(seu_subclone.list, function(subclone) {
    nmf_subclone <- get_nmf_programs(subclone, rank = 4:9)
    return(nmf_subclone)
})

nmf_programs <- unlist(nmf_programs, recursive = FALSE)
print(paste("NMF computed for each subclone in study", study_name))

saveRDS(
    nmf_programs,
    paste0(
        "nmf_study_prefilter/",
        study_name,
        "_nmfprograms_subclonewise.rds"
    )
)