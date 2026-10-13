library(dplyr)
library(tidyverse)
library(UCell)
library(Matrix)
library(matrixStats)
library(data.table)
library(Seurat)

# Sample name as input
args <- commandArgs(trailingOnly = TRUE)
s <- args[1]
if (is.na(s) || s == "") stop("Sample name not indicated")

setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load data
seu_hq_malignant <- readRDS("single_cell/seurat/tcca/seu_hq_malignant.rds")
expr <- GetAssayData(seu_hq_malignant, assay = "RNA", layer = "data")

cell_meta <- seu_hq_malignant@meta.data %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(cell, study_sample)

# Load list of meta-programs
meta_programs <- readRDS("functional_nmf/sample_wise/metaprograms_cpm/mp_list_reordered.rds")

# Define parameters for the analysis
ncores <- 4
maxrank <- 1500
B <- 50
n_bins <- 25
min_genes <- 30
frac_min <- 0.40
min_cells <- 100

out.dir <- "functional_nmf/therapeutic_from_mps/mp_high_cells"

message("[", Sys.time(), "] Processing sample: ", s)
fout <- file.path(out.dir, paste0("mp_", make.names(s), ".rds"))
if (file.exists(fout)) {
  message("It already exists: ", fout)
  quit(save = "no")
}


############# Statistical method to robustly select high-MP cells ##############
# 1. Function to generate control sets based on gene expression bins
make_control_sets <- function(mp_genes, gene_bin, B) {
    bins_needed <- table(gene_bin[mp_genes])
    pools <- lapply(names(bins_needed), function(bn) {
        setdiff(names(gene_bin)[gene_bin == as.integer(bn)], mp_genes)
    })
    names(pools) <- names(bins_needed)

    lapply(seq_len(B), function(b) {
        unlist(lapply(names(bins_needed), function(bn) {
            n_need <- bins_needed[[bn]]
            pool <- pools[[bn]]
            if (length(pool) == 0) {
                return(NULL)
            }
            if (length(pool) < n_need) {
                sample(pool, n_need, replace = TRUE)
            } else {
                sample(pool, n_need)
            }
        }), use.names = FALSE)
    })
}

# Function to get UCell scores for MP and control gene sets per sample
process_sample <- function(s) {

    cells_s <- cell_meta[cell_meta$study_sample == s, "cell"]
    if (length(cells_s) < min_cells) return(NULL)

    expr_s <- expr[, cells_s, drop = FALSE]
    expr_s <- as(expr_s, "dgCMatrix")

    # Detected genes in the sample
    keep <- Matrix::rowSums(expr_s > 0) > 0
    expr_s <- expr_s[keep, , drop = FALSE]

    if (nrow(expr_s) < 2000) return(NULL)

    # Specific binning in that sample
    gene_mean <- Matrix::rowMeans(expr_s)
    gene_bin <- cut(rank(gene_mean, ties.method = "random"),
                    breaks = n_bins, labels = FALSE)
    names(gene_bin) <- rownames(expr_s)

    # Ranking: once per sample
    ranks_s <- StoreRankings_UCell(expr_s, maxRank = maxrank, ncores = ncores)

    res_list <- list()

    for (mp_name in names(meta_programs)) {

        mp_genes <- intersect(meta_programs[[mp_name]], rownames(expr_s))
        if (length(mp_genes) < min_genes) next

        controls <- make_control_sets(mp_genes, gene_bin, B)
        controls <- Filter(function(x) length(x) >= min_genes, controls)

        if (length(controls) < 20) next
        names(controls) <- sprintf("ctrl%03d", seq_along(controls))

        sigs <- c(list(observed = mp_genes), controls)
        scores <- ScoreSignatures_UCell(precalc.ranks = ranks_s,
                                        features = sigs,
                                        ncores = ncores)

        obs <- scores[, "observed_UCell"]
        null_mat <- scores[, grep("^ctrl", colnames(scores)), drop = FALSE]

        null_mu <- rowMeans(null_mat)
        null_sd <- rowSds(null_mat)
        null_sd[null_sd < 1e-8] <- 1e-8

        z <- (obs - null_mu) / null_sd

        r_sub <- ranks_s[mp_genes, , drop = FALSE]
        frac_detected <- Matrix::colMeans(r_sub > 0 & r_sub <= maxrank)

        res_list[[mp_name]] <- data.table(
            cell_name = rownames(scores),
            sample = s,
            mp = mp_name,
            ucell = obs,
            z = z,
            pval = pnorm(z, lower.tail = FALSE),
            frac_detected = frac_detected
        )
    }

    rm(ranks_s, expr_s); gc()
    if (length(res_list) == 0) return(NULL)
    rbindlist(res_list)
}


# Run the functions
res <- process_sample(s)
if (!is.null(res)) saveRDS(res, fout)

message("[", Sys.time(), "] Completed: ", s)