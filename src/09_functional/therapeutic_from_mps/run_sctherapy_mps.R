library(dplyr)
library(tidyverse)
library(Seurat)
library(data.table)
library(httr)
library(logger)
library(jsonlite)
library(readr)

setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf")

# Load predict_compounds function from scTherapy
source("/home/mgonzalezb/bc-meta/scTherapy/predict_compounds.R")

# Receive the sample name as argument
args <- commandArgs(trailingOnly = TRUE)
s <- args[1]
if (is.na(s) || s == "") stop("Sample name not provided")
message("[", Sys.time(), "] Processing sample: ", s)

# Load gene info
gene_list <- "/home/mgonzalezb/bc-meta/scTherapy/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()

# Load MP-high cells
mp_high_cells <- readRDS("therapeutic_from_mps/mp_high_cells/mp_high_cells.rds")
mp_groups <- readRDS("therapeutic_from_mps/mp_high_cells/mp_groups_viable.rds")

# Filter MP groups for this sample
mp_groups_s <- mp_groups[mp_groups$sample == s, ]

if (nrow(mp_groups_s) == 0) {
    message("No viable MP groups for sample: ", s)
    quit(save = "no")
}

message("MPs to process: ", nrow(mp_groups_s))

# Function to compute DEGs
get_mp_deg <- function(
  seu_sample,
  mp_high_cell_names,
  padj_threshold = 0.05,
  lfc_threshold = 1,
  stable_lfc = 0.1
) {
    n_healthy <- sum(seu_sample$malignancy == "False")
    n_mp <- length(mp_high_cell_names)

    seu_sample$sctherapy_ident <- case_when(
        colnames(seu_sample) %in% mp_high_cell_names ~ "mp_high",
        seu_sample$malignancy == "False" ~ "healthy",
        TRUE ~ NA_character_
    )

    cells_use <- colnames(seu_sample)[!is.na(seu_sample$sctherapy_ident)]
    seu_use <- seu_sample[, cells_use]
    Idents(seu_use) <- "sctherapy_ident"

    deg <- tryCatch(
        {
            FindMarkers(
                object = seu_use,
                ident.1 = "mp_high",
                ident.2 = "healthy",
                min.pct = -Inf,
                logfc.threshold = -Inf,
                min.cells.feature = 1,
                min.cells.group = 1,
                test.use = "wilcox"
            )
        },
        error = function(e) {
            message("FindMarkers error: ", e$message)
            return(NULL)
        }
    )

    if (is.null(deg) || nrow(deg) == 0) {
        return(NULL)
    }

    deg_filtered <- deg %>%
        mutate(gene_symbol = rownames(.)) %>%
        inner_join(gene_info, by = "gene_symbol") %>%
        filter(
            (p_val_adj <= padj_threshold & abs(avg_log2FC) > lfc_threshold) |
                abs(avg_log2FC) < stable_lfc
        )

    if (nrow(deg_filtered) == 0) {
        return(NULL)
    }

    setNames(as.list(deg_filtered$avg_log2FC), deg_filtered$gene_symbol)
}

# Load Seurat and metadata
seu_lvl2 <- readRDS("../single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
metadata_hq <- read.table(
    "../single_cell/seurat/tcca/tcca_metadata_h5ad_sctherapy_cells.tsv",
    header = TRUE, sep = "\t"
) %>% mutate(study_sample = paste0(study, "_", sample))

# Subset Seurat for this sample and cells
meta_s <- metadata_hq %>%
    filter(study_sample == s) %>%
    column_to_rownames(var = "cell")
cells_s <- rownames(meta_s)

if (length(cells_s) == 0) {
    message("No cells found for sample: ", s)
    quit(save = "no")
}

seu_s <- subset(seu_lvl2, cells = cells_s)
seu_s <- NormalizeData(seu_s)
seu_s@meta.data <- meta_s[colnames(seu_s), ]

n_nonmal <- sum(seu_s$malignancy == "False")
# if (n_nonmal < 20) {
#   message("Insufficient non-malignant cells: ", n_nonmal)
#   quit(save = "no")
# }

message(
    "Total cells: ", ncol(seu_s),
    "| Non-malignant: ", n_nonmal,
    "| Malignant: ", sum(seu_s$malignancy == "True")
)


# Setup output
out_dir <- "therapeutic_from_mps/sctherapy_mps/"

# Process each MP for this sample
fout_sample <- paste0(out_dir, make.names(s), "_mps_drugs.tsv")
if (file.exists(fout_sample)) {
    message("[", Sys.time(), "] Already exists: ", s)
    quit(save = "no")
}

sample_tables <- list()

for (i in seq_len(nrow(mp_groups_s))) {
    mp_name <- mp_groups_s$mp[i]
    message("  [", mp_name, "] Processing...")

    # Get MP-high cells for this MP in this sample
    mp_high_names <- mp_high_cells %>%
        filter(mp_high_cells$sample == s, mp_high_cells$mp == mp_name) %>%
        pull("cell_name")

    if (length(mp_high_names) < 20) {
        message(
            "[", mp_name, "] Skipping - insufficient MP-high cells: ",
            length(mp_high_names)
        )
        next
    }

    message("[", mp_name, "] MP-high cells: ", length(mp_high_names))

    # Compute DEGs
    deg_list <- get_mp_deg(
        seu_sample = seu_s,
        mp_high_cell_names = mp_high_names
    )

    if (is.null(deg_list)) {
        message("[", mp_name, "] No DEGs found")
        next
    }

    # Predict drugs
    drugs <- tryCatch(
        {
            predict_drugs(
                degs_list = deg_list,
                exclude_low_confidence = TRUE
            )
        },
        error = function(e) {
            message(paste(Sys.time(), "Error in", s, mp_name, ":", e$message, "\n"))
            return(NULL)
        }
    )

    if (is.null(drugs) || nrow(drugs) == 0) {
        message("[", mp_name, "] No drugs predicted")
        next
    }

    drugs$sample <- s
    drugs$mp <- mp_name
    sample_tables <- c(sample_tables, list(drugs))
    message("[", mp_name, "] Done-drugs: ", nrow(drugs))
}

# Save the sample table with all the MPs
if (length(sample_tables) > 0) {
    sample_drug_table <- do.call(rbind, sample_tables)
    write.table(sample_drug_table, fout_sample,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    message(
        "[", Sys.time(), "] Completed: ", s,
        " | MPs: ", length(sample_tables),
        " | Drugs: ", nrow(sample_drug_table)
    )
} else {
    message("[", Sys.time(), "] No results for: ", s)
}
