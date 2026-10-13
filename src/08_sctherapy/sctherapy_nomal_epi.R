library(dplyr)
library(tidyverse)
library(Seurat)
library(httr)
library(logger)
library(jsonlite)
library(readr)
setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load predict_compounds function from scTherapy
source("/home/mgonzalezb/bc-meta/scTherapy/predict_compounds.R", local = FALSE)

# Function to compute the list of DEG per subclone
gene_list <- "/home/mgonzalezb/bc-meta/scTherapy/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()
get_subclone_deg <- function(
    study_sample, 
    min_healthy_cells = 20,
    min_malignant_cells = 10,
    padj_threshold = 0.05,
    lfc_threshold = 1,
    stable_lfc = 0.1){
    

    # Verify enough number of cells
    n_healthy <- sum(!study_sample$malignancy == "True")
    n_malignant <- sum(study_sample$malignancy == "True")

    if (n_healthy < min_healthy_cells) {
        message("Insufficient healthy cells: ", n_healthy, " < ", min_healthy_cells)
        return(NULL)
    }
  
    if (n_malignant < min_malignant_cells) {
        message("Insufficient malignant cells: ", n_malignant, " < ", min_malignant_cells)
        return(NULL)
    }

    # Asign identities to cells
    study_sample$sctherapy_ident <- ifelse(
        study_sample$malignancy == "False", "healthy", study_sample$scevan_subclone
    )
    Idents(study_sample) <- "sctherapy_ident"

    # Get subclones
    all_subclones <- study_sample@meta.data %>%
        filter(malignancy == "True") %>%
        pull(scevan_subclone) %>%
        unique()

    # Compute DEGs
    deg_all_subclones <- list()
    for (subclone in all_subclones){
        deg_subclone <- FindMarkers(
            object = study_sample,
            ident.1 = subclone,
            ident.2 = "healthy",
            min.pct = -Inf,
            logfc.threshold = -Inf,
            min.cells.feature = 1,
            min.cells.group = 1,
            test.use = "wilcox"
        )
        deg_subclone$gene_symbol <- rownames(deg_subclone)

        # Filter informative genes for scTherapy
        # Significantly altered genes (up or down)
        deg_subclone <- deg_subclone %>%
            mutate(gene_symbol = rownames(.)) %>% 
            inner_join(gene_info, by = "gene_symbol") %>%
            filter((p_val_adj <= padj_threshold & (abs(avg_log2FC) > lfc_threshold)) | (abs(avg_log2FC) < stable_lfc))

        deg_subclone <- setNames(as.list(deg_subclone$avg_log2FC), deg_subclone$gene_symbol)
        deg_all_subclones <- c(deg_all_subclones, list(deg_subclone))
    }

    names(deg_all_subclones) <- all_subclones
    return(deg_all_subclones)
}


# Load scTherapy cells (high-quality cells only)
seu_lvl2 <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
seu_lvl2 <- NormalizeData(seu_lvl2)

metadata_hq <- read.table(
    "single_cell/seurat/tcca/tcca_metadata_h5ad_sctherapy_cells.tsv",
    header = TRUE,
    sep = "\t"
) %>%
    mutate(study_sample = paste0(study, "_", sample))

samples_with_epithelial <- metadata_hq %>%
    filter(malignancy == "False" & cell_type_broad == "Epithelial") %>%
    group_by(study_sample) %>%
    summarise(n_epithelial = n()) %>%
    filter(n_epithelial >= 20) %>%
    pull(study_sample)

# Subset seurat objetct and metadata
new_metadata <- metadata_hq %>%
    filter(
        study_sample %in% samples_with_epithelial & 
        (malignancy == "True" | (malignancy == "False" & cell_type_broad == "Epithelial"))
        )

subset_sctherapy <- subset(seu_lvl2, cells = new_metadata$cell)
new_metadata <- new_metadata %>% column_to_rownames(var = "cell")
subset_sctherapy@meta.data <- new_metadata
subset_sctherapy <- JoinLayers(subset_sctherapy)

# Run scTherapy for each sample
out_dir <- "single_cell/sctherapy/nonmal_epithelial/"
error_log_file <- paste0(out_dir, "error_log.txt")
all_sample_tables <- list()

for (s in unique(new_metadata$study_sample)) {

    fout_deg <- paste0(out_dir, make.names(s), "_deg.rds")
    fout_drugs <- paste0(out_dir, make.names(s), "_drugs.tsv")

    # Checkpointing
    if (file.exists(fout_drugs)) {
        message("[", Sys.time(), "] Already exists: ", s)
    next
    }

    message("[", Sys.time(), "] Processing: ", s)

    # Select sample cells
    cells_s <- new_metadata %>%
        filter(study_sample == s) %>%
        rownames_to_column(var = "cell") %>%
        pull(cell)

    if (length(cells_s) == 0) next

    seu_s <- subset(subset_sctherapy, cells = cells_s)

    # Verify that sample have both cell types
    n_malignant  <- sum(seu_s$malignancy == "True")
    n_epithelial <- sum(seu_s$malignancy == "False" & seu_s$cell_type_broad == "Epithelial")

    if (n_malignant == 0 || n_epithelial < 20) {
        message("Skipping - malignant: ", n_malignant, " epithelial: ", n_epithelial)
    next
    }

    # Compute DEGs
    deg_all_subclones <- get_subclone_deg(seu_s)

    if (is.null(deg_all_subclones)) {
        message("No DEGs found for: ", s)
    next
    }

    saveRDS(deg_all_subclones, fout_deg)

    # Predict drugs per subclone
    drugs_per_subclone <- setNames(
        lapply(names(deg_all_subclones), function(subclone_name) {
        tryCatch({
            predict_drugs(
            degs_list = deg_all_subclones[[subclone_name]],
            exclude_low_confidence = TRUE
            )
        }, error = function(e) {
            message(paste(Sys.time(), "Error in subclone", subclone_name, ":", e$message))
            return(NULL)
        })
        }),
        names(deg_all_subclones)
    )

    # Remove NULL results
    drugs_per_subclone <- Filter(Negate(is.null), drugs_per_subclone)

    if (length(drugs_per_subclone) == 0) {
        message("No drugs predicted for: ", s)
        next
    }

    # Build drug table with sample and subclone info
    drugs_table <- do.call(rbind, lapply(names(drugs_per_subclone), function(subclone_name) {
        df <- drugs_per_subclone[[subclone_name]]
        df$sample <- s
        df$subclone <- subclone_name
        return(df)
    }))

    # Save sample results
    write.table(drugs_table, fout_drugs, sep = "\t", quote = FALSE, row.names = FALSE)

    all_sample_tables <- c(all_sample_tables, list(drugs_table))
    message("Done - subclones: ", length(drugs_per_subclone),
            " drugs: ", nrow(drugs_table))
}

# Build and save full table
if (length(all_sample_tables) > 0) {
    full_table <- do.call(rbind, all_sample_tables)
    write.table(
    full_table,
    paste0(out_dir, "full_table_drug_prediction_epithelial.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
    )
    message("Total samples processed: ", length(all_sample_tables))
    message("Total subclones: ", length(unique(full_table$subclone)))
    message("Total drug predictions: ", nrow(full_table))
} else {
    message("No results generated")
}
