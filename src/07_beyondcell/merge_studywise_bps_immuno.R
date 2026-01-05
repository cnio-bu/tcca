library(BPCells)
library(tidyverse)
library(Seurat)

setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/")

all_mats <- list.dirs(
    path = "studywise_bpcells",
    full.names = TRUE
    )


all_mats <- all_mats[2:length(all_mats)]

mats <- map(
    all_mats,
    open_matrix_dir
)

# All drug set
all_drugs <- unique(unlist(lapply(mats, rownames)))

# Align all matrices so they have the same set of drugs, filling in missing drugs 
# with 0 in studies where they are absent
align_sparse <- function(mat_dir, all_drugs) {
    mat <- as(mat_dir, "dgCMatrix") # extract dgCMatrix
    row_idx <- match(rownames(mat), all_drugs)

    Matrix::sparseMatrix(
        i = rep(row_idx, times = ncol(mat)),
        j = rep(1:ncol(mat), each = nrow(mat)),
        x = as.numeric(mat),
        dims = c(length(all_drugs), ncol(mat)),
        dimnames = list(all_drugs, colnames(mat))
    )
}

# Align all matrices
mats_aligned <- lapply(mats, align_sparse, all_drugs = all_drugs)

# Combine into a single matrix
full_mat <- do.call(cbind, mats_aligned)
full_mat[is.na(full_mat)] <- 0
full_mat <- as(full_mat, "sparseMatrix")

write_matrix_dir(
    mat = full_mat,
    dir = "./full_mat_beyondcell",
    overwrite = TRUE
    )

rm(full_mat, mats)
gc()

mat <- open_matrix_dir(dir = "full_mat_beyondcell")

## load metadata
all.meta <- list.files(
    "studywise_bpcells",
    pattern = "*.tsv",
    full.names = TRUE
    )

meta.data <- all.meta %>%
    map(read.table, row.names = 1, header = TRUE)

meta.data[[27]]$patient <- meta.data[[27]]$orig.ident

for(i in c(1:length(all.meta))){
    this_study <- all.meta[[i]]
    meta.data[[i]]$study <- this_study
}

meta.data_full <-  meta.data  %>%
    map(~.x %>%
             select(
                 nCount_RNA,
                 nFeature_RNA,
                 sample,
                 patient,
                 percent.mt,
                 percent.ribo,
                 S.Score,
                 G2M.Score,
                 Phase,
                 study
                 ) %>%
             mutate(
                 "sample" = as.character(sample),
                 "patient" = as.character(patient),
                 "sample" = str_trim(sample),
                 "patient" = str_trim(patient)
             )) %>%
    bind_rows() %>%
    mutate(
        study = basename(study),
        study = gsub(pattern = "*.tsv", replacement = "", x = study)
    ) %>%
    ## Exploit the fact that we loaded all the mats and meta.data in order
    ## and the indexes will match. Forego  the old barcodes, they are messy
    mutate(
        new_cell_id = c(1:ncol(mat))
    )


## get clinical data
clinical <- data.table::fread(
    "../single_cell/seurat/v5/clinical_metadata_v4_clean_new.tsv"
    )

## Add clinical metadata
meta.data_full_clinical <- meta.data_full %>%
    rownames_to_column("cell") %>%
    inner_join(
        y = clinical,
        by = c("sample" = "sample", "study" = "study")
    ) %>%
    select(-patient.x) %>%
    rename("patient" = patient.y)  %>%
    mutate(
        refined_tumor_site = case_when(
            refined_tumor_site == "" ~ "Unknown",
            TRUE ~ refined_tumor_site
        )
    ) %>%
    mutate(
        original_cell_id = gsub(
            pattern = "\\.\\.\\..*$", ## annoying ... by seurat
            replacement = "",
            x = cell
        )
    )

colnames(mat) <- c(1:ncol(mat))
mat2 <- mat[, meta.data_full_clinical$new_cell_id]

write_matrix_dir(
    mat = mat2,
    dir = "full_mat_beyondcell",
    overwrite = TRUE
    )

write_tsv(
    x = meta.data_full_clinical,
    file = "beyondcell_metadata_with_clinical.tsv"
    )
