library(BPCells)
library(ComplexHeatmap)
library(uwot)
library(tidyverse)
library(Seurat)

all_mats <- list.dirs(
    path = "results/beyondcell_bp",
    full.names = TRUE
    )


# yay!
cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

all_mats <- all_mats[2:36]

mats <- map(
    all_mats,
    open_matrix_dir
)

full_mat <- do.call(cbind.fill, mats)
full_mat[is.na(full_mat)] <- 0
full_mat <- as(full_mat, "sparseMatrix")

write_matrix_dir(
    mat = full_mat,
    dir = "results/beyondcell_bp/full_mat_beyondcell"
    )

rm(full_mat, mats)
gc()

mat <- open_matrix_dir(dir = "results/beyondcell_bp/full_mat_beyondcell/")

## load metadata
all.meta <- list.files(
    "results/beyondcell_bp/",
    pattern = "*.tsv",
    full.names = TRUE
    )


meta.data <- all.meta %>%
    map(read.table, row.names = 1, header = TRUE) 

meta.data[[26]]$patient <- meta.data[[26]]$orig.ident

for(i in c(1:length(all.meta))){
    print(i)
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
             )) %>%
    bind_rows() %>%
    mutate(
        study = basename(study),
        study = gsub(pattern="*.tsv", replacement="", x = study)
    )


## get clinical data
clinical <- data.table::fread("results/annotation/clinical_metadata_v2_clean.tsv")
clinical[
    clinical$sample == "T19" &
        clinical$study == "adrenalnb_rui_chong",
    "sample"
] <- "T19_1"

meta.data_full <- meta.data_full %>%
    mutate(
        sample = case_when(
            study == "adrenal_nb_rui_chong" & sample == "T19" ~ "T19_1",
            TRUE ~ sample
        )
    )

## Add clinical metadata
meta.data_full_clinical <- meta.data_full %>%
    rownames_to_column("cell") %>%
    left_join(
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
        ),
        new_cell_id = c(1:ncol(mat))
    )

write_tsv(
    x = meta.data_full_clinical, 
    file = "results/annotation/beyondcell_metadata_with_clinical.tsv"
    )
