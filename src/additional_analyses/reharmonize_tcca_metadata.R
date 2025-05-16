library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/")

metadata <- read.table("tcca_annotation_raw.tsv", header = TRUE, sep = "\t")

# Check whether sex column includes the inferred sex for unknown samples
seu.lvl2 <- readRDS("../v5/lvl2/seu_lvl2_sex_inferred.rds")
identical(metadata$cell, rownames(seu.lvl2@meta.data))

metadata$sex <- seu.lvl2@meta.data$sex

# Add subclonal info to all cells
clonality <- read.table("../../cna_metadata/full_clonality_table_lvl2.tsv", row.names = NULL)

duplicated <- clonality %>%
    filter(original_barcode %in% original_barcode[duplicated(original_barcode)])

clonality <- clonality %>%
    filter(!(original_barcode %in% duplicated$original_barcode)) %>% # Remove 16 cells with duplicated names
    dplyr::select(original_barcode, scevan_subclone)

metadata <- left_join(metadata, clonality, by = c("cell" = "original_barcode")) %>%
    mutate(
        scevan_subclone = ifelse(scevan_subclone != "non_tumor",
            paste(study, sample, scevan_subclone, sep = "."), scevan_subclone
        ),
        malignancy = ifelse(malignancy == "True", TRUE, FALSE)
    ) %>%
    relocate(scevan_subclone, .after = "tme_archetype")

clinical <- read.table("../v5/clinical_metadata_v4_new.tsv", header = TRUE, sep = "\t")

# An error detected with BRCA sample 1 biological replicate in Hugo Gonzalez et al.
# and there was a mismatch with some ages from this study
subset <- clinical[clinical$study == "brmets_hugo_gonzalez", ]
metadata_hugo <- metadata[metadata$study == "brmets_hugo_gonzalez", ] %>%
    left_join(subset, by = "sample", suffix = c("", ".new"))
cols_to_replace <- c(
    "patient", "sex", "age", "tumor_type", "tumor_subtype", "stage", "tnm", "tumor_site",
    "sample_type", "treated", "treatment_type", "treatment_response", "treatment_info",
    "os", "pfi", "sequencing_tech", "genome_assembly", "data_pmid", "refined_tumor_site"
)
metadata_hugo <- metadata_hugo %>%
    mutate(across(
        all_of(cols_to_replace),
        ~ coalesce(get(paste0(cur_column(), ".new")), .)
    )) %>%
    select(-ends_with(".new"))

metadata[metadata$study == "brmets_hugo_gonzalez", ] <- metadata_hugo

## Rename tumor type according to TCGA cancer types
# All low grade gliomas (oligodendroglioma and astrocytoma) are classified as GBMs
# in the current version. We must update them to OGD and A in the tumor type column
metadata <- metadata %>%
    mutate(tumor_type = ifelse(grepl("astrocitoma|astrocytoma", tumor_subtype), "A", tumor_type))


cancer_types <- list(
    "GBM" = c("GBM", "MB"),
    "LGG" = c("A", "OGD"),
    "NB" = c("NB", "GNB"),
    "LAML" = c("LAML"),
    "ALL" = c("ALL"),
    "CLL" = c("CLL"),
    "MM" = c("MM"),
    "SKCM" = c("SKCM", "SKAM"),
    "SKSC" = c("SKSC"),
    "BCC" = c("BCC"),
    "UVM" = c("UVM"),
    "SARC" = c("SARC", "GIST", "PLEU"),
    "MESO" = c("MESO"),
    "BRCA" = c("BRCA"),
    "LUAD" = c("LUAD", "NSCLC"),
    "LUSC" = c("LUSC"),
    "LCLC" = c("LCLC"),
    "SCLC" = c("SCLC"),
    "OV" = c("OV"),
    "COAD" = c("COAD"),
    "READ" = c("READ"),
    "CESC" = c("CESC"),
    "UCEC" = c("UCEC"),
    "UCS" = c("UCS"),
    "LIHC" = c("LIHC"),
    "CHOL" = c("CHOL"),
    "BLCA" = c("BLCA", "MISC", "KTCC"),
    "HNSC" = c("HNSC"),
    "PRAD" = c("PRAD"),
    "KIRC" = c("KIRC", "KRCC", "KIRCH"),
    "ESCA" = c("ESCA", "ESCC"),
    "PAAD" = c("PAAD"),
    "THCA" = c("THCA"),
    "STAD" = c("STAD")
)


cancer_types <- enframe(cancer_types, name = "refined_tumor_type", value = "tumor_type") %>%
    unnest()

metadata <- metadata %>%
    left_join(cancer_types, by = "tumor_type") %>%
    relocate(refined_tumor_type, .after = "tumor_type")

clinical <- clinical %>%
    left_join(cancer_types, by = "tumor_type") %>%
    relocate(refined_tumor_type, .after = "tumor_type")

write.table(metadata, "tcca_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(clinical, "../v5/clinical_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
