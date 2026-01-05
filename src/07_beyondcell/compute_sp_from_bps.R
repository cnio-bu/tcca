library(BPCells)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)

setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/")

##### FUNCTIONS ####
# Funciton to compute the Switch Point from BCS
compute_sp <- function(norm_bcs, scaled_bcs) {
    if (any(norm_bcs == 0)) {
        sp <- rep(which(norm_bcs == 0)[1], times = 2)
    } else {
        lower.bound <- which(norm_bcs == max(norm_bcs[norm_bcs <= 0]))[1]
        upper.bound <- which(norm_bcs == min(norm_bcs[norm_bcs >= 0]))[1]
        sp <- c(lower.bound, upper.bound)
    }
    sp_scaled <- round(sum(scaled_bcs[sp]) / 2, digits = 2)
    return(sp_scaled)
}

bpcell_mat <- list.dirs(
    path = "./studywise_bpcells",
    full.names = TRUE
)
bpcell_mat <- bpcell_mat[-1]
metadata_tsv <- list.files(
    path = "./studywise_bpcells",
    pattern = ".tsv",
    full.names = TRUE
)

# Total drugs
adrenalnb_rui_chong <- open_matrix_dir("studywise_bpcells/adrenalnb_rui_chong")
all_drugs <- rownames(adrenalnb_rui_chong)

sp_list <- list()
for (study in bpcell_mat) {
    study_name <- basename(study)
    bc_mat <- open_matrix_dir(study)
    bc_metadata <- read.table(paste0(study, ".tsv"))
    if (study_name == "pdac_junya_peng") {
        bc_metadata$patient <- bc_metadata$orig.ident
    }
    bc_mat <- as.matrix(bc_mat)
    bc_mat[is.na(bc_mat)] <- 0

    # Align all matrices so they have the same set of drugs, filling in missing drugs
    # with 0 in studies where they are absent
    row_idx <- match(rownames(bc_mat), all_drugs)

    mat_complete <- Matrix::sparseMatrix(
        i = rep(row_idx, times = ncol(bc_mat)),
        j = rep(1:ncol(bc_mat), each = nrow(bc_mat)),
        x = as.numeric(bc_mat),
        dims = c(length(all_drugs), ncol(bc_mat)),
        dimnames = list(all_drugs, colnames(bc_mat))
    )
    
    bc_mat <- mat_complete

    bc_metadata <- bc_metadata %>%
        mutate(study_sample = paste0(study_name, "_", sample))
    bc_mat <- as.data.frame(t(bc_mat))
    bc_mat$study_sample <- bc_metadata$study_sample
    
    # Compute mean BCS and residuals' mean for each sample
    # bc_mat <- bc_mat %>%
    #     pivot_longer(cols = starts_with("sig-"), names_to = "drug", values_to = "BCS") %>%
    #     group_by(study_sample, drug) %>%
    #     # mutate(BCS_scaled = scales::rescale(BCS, to = c(0, 1))) %>%
    #     summarise(
    #         mean_BCS = mean(BCS, na.rm = TRUE),
    #         mean_resid = mean(BCS - mean_BCS, na.rm = TRUE)
    #     ) %>%
    #     select(study_sample, drug, mean_resid) %>%
    #     pivot_wider(names_from = study_sample, values_from = mean_resid, values_fill = list(mean_resid = NA)) %>%
    #     column_to_rownames(var = "drug")
    # bc_mat <- as.matrix(bc_mat)

    # Compute SP for each sample
    bc_mat <- bc_mat %>%
        pivot_longer(cols = -study_sample, names_to = "drug", values_to = "BCS") %>%
        group_by(study_sample, drug) %>%
        mutate(BCS_scaled = scales::rescale(BCS, to = c(0, 1)))
    sp_scaled <- bc_mat %>%
        group_by(study_sample, drug) %>%
        summarise(
            sp_scaled = compute_sp(BCS, BCS_scaled),
            .groups = "drop"
        ) %>%
        pivot_wider(names_from = study_sample, values_from = sp_scaled, values_fill = list(sp_scaled = NA)) %>%
        column_to_rownames(var = "drug")
    
    sp_scaled <- as.matrix(sp_scaled)
    sp_list <- append(sp_list, list(study_name = sp_scaled))
    # rm(bc_mat, sp_scaled)
    rm(bc_mat)
    gc()
    print(study_name)
}

saveRDS(sp_list, "sp_drug_per_sample.rds")




source(file = "~/bc-meta/figures/TCCA_palette.R")
sp_list <- readRDS("sp_drug_per_sample.rds")
sp_mat <- do.call(cbind, sp_list)
colnames(sp_mat) <- gsub(" ", "", basename(colnames(sp_mat)))
# Plot drugs that do not have NA for any of the samples.
sp_mat <- sp_mat[apply(sp_mat, 1, function(x) !any(is.na(x))), ]
metadata <- read.table("beyondcell_metadata_with_clinical.tsv", sep = "\t", header = TRUE)
metadata <- metadata %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(
        study_sample,
        study, 
        sample, 
        patient, 
        sex, 
        age, 
        tumor_type, 
        refined_tumor_site, 
        sample_type, 
        treated
    ) %>%
    distinct()

seu <- readRDS("../single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
sex <- seu@meta.data %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(study_sample, sex) %>%
    distinct()

metadata <- metadata %>%
    left_join(sex, by = "study_sample")

translat_human_sites <- c(
    "adrenal_gland" = "Adrenal gland",
    "bladder" = "Bladder",
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "breast" = "Breast",
    "colon" = "Colon",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "kidney" = "Kidney",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue"
)


clinical_features <- metadata %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex.y == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
    )


sample_annot_df <- clinical_features %>%
    select(
        study_sample,
        sex,
        adult_pediatric,
        is_blood,
        sample_type,
        summarised_tumor_site,
        treated,
        study
    ) %>%
    as.data.frame()

sample_annot_df$summarised_tumor_site <- translat_human_sites[sample_annot_df$summarised_tumor_site]

rownames(sample_annot_df) <- gsub(" ", "", sample_annot_df$study_sample)
sample_annot_df$study_sample <- NULL
sp_mat <- sp_mat[, rownames(sample_annot_df)]
colnames(sample_annot_df) <- c(
    "Sex",
    "Age group",
    "Solid/Liquid",
    "Sample type",
    "Sample site",
    "Treatment",
    "Study"
)
study_colors <- scales::hue_pal()(36)
names(study_colors) <- unique(metadata$study)
pals <- list(
    "Sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Sample type" = pm_colors,
    "Sample site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "Study" = study_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = sample_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = c("Sample site" = FALSE, "Study" = FALSE)
)


## get drug names
drugs <- data.table::fread("~/bc-meta/reference/final_moas - Collapsed.tsv") %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    mutate(
        collapsed.MoAs = case_when(
            preferred.drug.names == "VANDETANIB" ~ "VEGFR inhibitor",
            preferred.drug.names == "DASATINIB" ~ "Kinase inhibitor",
            preferred.drug.names == "RIGOSERTIB" ~ "Other",
            preferred.drug.names == "SORAFENIB" ~ "Kinase inhibitor",
            TRUE ~ collapsed.MoAs
        )
    ) %>%
    distinct() %>%
    as.data.frame()

drugs <- drugs[drugs$IDs %in% rownames(sp_mat), ]
rownames(drugs) <- drugs$IDs

MoAs <- drugs[rownames(sp_mat), c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"


png(
    file = "results/heatmap_mean_bcs_per_sample.png",
    res = 500,
    width = 19,
    height = 14,
    units = "in"
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
    at = names(pals$`Sample site`),
    legend_gp = gpar(fill = pals$`Sample site`),
    ncol = 2, # Split Group legend into 2 columns
    gap = unit(10, "mm"),
    title = "Sample site"
)
study_legend <- Legend(
    at = names(pals$`Study`),
    legend_gp = gpar(fill = pals$`Study`),
    ncol = 4, # Split Group legend into 2 columns
    gap = unit(10, "mm"),
    title = "Study"
)


heat <- ComplexHeatmap::Heatmap(
    mat = sp_mat,
    #col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    top_annotation = top_annotation,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = drugs[rownames(sp_mat), "collapsed.MoAs"],
    column_order = rownames(sample_annot_df)[order(sample_annot_df$`Study`)],
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    column_split = sample_annot_df$`Study`,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    row_title = unique(drugs[rownames(sp_mat), "collapsed.MoAs"]),
    row_order = rownames(sp_mat)[order(drugs[rownames(sp_mat), "collapsed.MoAs"])],
    show_row_names = FALSE,
    column_names_side = "top",
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 8),
    column_title = NULL,
    heatmap_legend_param = list(title = "Switch Point"),
    heatmap_width = unit(17, "in"),
    heatmap_height = unit(11, "in")
)

ht_opt(
    "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
    "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top", annotation_legend_list = list(tumor_site_legend, study_legend))
dev.off()



