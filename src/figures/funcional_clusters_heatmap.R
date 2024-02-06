library(BPCells)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)

## Source TCCA palette
source(file = "src/figures/TCCA_palette.R")

sketched_5k <- open_matrix_dir(dir = "results/functional/sketch_mat_functional_5k/")

sketched_mat <- as.matrix(sketched_5k)
sketched_mat <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

## human readable origins
translat_human_sites <- c(
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "adrenal_gland" = "Adrenal gland",
    "breast" = "Breast",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue"
)

## functional annotation
fcs <- data.table::fread(
    "results/annotation/beyondcell_with_therapeutic_clusters.tsv"
) %>%
    filter(
        new_cell_id %in% colnames(sketched_mat)
    )


clinical_features <- fcs %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated, "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        #therapeutic_clusters_0.2 = as_factor(therapeutic_clusters_0.2)
    )



cells_annot_df <- clinical_features %>%
    select(
        new_cell_id,
        sex,
        adult_pediatric,
        is_blood,
        sample_type,
        summarised_tumor_site,
        treated,
        therapeutic_clusters_0.2) %>%
    as.data.frame()

cells_annot_df$summarised_tumor_site <-  translat_human_sites[cells_annot_df$summarised_tumor_site]

rownames(cells_annot_df) <- cells_annot_df$new_cell_id
cells_annot_df$new_cell_id <- NULL

colnames(cells_annot_df) <- c(
    "Chromosomal sex",
    "Age group",
    "Solid/Liquid",
    "Tumor type",
    "Origin",
    "Treatment",
    "Therapeutic Cluster"
)

pals = list(
    "Chromosomal sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Tumor type" = pm_colors,
    "Origin" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "Therapeutic Cluster" = tcs_colors
)

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df =  cells_annot_df,
    which = "row",
    col = pals,
    annotation_name_side = "top",
    annotation_name_rot = 45
)

top_rv <- matrixStats::rowVars(sketched_mat)
top_rv <- top_rv[top_rv >= 0.58] ## median

png(
    file = "results/figures/sketched_functional_with_tcs.png",
    res = 300,
    width = 14, 
    height = 18,
    units = "in"
)


heat <- ComplexHeatmap::Heatmap(
    mat = t(sketched_mat[names(top_rv), ]),
    right_annotation = right_annotation,
    #top_annotation = top_annotation,
    cluster_rows = FALSE,
    row_order = rownames(cells_annot_df[order(cells_annot_df$`Therapeutic Cluster`), ]),
    cluster_row_slices = TRUE,
    row_split = cells_annot_df$`Therapeutic Cluster`,
    row_title = NULL,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split = 4, 
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    #heatmap_width = unit(8, "in"),
    #heatmap_height = unit(14, "in")
)

draw(heat)
dev.off()
