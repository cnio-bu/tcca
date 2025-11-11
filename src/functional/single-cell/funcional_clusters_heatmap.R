library(BPCells)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)


setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional/")

## Source TCCA palette
source(file = "/storage/scratch01/users/mgonzalezb/bc-meta/TCCA_palette.R")

sketched_5k <- open_matrix_dir(dir = "results_ucell/sketch_mat_functional_5k/")

sketched_mat <- as.matrix(sketched_5k)
sketched_mat <- t(scale(x = t(sketched_mat), center = TRUE, scale = TRUE))

## human readable origins
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

## functional annotation
seu_ucell <- readRDS("results_ucell/ucell_pancancer.Rds")
fcs <- seu_ucell@meta.data %>%
    filter(
        new_cell_id %in% colnames(sketched_mat)
    )
# fcs <- data.table::fread(
#     "results/fcs.tsv"
# ) %>%
#     filter(
#         new_cell_id %in% colnames(sketched_mat)
#     )


clinical_features <- fcs %>%
    mutate(
        summarised_tumor_site = case_when(
        refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
        TRUE ~ "other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
        #functional_clusters_0.2 = as_factor(functional_clusters_0.2)
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
        # functional_clusters_0.2
        ) %>%
    as.data.frame()

cells_annot_df$summarised_tumor_site <-  translat_human_sites[cells_annot_df$summarised_tumor_site]

rownames(cells_annot_df) <- cells_annot_df$new_cell_id
cells_annot_df$new_cell_id <- NULL

colnames(cells_annot_df) <- c(
    "Sex",
    "Age group",
    "Solid/Liquid",
    "Tumor type",
    "Tumor site",
    "Treatment"
    #"Functional Cluster"
)

pals <- list(
    "Sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Tumor type" = pm_colors,
    "Tumor site" = tumor_sites_colors,
    "Treatment" = treatment_colors
    #"Functional Cluster" = tcs_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df =  cells_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = c("Tumor site" = FALSE))

top_rv <- matrixStats::rowVars(sketched_mat)
top_median <- median(top_rv)
top_rv <- top_rv[top_rv >= top_median] ## median

png(
    file = "results_ucell/sketched_functional_ucell.png",
    res = 300,
    width = 18, 
    height = 18,
    units = "in"
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Tumor site`),
  legend_gp = gpar(fill = pals$`Tumor site`),
  ncol = 2,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Tumor site"
)


heat <- ComplexHeatmap::Heatmap(
    mat = sketched_mat[names(top_rv), ],
    #right_annotation = right_annotation,
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    #column_order = rownames(cells_annot_df[order(cells_annot_df$`Functional Cluster`), ]),
    cluster_row_slices = TRUE,
    # row_split = cells_annot_df$`Therapeutic Cluster`,
    row_title = NULL,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_names_rot = 0,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_legend_param = list(title = "UCell"),
    heatmap_width = unit(8, "in"),
    heatmap_height = unit(14, "in")
)

ht_opt("ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), "HEATMAP_LEGEND_PADDING" = unit(15, "cm"))
draw(heat, annotation_legend_side = "top",  annotation_legend_list = list(tumor_site_legend))
dev.off()
