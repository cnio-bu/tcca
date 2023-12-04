library(BPCells)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)

## Source TCCA palette
source(file = "src/figures/TCCA_palette.R")

sketched_mat <- open_matrix_dir(dir = "results/beyondcell_bp/sketch_mat_beyondcell")
sketched_5k <- open_matrix_dir(dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k/")

sketched_mat <- as.matrix(sketched_5k)
sketched_mat <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

## therapeutic clusters
tcs <- data.table::fread(
    "results/annotation/beyondcell_with_therapeutic_clusters.tsv"
    ) %>%
    filter(
        new_cell_id %in% colnames(sketched_mat)
    )

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

clinical_features <- tcs %>%
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
        therapeutic_clusters_0.2 = as_factor(therapeutic_clusters_0.2)
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
top_rv <- top_rv[top_rv >= 2]


## get drug names
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
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

drugs <- drugs[drugs$IDs %in% names(top_rv), ]
rownames(drugs) <- drugs$IDs

MoAs <- drugs[names(top_rv), c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"

moa_pals <- list(
    "Mechanism of action" = MoAs_colors
)


top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = MoAs,
    which = "column",
    col = moa_pals,
    show_annotation_name = FALSE
)

png(
    file = "results/figures/sketched_beyondcell_with_tcs.png",
    res = 300,
    width = 14, 
    height = 18,
    units = "in"
    )


test <- ComplexHeatmap::Heatmap(
    mat = t(sketched_mat[names(top_rv),]),
    #mat = t(sketched_mat),
    right_annotation = right_annotation,
    top_annotation = top_annotation,
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
    column_labels = drugs[names(top_rv), "preferred.drug.names"],
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_width = unit(8, "in"),
    heatmap_height = unit(14, "in")
)

draw(test)
dev.off()

