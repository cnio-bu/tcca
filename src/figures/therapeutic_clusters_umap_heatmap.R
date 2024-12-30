library(BPCells)
library(ComplexHeatmap)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta/")

## Source TCCA palette
source(file = "~/figures/TCCA_palette.R")

sketched_mat <- open_matrix_dir(dir = "beyodncell/results/sketch_mat_beyondcell")
sketched_5k <- open_matrix_dir(dir = "beyodncell/results/sketch_mat_beyondcell_5k/")

sketched_mat <- as.matrix(sketched_5k)
sketched_mat <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

## therapeutic clusters
tcs <- data.table::fread(
    "results/tcs.tsv"
) %>%
    filter(
        new_cell_id %in% colnames(sketched_mat)
    )

## extract sex inferred
seu <- readRDS("./single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
seu <- subset(seu, subset = malignancy == TRUE)
colnames(seu) <- paste0("c", c(1:ncol(seu)))

tcs$sex <- seu@meta.data[tcs$new_cell_id, "sex"]

tcs <- tcs %>%
    filter(
        new_cell_id %in% colnames(sketched_mat)
    )

## human readable origins
translat_human_sites2 <- c(
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "adrenal_gland" = "Adrenal gland",
    "breast" = "Breast",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue",
    "bladder" = "Bladder",
    "colon" = "Colon",
    "kidney" = "Kidney"
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
        therapeutic_clusters_k.300.res.0.5 = as_factor(therapeutic_clusters_k.300.res.0.5)
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
        therapeutic_clusters_k.300.res.0.5) %>%
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

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df =  cells_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = c("Sample site" = FALSE)
)

top_rv <- matrixStats::rowVars(sketched_mat)
top_rv <- top_rv[top_rv >= 2]

## get drug names
drugs <- data.table::fread("~/reference/final_moas - Collapsed.tsv") %>%
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


right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = MoAs,
    which = "row",
    col = moa_pals,
    show_annotation_name = FALSE
)

png(
    file = "results/sketched_beyondcell_with_tcs.png",
    res = 300,
    width = 19, 
    height = 14,
    units = "in"
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  ncol = 2,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Sample site"
)

test <- ComplexHeatmap::Heatmap(
    mat = sketched_mat[names(top_rv),],
    #mat = t(sketched_mat),
    right_annotation = right_annotation,
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    cluster_row_slices = TRUE,
    row_split = 5,
    row_title = NULL,
    column_order = rownames(cells_annot_df[order(cells_annot_df$`Therapeutic Cluster`), ]),
    cluster_columns = FALSE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split =  cells_annot_df$`Therapeutic Cluster`, 
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    row_labels = drugs[names(top_rv), "preferred.drug.names"],
    show_row_names = TRUE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_legend_param = list(title = "BCS score"),
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(8, "in")
)
ht_opt("ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), 
       "HEATMAP_LEGEND_PADDING" = unit(1, "cm"), 
       "legend_gap" = unit(1, "cm"))
draw(heat, 
     annotation_legend_side = "top",  
     annotation_legend_list = list(tumor_site_legend))
draw(test)
dev.off()


## Generate UMAP plot from raw UMAP
tcs_umap <- readRDS("results/tcs_umap.rds")

tcs_umap_clean <- tcs_umap +
    ggtitle("") +
    scale_color_manual(
        name = "Therapeutic cluster",
        values = tcs_colors
        ) +
    scale_shape_manual() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )

tcs_umap_clean$layers[[1]]$aes_params$size <- 0.2
tcs_umap_clean$layers[[1]]$aes_params$alpha <- 0.9

ggsave(
    plot = tcs_umap_clean,
    filename = "results/umap_tcs.png",
    dpi = 500,
    height = 7,
    width = 7
    )
