library(BPCells)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(igraph)
library(patchwork)
library(tidyverse)

## load color pal
source(file = "src/figures/TCCA_palette.R")

## load drug data
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    distinct() %>%
    as.data.frame()

rownames(drugs) <- drugs$IDs

## load annotations
cell_annot_df <- read.table(
    file = "results/paad/seu_metadata_malignants.tsv"
) 

## load sketched mat for heatmaps
sketched_mat <- read.table("results/paad/metacom_mat_2k_sketch.tsv") %>%
    as.matrix()

cell_annot_sketch <- cell_annot_df[rownames(sketched_mat), ]

metacom_human_names <- c(
    paste(
        "Metacommunity",
        rep(c("untreated"), times = 6),
        c(1:6)
    )
)

cell_patient_order <- cell_annot_sketch[rownames(sketched_mat), ]
cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$patient), ])


##continuous heat pal
## load gemcitabine sensitivity
bc_mat <- readRDS("results/paad/bc_normalized_mat_rgrss.rds")
gemcitabine_gs <- bc_mat["sig-20902", ]

cell_annot_sketch$gemcitabine <- gemcitabine_gs[rownames(cell_annot_sketch)]
cell_annot_sketch$gemcitabine <- scale(
    cell_annot_sketch$gemcitabine,
    center = TRUE,
    scale = TRUE
    )

top_col <- list(
    "Primary or metastasis" = pm_colors,
    "Therapeutic cluster" = 
    "Gemcitabine - GDSC" = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
 #   "Bortezomib - CTRP" = colorRamp2(c(-6, 0, 6), c("blue", "white", "red")),
 #   "Bortezomib - PRISM" = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
)

cell_annot_sketch$sample_type <- as.factor(cell_annot_sketch$sample_type)
levels(cell_annot_sketch$sample_type) <- c("Metastasis", "Primary") 

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Primary or metastasis" = cell_annot_sketch[rownames(sketched_mat), ]$sample_type,
    "Gemcitabine - GDSC" = cell_annot_sketch[rownames(sketched_mat), ]$gemcitabine,
    col = top_col,
    which = "column",
    show_annotation_name = TRUE,
    annotation_legend_param = list(
    "Gemcitabine - GDSC" = list(at = c(-4, -2, 0, 2, 4))
    )
)



sketched_mat_heat <- ComplexHeatmap::Heatmap(
    name = "Module score",
    mat = t(scale(sketched_mat, center = TRUE, scale = TRUE)),
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    cluster_row_slices = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    #column_order = cell_patient_order,
    column_split = data.frame(
        cell_annot_sketch[rownames(sketched_mat), ]$patient,
        cell_annot_sketch[rownames(sketched_mat), ]$tumor_site
    ),
    cluster_column_slices = TRUE,
    clustering_distance_columns = "pearson",
    row_labels = metacom_human_names,
    row_names_gp = gpar(fontsize = 8),
    show_row_names = TRUE,
    column_names_side = "bottom",
    top_annotation = top_annotation,
    column_title_gp = gpar(fontsize = 8),
    column_title_rot = 45
)

png(
    filename = "results/figures/paad_primary_met_clustered.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(sketched_mat_heat)

dev.off()

sketched_mat_heat_sorted <- ComplexHeatmap::Heatmap(
    name = "Module score",
    mat = t(scale(sketched_mat, center = TRUE, scale = TRUE)),
    cluster_rows = TRUE,
    clustering_distance_rows = "pearson",
    cluster_row_slices = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_order = cell_patient_order,
    column_split = data.frame(
        cell_annot_sketch[rownames(sketched_mat), ]$patient,
        cell_annot_sketch[rownames(sketched_mat), ]$tumor_site
    ),
    cluster_column_slices = TRUE,
    clustering_distance_columns = "pearson",
    row_labels = metacom_human_names,
    row_names_gp = gpar(fontsize = 8),
    show_row_names = TRUE,
    column_names_side = "bottom",
    top_annotation = top_annotation,
    column_title_gp = gpar(fontsize = 8),
)

png(
    filename = "results/figures/paad_primary_met_sorted.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(sketched_mat_heat_sorted)

dev.off()
