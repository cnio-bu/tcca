library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

# Load the Beyondcell object with immunotherapy sensitivity predicted
bc <- readRDS("beyondcell_pancancer_immuno.Rds")

# Load cells from subclones with sctherapy predictions
seu_subclones <- readRDS("../single_cell/sctherapy/results/marker_genes/seu_cluster_markers.rds")

# We need based on the "c[digit]" cell id used in the Beyondcell analyses
metadata <- read.table("../single_cell/seurat/tcca/tcca_metadata.tsv",
    sep = "\t",
    header = TRUE
)

metadata <- metadata %>%
    filter(malignancy == "True") %>%
    mutate(cell_id = paste0("c", row_number()))


bc@meta.data <- bc@meta.data %>%
    rownames_to_column(var = "cell_id") %>%
    select(-cell) %>%
    left_join(select(metadata, cell, cell_id), by = "cell_id") %>%
    column_to_rownames(var = "cell_id")
colnames(bc) <- bc@meta.data$cell


# Add beyondcell sensitivity scores to the subclones metadata for dotplot
common_cells <- intersect(colnames(bc), colnames(seu_subclones))
bc <- bc[rownames(bc)[582:nrow(bc)], common_cells]


# Create a Complex heatmap
bc <- bc[rownames(bc)[582:nrow(bc)], common_cells]
bc <- Seurat::SketchData(
    object = bc,
    assay = "RNA",
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch_5k_2"
)

sketch_mat <- bc[["sketch_5k_2"]]$data
sketch_mat <- scale(x = sketch_mat, center = TRUE, scale = TRUE)

cells_annot_df <- as.data.frame(t(sketch_mat)) %>%
    rownames_to_column(var = "cell") %>%
    select(cell) %>%
    left_join(select(metadata, cell, tme_archetype, scTherapy_cluster), by = "cell") %>%
    column_to_rownames(var = "cell")
    
colnames(cells_annot_df) <- c("TME archetype", "Cluster")


pals <- list(
   "TME archetype" = tme_colors,
    "Cluster" = sctherapy_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = cells_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = c("Cluster" = FALSE)
)

tme_legend <- Legend(
    at = unique(pals$`TME archetype`),
    legend_gp = gpar(fill = pals$`TME archetype`),
    ncol = 2, # Split Group legend into 2 columns
    gap = unit(10, "mm"),
    title = "TME archetype",
    title_gp = gpar(fontsize = 12, fontface = "bold"), # Title font
    labels_gp = gpar(fontsize = 12)
)

heat <- ComplexHeatmap::Heatmap(
    mat = sketch_mat,
    top_annotation = top_annotation,
    cluster_rows = TRUE,
    cluster_row_slices = TRUE,
    column_order = rownames(cells_annot_df[order(cells_annot_df$Cluster), ]),
    cluster_columns = FALSE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split = cells_annot_df$Cluster,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_legend_param = list(title = "BCS score"),
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(8, "in")
)

png(
    file = "../single_cell/sctherapy/results/tme_immunotherapy_heatmap.png",
    width = 14,
    height = 14,
    units = "in",
    res = 500
)

ht_opt(
    "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
    "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
    "legend_gap" = unit(1, "cm")
)

draw(heat,
    annotation_legend_side = "top",
    heatmap_legend_side = "bottom",
    annotation_legend_list = list(
    tme_legend
    )
)
dev.off()

scaled_immuno_bc <- round(t(apply(bcs_immuno, 1, scales::rescale, to = c(0, 1))), digits = 2)
scaled_immuno_bc <- as.data.frame(t(scaled_immuno_bc))


# 1. Pivot normalized scores (transposed so cells are rows)
norm_mat <- as.matrix(bc[["RNA"]]$data)
scaled_mat <- t(scale(x = norm_mat, center = TRUE, scale = TRUE))

norm_df <- as.data.frame(t(norm_mat)) %>%
    rownames_to_column("cell") %>%
    pivot_longer(
        cols = -cell,
        names_to = "immunotherapy",
        values_to = "norm_score"
    )

# 2. Pivot scaled scores (already cells as rows)
scaled_df <- as.data.frame(scaled_mat) %>%
    rownames_to_column("cell") %>%
    pivot_longer(
        cols = -cell,
        names_to = "immunotherapy",
        values_to = "scaled_score"
    )

# 3. Join normalized + scaled + clusters
bubble_df <- norm_df %>%
    left_join(scaled_df, by = c("cell", "immunotherapy")) %>%
    left_join(select(metadata, cell, scTherapy_cluster), by = "cell")

# 4. Summarize per cluster and therapy
bubble_df_summary <- bubble_df %>%
    group_by(scTherapy_cluster, immunotherapy) %>%
    summarise(
        pct_positive = mean(norm_score > 0, na.rm = TRUE) * 100,
        mean_scaled_score = mean(scaled_score, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(scTherapy_cluster = factor(scTherapy_cluster, levels = c(1:10)))

# Compute SP per drug
bubble_mean <- ggplot(bubble_df_summary, aes(x = immunotherapy, y = scTherapy_cluster)) +
    geom_point(aes(size = pct_positive, color = mean_scaled_score)) +
     scale_colour_gradientn(colors = c("#1D61F2", "#83A8F7", "#F7F7F7", "#FF9CBB", "#DA0078"), 
                                        limits = c(-0.1, 0.2),
                                        breaks = seq(-0.1, 0.2, by = 0.1),
                                        oob = scales::squish)  +
    scale_size(range = c(1, 6)) +
    theme_minimal() +
    labs(
        color = "Mean Scaled BCS",
        size = "% Positive BCS cells",
        x = "Immunotherapy", y = "Clusters"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12)
    )

ggsave("../single_cell/sctherapy/results/tme_immunotherapy.png",
    plot = bubble_mean,
    width = 10,
    height = 10
)

pdf("sample_wise/figures/bubble_mean_zscaled.pdf",
    width = 15,
    height = 4
)
bubble_mean
dev.off()