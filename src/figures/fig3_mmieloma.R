library(ComplexHeatmap)
library(circlize)
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
    file = "results/mmieloma/cell_annotation_10k.tsv"
    )

cell_annot_df$new_time <- as.factor(cell_annot_df$new_time)
cell_annot_df$new_time <- fct_relevel(cell_annot_df$new_time, "pre", "post", "post_2")
levels(cell_annot_df$new_time) <- c("Pretreated", "First line", "Second line")

cell_annot_df$sc_gain_1q <- as.factor(cell_annot_df$sc_gain_1q)
levels(cell_annot_df$sc_gain_1q) <- c("WT", "Amplified")

## load matrices
meki_mat <- read.table("results/mmieloma/meki_mat_enrichment.tsv") %>%
    as.matrix()

imid_mat <- read.table("results/mmieloma/imid_mat_enrichment.tsv") %>%
    as.matrix()

pi_mat <- read.table("results/mmieloma/pi_mat_enrichment.tsv") %>%
    as.matrix()


## color palette
timepoints_col <- c(sl_colors, treatment_colors[[1]])
names(timepoints_col) <- levels(cell_annot_df$new_time)

chr1_q_amp_col <- c("#BBB9B7", "#ff4430")
names(chr1_q_amp_col) <- levels(cell_annot_df$sc_gain_1q)

top_col <- list("Timepoint" = timepoints_col, "1q amplification" = chr1_q_amp_col)

draw_and_save_heat <- function(mat){
    
    ## Top annotation
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        "Timepoint" = cell_annot_df[rownames(mat), c("new_time")],
        "1q amplification" = cell_annot_df[rownames(mat), c("sc_gain_1q")],
        col = top_col,
        which = "column"
    )
    
    cell_patient_order <- cell_annot_df[rownames(mat), ]
    cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$PID_new), ])

    b <- ComplexHeatmap::Heatmap(
        name = "Module score",
        mat = scale(t(mat), center = TRUE, scale = TRUE),
        #    col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_order = cell_patient_order,
        cluster_column_slices = TRUE,
        cluster_row_slices = TRUE,
        clustering_distance_columns = "pearson",
        column_split = data.frame(
            cell_annot_df[rownames(mat), ]$drug_t1_response,
            cell_annot_df[rownames(mat), ]$PID_new
        ),
        column_gap = unit(2, "mm"),
        show_column_names = FALSE,
        row_labels = paste0("Metacommunity ", c(1:6)),
        row_names_side = "left",
        top_annotation = top_annotation
    )
    
    return(b)
    
}


## Generate MEKI heat and save 
png(
    filename = "results/figures/mmieloma_meki_heat.png",
    width = 19,
    height = 2.5,
    units = "in",
    res = 100
    )
draw(meki_heat)

dev.off()

## GENERATE PI heat and save
pi_heat <- draw_and_save_heat(pi_mat)

png(
    filename = "results/figures/mmieloma_pi_heat.png",
    width = 19,
    height = 2.5,
    units = "in",
    res = 100
)
draw(pi_heat)

dev.off()

## Generate PIMID heat and save
imid_heat <- draw_and_save_heat(imid_mat)

png(
    filename = "results/figures/mmieloma_imid_heat.png",
    width = 48,
    height = 2.5,
    units = "in",
    res = 100
)
draw(imid_heat)

dev.off()
