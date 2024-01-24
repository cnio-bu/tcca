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
    file = "results/aml/bc_metadata_annotated.tsv"
) %>%
    select(cell_id:nFeature_RNA)

cell_annot_df$Subgroup <- as.factor(cell_annot_df$Subgroup)
cell_annot_df$Patient_ID <- as.factor(cell_annot_df$Patient_ID)
cell_annot_df$Known_CNVs <- as.factor(cell_annot_df$Known_CNVs)
cell_annot_df$Treatment_Outcome <- as.factor(cell_annot_df$Treatment_Outcome)
cell_annot_df$Patient_Sample <- as.factor(cell_annot_df$Patient_Sample)

## load sketched mat for heatmaps
sketched_mat <- read.table("results/aml/metacom_mat_20k_sketch.tsv") %>%
    as.matrix()

cell_annot_sketch <- cell_annot_df[rownames(sketched_mat), ]

sketched_subset <- sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "RUNX", "cell_id"], 
]

## test
## Top annotation
top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Timepoint" = cell_annot_sketch[rownames(sketched_subset), ]$Patient_Sample,
    "Treatment outcome" = cell_annot_sketch[rownames(sketched_subset), ]$Treatment_Outcome,
   # col = top_col,
    which = "column"
)



cell_patient_order <- cell_annot_sketch[rownames(sketched_subset), ]
cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$Patient_ID), ])


b <- ComplexHeatmap::Heatmap(
    name = "Module score",
    mat = t(sketched_subset),
    #    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    column_order = cell_patient_order,
    column_split = data.frame(
        cell_annot_sketch[rownames(sketched_subset), ]$Patient_Sample,
        cell_annot_sketch[rownames(sketched_subset), ]$Patient_ID
    ),
    top_annotation = top_annotation
)
