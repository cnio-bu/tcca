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

## load bortezomib and sorafenib susceptibility data
bc_seu <- open_matrix_dir("results/aml/aml_beyondcell_mat")

## sig-21115 = SORAFENIB - CTRP
## sig-21377 = BORTEZOMIB  - CTRP
## sig-21378 = BORTEZOMIB - PRISM

drugs_to_keep <- bc_seu[c("sig-21115", "sig-21377", "sig-21378"), ]
drugs_to_keep <- as.matrix(drugs_to_keep)
drugs_annot <- as.data.frame(t(drugs_to_keep))
colnames(drugs_annot) <- c("sorafenib", "bortezomib_ctrp", "bortezomib_prism")

drugs_annot$sorafenib <- scale(
    x = drugs_annot$sorafenib,
    center = TRUE,
    scale = TRUE
    )

drugs_annot$bortezomib_ctrp <- scale(
    x = drugs_annot$bortezomib_ctrp,
    center = TRUE,
    scale = TRUE
    )

drugs_annot$bortezomib_prism <- scale(
    x = drugs_annot$bortezomib_prism,
    center = TRUE,
    scale = TRUE
    )

draw_and_save_heat <- function(mat){
    
    ## Top annotation
    ## color palette
    timepoints_col <- c(sl_colors, treatment_colors[[1]])
    names(timepoints_col) <- levels(cell_annot_sketch[rownames(mat), ]$Patient_Sample)
    
    ##continuous heat pal
    
    top_col <- list(
        "Timepoint" = timepoints_col,
        "Sorafenib - CTRP" = colorRamp2(c(-6, 0, 6), c("blue", "white", "red")),
        "Bortezomib - CTRP" = colorRamp2(c(-6, 0, 6), c("blue", "white", "red")),
        "Bortezomib - PRISM" = colorRamp2(c(-6, 0, 6), c("blue", "white", "red"))
        )
    
    
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        "Timepoint" = cell_annot_sketch[rownames(mat), ]$Patient_Sample,
        "Treatment outcome" = cell_annot_sketch[rownames(mat), ]$Treatment_Outcome,
        "Sorafenib - CTRP" =  drugs_annot[rownames(mat), ]$sorafenib,
        "Bortezomib - CTRP" =  drugs_annot[rownames(mat), ]$bortezomib_ctrp,
        "Bortezomib - PRISM" =  drugs_annot[rownames(mat), ]$bortezomib_prism,
        col = top_col,
        which = "column",
        show_annotation_name = TRUE,
        annotation_legend_param = list("Sorafenib - CTRP" = list(
            at = c(-6, 0, 6)
            ),
            "Bortezomib - CTRP" = list(at = c(-6, 0, 6)),
            "Bortezomib - PRISM" = list(at = c(-6, 0, 6))
        )
    )
    
    
    
    cell_patient_order <- cell_annot_sketch[rownames(mat), ]
    cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$Patient_ID), ])
    
    metacom_human_names <- c(
        paste(
            "Metacommunity",
            rep(c("untreated"), times = 6),
            c(1:6)
        ),
        paste(
            "Metacommunity",
            rep(c("treated"), times = 6), 
            c(1:6)
        )
    )
    
    b <- ComplexHeatmap::Heatmap(
        name = "Module score",
        mat = t(scale(mat, center = TRUE, scale = TRUE)),
        # col = col_fun,
        cluster_rows = TRUE,
        clustering_distance_rows = "pearson",
        cluster_row_slices = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_order = cell_patient_order,
        column_split = data.frame(
            cell_annot_sketch[rownames(mat), ]$Patient_ID,
            cell_annot_sketch[rownames(mat), ]$Patient_Sample
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
    
    
}

## Generate RUNX heat and save 
runx_heat <- draw_and_save_heat(sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "RUNX", "cell_id"], 
]
)
## Generate FLT heat and save 
flt_heat <- draw_and_save_heat(sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "FLT", "cell_id"], 
]
)
## Generate CBFB heat and save 
cbfb_heat <- draw_and_save_heat(sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "CBFB", "cell_id"], 
]
)
## Generate KMT2A heat and save 
kmt_heat <- draw_and_save_heat(sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "KMT2A", "cell_id"], 
]
)
## Generate Other heat and save 
other_abs_heat <- draw_and_save_heat(sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Subgroup == "Other", "cell_id"], 
]
)

png(
    filename = "results/figures/aml_runx_metacom_heat.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(runx_heat)

dev.off()

png(
    filename = "results/figures/aml_flt_metacom_heat.png",
    width = 20,
    height = 4,
    units = "in",
    res = 200
)
draw(flt_heat)

dev.off()



png(
    filename = "results/figures/aml_cbfb_metacom_heat.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(cbfb_heat)

dev.off()

png(
    filename = "results/figures/aml_kmt_metacom_heat.png",
    width = 24,
    height = 8,
    units = "in",
    res = 200
)
draw(kmt_heat)

dev.off()


png(
    filename = "results/figures/aml_other_metacom_heat.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(other_abs_heat)

dev.off()

## test; get only diagnostic samples
diax_mat <- sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Patient_Sample == "Diagnosis", "cell_id"], 1:6]

## test; get only relapse samples
relapse_mat <- sketched_mat[cell_annot_sketch[
    cell_annot_sketch$Patient_Sample == "Relapse", "cell_id"], 1:6]


draw_and_save_heat_cond2 <- function(mat){
    
    ## Top annotation
    ## color palette
    driver_col <- MoAs_colors[1:length(unique(cell_annot_sketch$Subgroup))]
    names(driver_col) <- levels(cell_annot_sketch[rownames(mat), ]$Subgroup)
    
    top_col = list(
        "Driver subgroup" = driver_col,
        "Sorafenib - CTRP" = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
        "Bortezomib - CTRP" = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
        "Bortezomib - PRISM" = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
        )

    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
        "Driver subgroup" = cell_annot_sketch[rownames(mat), ]$Subgroup,
        "Treatment outcome" = cell_annot_sketch[rownames(mat), ]$Treatment_Outcome,
        "Sorafenib - CTRP" =  drugs_annot[rownames(mat), ]$sorafenib,
        "Bortezomib - CTRP" =  drugs_annot[rownames(mat), ]$bortezomib_ctrp,
        "Bortezomib - PRISM" =  drugs_annot[rownames(mat), ]$bortezomib_prism,
        col = top_col,
        which = "column",
        show_annotation_name = TRUE,
        annotation_legend_param = list("Sorafenib - CTRP" = list(
            at = c(-4, 0, 4)
        ),
        "Bortezomib - CTRP" = list(at = c(-4, 0, 4)),
        "Bortezomib - PRISM" = list(at = c(-4, 0, 4))
        )
    )
    
    cell_patient_order <- cell_annot_sketch[rownames(mat), ]
    cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$Patient_ID), ])
    
    metacom_human_names <- paste(
            "Metacommunity",
            rep(c("untreated"), times = 6),
            c(1:6)
        )
    
    
    b <- ComplexHeatmap::Heatmap(
        name = "Module score",
        mat = t(scale(mat, center = TRUE, scale = TRUE)),
        # col = col_fun,
        cluster_rows = TRUE,
        clustering_distance_rows = "pearson",
        cluster_row_slices = TRUE,
        cluster_columns = FALSE,
        cluster_column_slices = TRUE,
        clustering_distance_columns = "pearson",
        show_column_names = FALSE,
        column_order = cell_patient_order,
        column_split = data.frame(
            cell_annot_sketch[rownames(mat), ]$Subgroup,
            cell_annot_sketch[rownames(mat), ]$Patient_ID
        ),
        row_labels = metacom_human_names,
        row_names_gp = gpar(fontsize = 8),
        show_row_names = TRUE,
        column_names_side = "bottom",
        top_annotation = top_annotation,
        column_title_gp = gpar(fontsize = 8),
        column_title_rot = 45
    )
    
    
}

diax_heat <- draw_and_save_heat_cond2(diax_mat)
relap_heat <- draw_and_save_heat_cond2(relapse_mat)

png(
    filename = "results/figures/aml_diagnosis_metacom_heat.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(diax_heat)

dev.off()

png(
    filename = "results/figures/aml_relapse_metacom_heat.png",
    width = 19,
    height = 4,
    units = "in",
    res = 200
)
draw(relap_heat)

dev.off()


## overlap analysis
metacom_untreated <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
    )

studies <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    select(IDs, studies) %>%
    deframe()

metacom_untreated$study <- studies[metacom_untreated$signature]
write.table(x = metacom_untreated, file = "annot.tsv")

metacom_untreated_drugs <- split(
    metacom_untreated$signature,
    metacom_untreated$meta_community
)

metacom_treated <- read.table(
    file = "results/modules/annotated/metagroup_patients_treated_consensus_drugs.tsv"
    )

metacom_treated_drugs <- split(
    metacom_treated$signature,
    metacom_treated$meta_community
)

## unt6 tt4
common_drugs <- intersect(metacom_untreated_drugs[[5]], metacom_treated_drugs[[5]])


