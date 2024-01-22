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
    file = "results/mmieloma/cell_annotation.tsv"
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
meki_heat <- draw_and_save_heat(meki_mat)

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

## Calculate the shift in metacom. enrichment between responses
full_mat <- read.table("results/mmieloma/full_mat_metacom.tsv") %>%
    rename("cell_id" = cell_barcode) %>%
    as.data.frame()

cell_annot_df <- cell_annot_df %>%
    rownames_to_column("cell_id")

full_mat_annotated <- full_mat %>%
    left_join(y = cell_annot_df, by = "cell_id") %>%
    pivot_longer(
        cols = metacom_1:metacom_6,
        names_to = "community",
        values_to = "enrichment"
        ) %>%
    filter(new_time %in% c("Pretreated", "First line")) %>%
    mutate(
        drug_t1_response = as_factor(drug_t1_response)
    )


full_mat_annotated$drug_t1_response <- fct_relevel(
    full_mat_annotated$drug_t1_response,
    "CR", "VGPR", "PR", "MR", "SD"
    )

## Add previous therapies n.
additional_meta <- read_tsv("reference/additional_metadata_mmieloma.tsv") %>%
    select(sample_id, prev_therapies)

full_mat_annotated <- full_mat_annotated %>%
    left_join(additional_meta, by = c("PID_new" = "sample_id"))

full_mat_annotated$treatment_group <- as.factor(full_mat_annotated$treatment_group)
full_mat_annotated$prev_therapies <- as.factor(full_mat_annotated$prev_therapies)

disp_plot <- ggplot(
    data = full_mat_annotated[!(is.na(full_mat_annotated$drug_t1_response)), ],
    aes(x = community, y = enrichment, fill = new_time)) +
    geom_boxplot() +
    scale_fill_discrete(name = "Timepoint") +
    scale_x_discrete(labels = paste0("Metacomunity ", rep(1:6))) +
    ylab("") +
    xlab("") +
    facet_wrap(treatment_group ~ drug_t1_response) + 
    stat_compare_means(method = "wilcox.test", na.rm = TRUE, label = "p.signif") +
    theme_bw() + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        legend.title = element_text(face = "bold")
    )
        

ggsave(
    filename = "results/figures/mmieloma_metacom_timepoint_by_response_timegroup.png",
    plot = disp_plot,
    dpi = 100,
    width = 14,
    height = 14
    )

disp_plot_groups <- ggplot(
    data = full_mat_annotated[!(is.na(full_mat_annotated$drug_t1_response)), ],
    aes(x = community, y = enrichment, fill = new_time)) +
    geom_boxplot() +
    scale_x_discrete(name = "", labels = paste0("Meta community ", rep(1:6))) +
    scale_fill_discrete(name = "Timepoint") +
    scale_y_continuous(name = "", limits = c(-10, 10), n.breaks = 10) +
    facet_wrap(~treatment_group) + 
    stat_compare_means(
        method = "wilcox.test",
        na.rm = TRUE,
        label = "p.signif",
        ) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        legend.title = element_text(face = "bold")
    )


ggsave(
    filename = "results/figures/mmieloma_metacom_treatments_by_group.png",
    plot = disp_plot_groups,
    dpi = 100,
    width = 21,
    height = 14
)


## study metacom change by 1a amp.
disp_plot_amp <- ggplot(
    data = full_mat_annotated[!(is.na(full_mat_annotated$drug_t1_response)), ],
    aes(x = community, y = enrichment, fill = sc_gain_1q)) +
    geom_boxplot() +
    scale_fill_discrete(name = "1q. status") +
    scale_x_discrete(labels = paste0("Metacomunity ", rep(1:6))) +
    ylab("") +
    xlab("") +
    facet_wrap(treatment_group ~ new_time, ncol = 2) + 
    stat_compare_means(method = "wilcox.test", na.rm = TRUE, label = "p.signif") +
    theme_bw() + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        #panel.background = element_blank(),
        #axis.line = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        legend.title = element_text(face = "bold")
    )

ggsave(
    plot = disp_plot_amp,
    filename = "results/figures/mmieloma_metacom_timepoint_amp.png",
    dpi = 100,
    width = 7,
    height = 14
    )
