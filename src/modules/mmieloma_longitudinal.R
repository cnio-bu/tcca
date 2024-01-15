library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggpubr)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = 'v5')

# 
# mmieloma_communities <- data.table::fread(
#     "results/modules/annotated/MM_treated_communities.tsv"
#     )
# 
# clinical_data = data.table::fread("results/annotation/clinical_metadata_v3_clean.tsv")
# 
# 
# mielomas <- clinical_data %>%
#     filter(tumor_type == "MM", treated == "t" )
# 
# rm006 <- clinical_data %>%
#     filter(patient == "RRMM06")
# 
# samples_to_use <- rm006 %>%
#     pull(sample)
# 
# drug_data <- data.table::fread("reference/final_moas - Collapsed.tsv")
# 
# ## WHICH COMMS ARE RRMM06
# rm006_communities <- mmieloma_communities %>%
#     filter(sample %in% samples_to_use) %>%
#     group_by(sample) %>%
#     reframe(communities = community) %>%
#     distinct()
# 
# 
# comm_summary <- mmieloma_communities %>%
#     group_by(community, signature) %>%
#     mutate(
#         recurrency = n()
#     ) %>%
#     filter(community == 5) %>%
#     left_join(y = drug_data[,c("IDs", "studies", "preferred.drug.names")],
#               by = c("signature" = "IDs"), multiple = "first"
#               )
# 

## Community 5 top Hit is a ROSN inhibitor, known to happen in MM (ROS)
## It is also prominent in RM006 T3. Habria que ver la historia BRD-K71935468


#### NECESITO
###
### Tratamientos anotados del paciente
### Lineas de respuesta
### beyondcell completo para rrmm006 aunque NO esten todas las celulas en bc-meta

#Preguntas: 
#    1. Que modulos estan representados en cada sample de RRMM066
#    2. Miramos las drogas. Casan ?
#    3. Superimponer los modulos.
#    4. Superimponer metacomunidades

drugs_metacommunities_treated <- read.table(
    "results/modules/annotated/metagroup_patients_treated_consensus_drugs.tsv"
    )

meta_coms_set <- split(
    drugs_metacommunities_treated$signature,
    drugs_metacommunities_treated$meta_community
    )


bc <- readRDS("results/mmieloma/bc_seu.Rds")

## Manually manipulate pointer
bc@assays$RNA$counts@matrix@matrix@dir <- c("./results/mmieloma/mmieloma_bc")
bc@assays$RNA$data@matrix@matrix@dir <- c("./results/mmieloma/mmieloma_bc/")

DefaultAssay(bc) <- "RNA"

bc <- AddModuleScore(
    object = bc,
    features = meta_coms_set,
    name = "metacom_",
    seed = 120394,
    slot = "data",
    ctrl = 10
    )

bcrm13 <- subset(bc, subset = PID_new == "RRMM13")

bc13mat <- as.matrix(bcrm13@assays$RNA$data)
bc13mat <- scale(bc13mat, center = TRUE, scale = TRUE)

top_rv <- matrixStats::rowVars(bc13mat)
top_rv <- top_rv[top_rv >= 2]


## load drug data
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    distinct() %>%
    as.data.frame()

rownames(drugs) <- drugs$IDs

cell_annot_df <- bcrm13@meta.data[, c(
    "timepoint",
    "sample_id",
    "sc_gain_1q",
    "metacom_1",
    "metacom_2",
    "metacom_3",
    "metacom_4",
    "metacom_5",
    "metacom_6"
    )]

cell_annot_df$bortezo <- bc13mat["sig-21377", rownames(cell_annot_df)]
cell_annot_df$timepoint <- as.factor(cell_annot_df$timepoint)
cell_annot_df$new_time <- fct_relevel(cell_annot_df$timepoint, "pre", "post")

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "timepoint" = cell_annot_df[, c("new_time")],
    "bortezo" = anno_barplot(cell_annot_df[, "bortezo"]),
    "1q amplification" = cell_annot_df[, c("sc_gain_1q")],
    which = "row",
   # col = pals,
 #   annotation_name_side = "top",
    annotation_name_rot = 45
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "meta_community" = as.factor(drugs_metacommunities_treated[, c("meta_community")]),
    which = "column",
    annotation_name_side = "left"
)

a <- ComplexHeatmap::Heatmap(
    mat = t(bc13mat[drugs_metacommunities_treated$signature,]),
    cluster_rows = FALSE,
 #   row_order = rownames(cell_annot_df[order(cell_annot_df$new_time), ]),
    row_split = cell_annot_df$new_time,
    cluster_row_slices = TRUE,
    clustering_distance_rows = "pearson",
    #cluster_columns = TRUE,
    cluster_columns =  FALSE,
    cluster_column_slices = TRUE,
    column_split = drugs_metacommunities_treated$meta_community,
    clustering_distance_columns = "pearson",
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 4),
    column_names_side = "top",
    column_title = NULL,
    heatmap_width = unit(8, "in"),
    heatmap_height = unit(14, "in"),
    column_labels = drugs_metacommunities_treated$preferred.drug.names,
    right_annotation = right_annotation,
    top_annotation = top_annotation
)

png(filename = "results/rrmm_13_top40.png",
    res = 300,
    width = 14, 
    height = 18,
    units = "in",
    )
a
dev.off()

bortezo_change <- ggplot(
    data = cell_annot_df,
    aes(y = bortezo, fill = as.factor(sc_gain_1q))
    ) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_discrete(name = "1q amplified", labels = c("No", "Yes")) +
    scale_y_continuous(name = "Bortezomib bcScore", n.breaks = 10) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

ggsave(
    plot = bortezo_change,
    filename = "results/figures/bortezomib_rm13_1q_amplification.png",
    dpi = 100,
    height = 7,
    width = 7
    )


metacom_enrichment = cell_annot_df %>%
    select(new_time, sample_id, metacom_1:metacom_6) %>%
    pivot_longer(
        cols = metacom_1:metacom_6,
        names_to = "metacommunity",
        values_to = "enrichment"
        )

metacom_change <- ggplot(
    data = metacom_enrichment,
    aes(fill = new_time, y = enrichment, x = metacommunity)
) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test") +
    scale_y_continuous(limits = c(-6,6)) +
    scale_x_discrete(
        name = "",
        labels = paste0("Meta community", " ", rep(1:6))
        ) +
    scale_fill_discrete(
        name = "Timepoint",
        labels = c("pre-treatment", "post-treatment")
        ) +
    ylab("Module enrichment score") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )

ggsave(
    plot = metacom_change,
    filename = "results/figures/treated_metacommunities_rm013_module.png",
    dpi = 100,
    height = 14,
    width = 14
    )

## HEATMAP DE LOS MODULE SCORES!
module_mat <- bcrm13@meta.data %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, metacom_1:metacom_6) %>%
    as.data.frame()

rownames(module_mat) <- module_mat$cell_barcode
module_mat$cell_barcode <- NULL
module_mat <- as.matrix(module_mat)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "timepoint" = cell_annot_df[, c("new_time")],
    "1q amplification" = cell_annot_df[, c("sc_gain_1q")],
    which = "column"
    # col = pals,
    #   annotation_name_side = "top",
)


b <- ComplexHeatmap::Heatmap(
    mat = t(module_mat),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_order = rownames(cell_annot_df[order(cell_annot_df$new_time), ]),
    cluster_column_slices = TRUE,
    clustering_distance_columns = "pearson",
    column_split = cell_annot_df$new_time,
    show_column_names = FALSE,
    top_annotation = top_annotation
)

## All patients sketch by metacom variability
bc <- SketchData(
    object = bc,
    assay = "RNA",
    ncells = 10000,
    sketched.assay = "sketch_10k_new",
    method = "LeverageScore",
    var.name = "leveragev2",
    seed = 120394,
    verbose = TRUE,
    over.write = TRUE
)

DefaultAssay(bc) <- "sketch_10k_new"


## Add treatment groups to bc
additional_meta <- data.table::fread(
    "reference/additional_metadata_mmieloma.tsv"
    ) %>%
    mutate(
        treatment_group = replace_na(treatment_group, "None") 
    ) %>%
    select(-"1q_state") %>%
    as.data.frame()

bc_annotated_meta <- bc@meta.data %>%
    rownames_to_column("rownames") %>%
    left_join(y = additional_meta, by = c("PID_new" = "sample_id")) %>%
    as.data.frame()

rownames(bc_annotated_meta) <- bc_annotated_meta$rownames
bc_annotated_meta$rownames <- NULL
bc@meta.data <- bc_annotated_meta

table(bc@meta.data[colnames(bc@assays$sketch_10k_new), "treatment_group"])

## meki patients
## HEATMAP DE LOS MODULE SCORES!
module_mat <- bc@meta.data[colnames(bc@assays$sketch_10k_new$counts), ] %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, metacom_1:metacom_6, treatment_group) %>%
    as.data.frame()

meki_mat <- module_mat %>%
    filter(treatment_group == "PI") %>%
    select(-treatment_group) %>%
    as.data.frame()

rownames(meki_mat) <- meki_mat$cell_barcode
meki_mat$cell_barcode <- NULL
meki_mat <- as.matrix(meki_mat)

## MEKI multi-patient HEAT
cell_annot_df <- bc@meta.data[colnames(bc@assays$sketch_10k_new$counts), c(
    "PID_new",
    "timepoint",
    "sc_gain_1q",
    "drug_t1_response"
)]

cell_annot_df$timepoint <- as.factor(cell_annot_df$timepoint)
cell_annot_df$new_time <- fct_relevel(cell_annot_df$timepoint, "pre", "post", "post_2")

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Timepoint" = cell_annot_df[rownames(meki_mat), c("new_time")],
    "1q amplification" = cell_annot_df[rownames(meki_mat), c("sc_gain_1q")],
    "Treatment response" = cell_annot_df[rownames(meki_mat), c("drug_t1_response")],
    which = "column"
    # col = pals,
    #   annotation_name_side = "top",
)

cell_patient_order <- cell_annot_df[rownames(meki_mat), ]
cell_patient_order <- rownames(cell_patient_order[order(cell_patient_order$PID_new), ])

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

b <- ComplexHeatmap::Heatmap(
    mat = t(meki_mat),
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_order = cell_patient_order,
    cluster_column_slices = FALSE,
    cluster_row_slices = TRUE,
    clustering_distance_columns = "pearson",
    column_split = cell_annot_df[rownames(meki_mat), ]$PID_new,
    column_gap = unit(2, "mm"),
    show_column_names = FALSE,
  #  row_labels = paste0("Metacommunity ", c(1:6)),
    row_names_side = "left",
    top_annotation = top_annotation
)

