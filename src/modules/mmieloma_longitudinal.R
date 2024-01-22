library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggpubr)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")

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


## Test, for each patient, the overall number of features and cells
sample_mat_averages <- bc@assays$RNA$data
sample_mat_averages <- colMeans(sample_mat_averages)

sample_mat_averages_annotated <- sample_mat_averages %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    rename("average_bcscore" = ".") %>%
    left_join(y = bc@meta.data[, c(
        "Cell_barcode",
        "PID_new",
        "PID_sample_new",
        "timepoint")],
        by = c("cell_id" = "Cell_barcode")
        )

## Add treatment groups to bc
additional_meta <- data.table::fread(
    "reference/additional_metadata_mmieloma.tsv"
) %>%
    mutate(
        treatment_group = replace_na(treatment_group, "None") 
    ) %>%
    select(-"1q_state") %>%
    as.data.frame()

sample_mat_averages_annotated <- sample_mat_averages_annotated %>%
    left_join(y = additional_meta, by = c("PID_new" = "sample_id"))


averages_bscore_ggplot <- ggplot(
    data = sample_mat_averages_annotated,
    aes(
        x = PID_new,
        y = average_bcscore,
        fill = timepoint
        )
    ) +
    geom_boxplot() +
    facet_grid(~treatment_group) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        legend.title = element_text(face = "bold")
    )

ggsave(
    filename = "results/mmieloma/raw_bscores_averages.png",
    averages_bscore_ggplot,
    dpi = 100, height = 7,
    width = 14
    )

## Regress out patient specific effects
bc_mat <- bc@assays$RNA$data

bc <- Seurat::ScaleData(object = bc, vars.to.regress = "PID_new")

bc <- AddModuleScore(
    object = bc,
    features = meta_coms_set,
    name = "metacom_",
    seed = 120394,
    slot = "scale.data",
    ctrl = 10
)


## load drug data
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    distinct() %>%
    as.data.frame()

rownames(drugs) <- drugs$IDs

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

## Scale and regress the sketch too
bc <- Seurat::ScaleData(object = bc, vars.to.regress = "PID_new")

## Add annotation of treatment groups to bc object
bc_annotated_meta <- bc@meta.data %>%
    rownames_to_column("rownames") %>%
    left_join(y = additional_meta, by = c("PID_new" = "sample_id")) %>%
    as.data.frame()

rownames(bc_annotated_meta) <- bc_annotated_meta$rownames
bc_annotated_meta$rownames <- NULL
bc@meta.data <- bc_annotated_meta

## Export sketched mats for each treatment group to draw a heatmap
module_mat <- bc@meta.data[colnames(bc@assays$sketch_10k_new$scale.data), ] %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, metacom_1:metacom_6, treatment_group) %>%
    as.data.frame()

meki_mat <- module_mat %>%
    filter(treatment_group == "MEKi") %>%
    select(-treatment_group) %>%
    as.data.frame()

rownames(meki_mat) <- meki_mat$cell_barcode
meki_mat$cell_barcode <- NULL
meki_mat <- as.matrix(meki_mat)

imid_mat <- module_mat %>%
    filter(treatment_group == "IMiD") %>%
    select(-treatment_group) %>%
    as.data.frame()

rownames(imid_mat) <- imid_mat$cell_barcode
imid_mat$cell_barcode <- NULL
imid_mat <- as.matrix(imid_mat)

pi_mat <- module_mat %>%
    filter(treatment_group == "PI") %>%
    select(-treatment_group) %>%
    as.data.frame()

rownames(pi_mat) <- pi_mat$cell_barcode
pi_mat$cell_barcode <- NULL
pi_mat <- as.matrix(pi_mat)


## MEKI multi-patient HEAT
cell_annot_df <- bc@meta.data[colnames(bc@assays$RNA$counts), c(
    "PID_new",
    "timepoint",
    "sc_gain_1q",
    "drug_t1_response"
)]

cell_annot_df$timepoint <- as.factor(cell_annot_df$timepoint)
cell_annot_df$new_time <- fct_relevel(
    cell_annot_df$timepoint,
    "pre", "post", "post_2"
    )

# Write annotation to disk
write.table(x = cell_annot_df, file = "results/mmieloma/cell_annotation.tsv")

# Write matrices to disk
write.table(x = meki_mat, file = "results/mmieloma/meki_mat_enrichment.tsv")
write.table(x = pi_mat, file = "results/mmieloma/pi_mat_enrichment.tsv")
write.table(x = imid_mat, file = "results/mmieloma/imid_mat_enrichment.tsv")

# Export full metacom mat with 90k cells too
full_mat_90k <-  bc@meta.data[colnames(bc@assays$RNA$scale.data), ] %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, metacom_1:metacom_6, treatment_group) %>%
    as.data.frame()

write.table(x = full_mat_90k, file = "results/mmieloma/full_mat_metacom.tsv")

## Generate one last plot after regressing out the patient effect
sample_mat_averages_rgrss <- bc@assays$RNA$scale.data
sample_mat_averages_rgrss <- colMeans(sample_mat_averages_rgrss)

sample_mat_averages_annotated_rgrss <- sample_mat_averages_rgrss %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    rename("average_bcscore" = ".") %>%
    left_join(y = bc@meta.data[, c(
        "Cell_barcode",
        "PID_new",
        "PID_sample_new",
        "timepoint",
        "treatment_group")],
        by = c("cell_id" = "Cell_barcode")
    )


averages_bscore_ggplot <- ggplot(
    data = sample_mat_averages_annotated_rgrss,
    aes(
        x = PID_new,
        y = average_bcscore,
        fill = timepoint
    )
) +
    geom_boxplot() +
    facet_grid(~treatment_group) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        legend.title = element_text(face = "bold")
    )

ggsave(
    filename = "results/mmieloma/raw_bscores_averages_rgrss.png",
    averages_bscore_ggplot,
    dpi = 100, height = 7,
    width = 14
)
