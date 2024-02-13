library(BPCells)
library(ggpubr)
library(Seurat)
library(tidyverse)

bc <- readRDS("results/beyondcell_bp/beyondcell_pancancer.Rds")

metacom_untreated <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
    )

meta_coms_set <- split(
    metacom_untreated$signature,
    metacom_untreated$meta_community
)


bc <- AddModuleScore(
    object = bc,
    features = meta_coms_set,
    seed = 120394,
    slot = "data",
    name = "metacom_untreated_",
    ctrl = 20
)


meta.data <- bc@meta.data

metacom_types <- meta.data %>%
    select(
        cell,
        sample,
        tumor_type,
        sample_type,
        treated,
        study,
        metacom_untreated_1:metacom_untreated_6
        ) %>%
    filter(sample_type == "p" & treated == FALSE) %>%
    group_by(sample, study) %>%
    summarise(
        avg_metacom1 = median(metacom_untreated_1),
        tumor_type = tumor_type
    ) %>%
    ungroup() %>%
    distinct()

## tumors to keep
tumors_to_keep <- table(metacom_types$tumor_type)
tumors_to_keep <- tumors_to_keep[tumors_to_keep > 3]

metacom_types <- metacom_types %>%
    filter(tumor_type %in% names(tumors_to_keep)) %>%
    mutate(tumor_type = as_factor(tumor_type)) 

metacom_types$tumor_type <- fct_reorder(metacom_types$tumor_type, metacom_types$avg_metacom1, .desc = TRUE)

metacom_box <- ggplot(
    data = metacom_types,
    aes(y = avg_metacom1, x = tumor_type)) +
    geom_point() +
    geom_boxplot() +
    scale_y_continuous(n.breaks = 10, limits = c(-1, 1)) +
    stat_compare_means() +
    theme_bw()
