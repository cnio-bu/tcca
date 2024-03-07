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
    filter(
        sample_type == "p" &
            treated == FALSE &
            study != "cell_lines_gabriella_kinker"
        ) %>%
    pivot_longer(
        cols = metacom_untreated_1:metacom_untreated_6,
        names_to = "metacommunity",
        values_to = "cell_enrichment"
        ) %>%
    group_by(sample, study, metacommunity) %>%
    reframe(
        sample_enrichment = median(cell_enrichment),
        tumor_type = tumor_type
    ) %>%
    distinct()
    

metacom_types <- metacom_types %>%
    group_by(tumor_type, metacommunity) %>%
    mutate(n_samples = n()) %>%
    filter(
        n_samples >= 3
    ) %>%
    mutate(tumor_type = as_factor(tumor_type)) 

metacom_summary <- metacom_types %>%
    group_by(metacommunity) %>%
    mutate(
        pancancer_median = round(median(sample_enrichment), digits = 3)
    ) %>%
    group_by(tumor_type, metacommunity) %>%
    reframe(
        average_enrichment = round(median(sample_enrichment), digits = 3),
        enrichment_disp = round(sd(sample_enrichment), digits = 3),
        pancancer_median = pancancer_median
    ) %>%
    distinct() %>%
    arrange(
        desc(average_enrichment)
    )


write_tsv(
    x = metacom_summary,
    file = "results/modules/annotated/metacommunity_untreated_enrichment_by_cancer.tsv"
    )


## metacom types 2
metacom_types_highest_cell <- meta.data %>%
    select(
        cell,
        sample,
        tumor_type,
        sample_type,
        treated,
        study,
        metacom_untreated_1:metacom_untreated_6
    ) %>%
    filter(
        sample_type == "p" &
            treated == FALSE &
            study != "cell_lines_gabriella_kinker"
    ) %>%
    pivot_longer(
        cols = metacom_untreated_1:metacom_untreated_6,
        names_to = "metacommunity",
        values_to = "cell_enrichment"
    ) %>%
    group_by(cell) %>%
    slice_max(order_by = cell_enrichment, n = 1)


metacom_proportions <- metacom_types_highest_cell %>%
    group_by(sample, study, metacommunity) %>%
    summarise(
        n.cells = n()
    ) %>%
    group_by(sample) %>%
    mutate(
        n.total = sum(n.cells),
        n.prop = n.cells / n.total * 100,
        n.prop = round(n.prop, digits = 2)
    ) %>%
    arrange(desc(n.prop)) %>%
    group_by(sample) %>%
    mutate(
        best_metacom = head(metacommunity, n = 1)
    )


## add cancer types
metacom_proportions_annotated <- metacom_proportions %>%
    left_join(
        y = meta.data[, c("sample", "study", "sample_type", "tumor_type")],
        by = c("sample", "study")
    ) %>%
    distinct()
    
