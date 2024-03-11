library(BPCells)
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


## Somehow there are ~10 NA rows maybe after binding
meta.data <- bc@meta.data %>%
    filter(!is.na(sample)) %>%
    mutate(
        sample_study = paste0(sample, study, sep = "_")
    )

metacom_types <- meta.data %>%
    select(
        cell,
        new_cell_id,
        sample,
        tumor_type,
        sample_type,
        treated,
        study,
        sample_study,
        metacom_untreated_1:metacom_untreated_6
        )

## Export a reduced matrix for potential analyses
write_tsv(
    x = metacom_types,
    file = "results/modules/annotated/pancancer_modules_mt1.tsv"
    )

## Aggregate all enrichment scores by sample
metacom_sample_agg <- meta.data %>%
    group_by(sample_study) %>%
    mutate(
        metacom_1_overall = sum(metacom_untreated_1),
        metacom_2_overall = sum(metacom_untreated_2),
        metacom_3_overall = sum(metacom_untreated_3),
        metacom_4_overall = sum(metacom_untreated_4),
        metacom_5_overall = sum(metacom_untreated_5),
        metacom_6_overall = sum(metacom_untreated_6)
        
    )

## Calculate highest metacom for each cell
metacom_types_highest_cell <- metacom_sample_agg %>%
    select(
        cell,
        new_cell_id,
        sample,
        tumor_type,
        sample_type,
        treated,
        study,
        sample_study,
        metacom_untreated_1:metacom_untreated_6,
    ) %>%
    pivot_longer(
        cols = metacom_untreated_1:metacom_untreated_6,
        names_to = "metacommunity",
        values_to = "cell_enrichment"
    ) %>%
    group_by(new_cell_id) %>%
    slice_max(order_by = cell_enrichment, n = 1)


metacom_proportions <- metacom_types_highest_cell %>%
    group_by(sample_study, metacommunity) %>%
    reframe(
        sample = sample,
        study = study,
        n.cells = n()
    ) %>%
    distinct() %>%
    group_by(sample_study) %>%
    mutate(
        n.total = sum(n.cells),
        n.prop = n.cells / n.total,
        n.prop = round(n.prop, digits = 3)
    ) %>%
    arrange(desc(n.prop)) %>%
    group_by(sample_study) %>%
    mutate(
        best_metacom = head(metacommunity, n = 1)
    )

metacom_sample_agg <- metacom_sample_agg %>%
    select(sample, study, sample_study, metacom_1_overall:metacom_6_overall) %>%
    distinct()

## add cancer types
cancer_info <- meta.data %>%
    select("sample_study", "sample_type", "treated", "tumor_type", "tumor_subtype") %>%
    distinct()

metacom_proportions_annotated <- metacom_proportions %>%
    left_join(
        y = cancer_info,
        by = "sample_study"
    ) %>%
    distinct()

## sample sum information
metacom_sums_by_sample <- metacom_sample_agg %>%
    select(sample_study, metacom_1_overall:metacom_6_overall)

## excelize for A.F. and add sample enrichments
metacom_proportions_annotated_wide <- metacom_proportions_annotated %>%
    pivot_wider(
        id_cols = c(
            "sample",
            "study",
            "sample_study",
            "tumor_type",
            "tumor_subtype",
            "sample_type",
            "treated",
            "best_metacom",
            "n.total"
            ),
        names_from = c("metacommunity"),
        values_from = c("n.cells", "n.prop")
        ) %>%
    left_join(
        y = metacom_sums_by_sample,
        by = "sample_study"
    ) 

## get the name of the metacom with the highest enrichment
metacom_sample_agg_long <- metacom_sample_agg %>%
    pivot_longer(
        cols = metacom_1_overall:metacom_6_overall,
        names_to = "metacommunity",
        values_to = "enrichment"
    ) %>%
    group_by(sample_study) %>%
    slice_max(order_by = enrichment, n = 1)

metacom_proportions_annotated_wide <- metacom_proportions_annotated_wide %>%
    left_join(
        y = metacom_sample_agg_long[, c("sample_study", "metacommunity")],
        by = "sample_study"
    ) %>%
    dplyr::rename(
        "best_metacom_enrichment" = metacommunity
    )

write_tsv(
    x = metacom_proportions_annotated_wide,
    file = "results/modules/annotated/metacom_proportions_primary_wide.tsv"
    )


## metacom best metacom by cancer type and sample
metacom_proportions_annotated_sum <- metacom_proportions_annotated_wide %>%
    select(
        sample,
        tumor_type,
        tumor_subtype,
        sample_type,
        study,
        best_metacom,
        best_metacom_enrichment
        ) %>%
    distinct()

metacom_sample_best_by_freq <- table(
    metacom_proportions_annotated_sum$tumor_type,
    metacom_proportions_annotated_sum$sample_type,
    metacom_proportions_annotated_sum$best_metacom
    )

metacom_sample_best_by_freq <- as.data.frame(metacom_sample_best_by_freq)
colnames(metacom_sample_best_by_freq) <- c(
    "Tumor type",
    "Sample type",
    "Metacommunity",
    "N.samples"
    )

metacom_sample_best_wide <- metacom_sample_best_by_freq %>%
    pivot_wider(
        id_cols = c("Tumor type", "Sample type"),
        names_from = "Metacommunity",
        values_from = "N.samples"
    ) 

write_tsv(
    x = metacom_sample_best_wide,
    file = "results/modules/annotated/metacom_bests_primary_by_freq.tsv"
    )

brca_tumor_best_metacoms <- metacom_proportions_annotated_sum %>%
    filter(tumor_type == "BRCA")

brca_tumor_best_metacoms <- table(
    brca_tumor_best_metacoms$tumor_subtype,
    brca_tumor_best_metacoms$best_metacom
)

brca_tumor_best_metacoms <- as.data.frame(brca_tumor_best_metacoms)
colnames(brca_tumor_best_metacoms) <- c("Subtype", "Best metacommunity", "N.samples")

brca_tumor_best_metacoms_wide <- brca_tumor_best_metacoms %>%
    pivot_wider(id_cols = "Subtype", names_from = "Best metacommunity", values_from = "N.samples")

write_tsv(
    x = brca_tumor_best_metacoms_wide,
    file = "results/modules/annotated/brca_subtypes_metacommunities_untreated.tsv"
    )
