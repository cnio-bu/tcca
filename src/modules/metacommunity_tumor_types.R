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
        y = meta.data[, c("sample", "study", "sample_type", "tumor_type", "tumor_subtype")],
        by = c("sample", "study")
    ) %>%
    distinct() %>%
    group_by(study, sample) %>%
    mutate(
        metacom_sd = sd(n.prop)
    )

metacom_dispersion_ctype <- metacom_proportions_annotated %>%
    group_by(tumor_type, metacommunity) %>%
    summarise(
        avg_prop = median(n.prop),
        prop_sd = sd(n.prop),
        tumor_subtype = tumor_subtype
    )

metacom_densities <- ggplot(
    data = metacom_dispersion_ctype,
    aes(x = prop_sd, color = metacommunity)
    ) +
    geom_density(alpha = 0.3) +
    scale_x_continuous(name = "Standard deviation of proportions") +
    ylab("") +
    scale_color_discrete(labels = paste("Metacommunity", c(1:6), sep = " ")) +
    theme_minimal() +
    ggtitle("Pancancer dispersion of cell proportions by metacommunity")


metacom_averages <- ggplot(
    data = metacom_dispersion_ctype,
    aes(x = avg_prop, color = metacommunity)
) +
    geom_density(alpha = 0.3) +
    scale_x_continuous(
        name = "Average proportions",
       # labels = scales::percent(100, scale = 1)
        ) +
    ylab("") +
    scale_color_discrete(labels = paste("Metacommunity", c(1:6), sep = " ")) +
    theme_minimal() +
    ggtitle("Pancancer averages of cell proportions by metacommunity")


## excelize for A.F.
metacom_proportions_annotated_wide <- metacom_proportions_annotated %>%
    pivot_wider(
        id_cols = c("sample", "study", "tumor_type", "tumor_subtype"),
        names_from = c("metacommunity"),
        values_from = c("n.cells", "n.total", "n.prop", "metacom_sd")
        )


write_tsv(
    x = metacom_proportions_annotated_wide,
    file = "results/modules/annotated/metacom_proportions_primary_wide.tsv"
    )


## metacom best metacom by cancer type and sample
metacom_proportions_annotated_sum <- metacom_proportions_annotated %>%
    select(sample, tumor_type, tumor_subtype, study, best_metacom) %>%
    distinct()

metacom_sample_best <- table(
    metacom_proportions_annotated_sum$tumor_type,
    metacom_proportions_annotated_sum$best_metacom
    )

metacom_sample_best <- as.data.frame(metacom_sample_best)
colnames(metacom_sample_best) <- c("Tumor type", "Metacommunity", "N.samples")

metacom_sample_best_wide <- metacom_sample_best %>%
    pivot_wider(
        id_cols = c("Tumor type"),
        names_from = "Metacommunity",
        values_from = "N.samples"
    )

write_tsv(
    x = metacom_sample_best_wide,
    file = "results/modules/annotated/metacom_bests_primary.tsv"
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
