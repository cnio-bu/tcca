library(vegan)
library(tidyverse)

ith_by_sample <- data.table::fread(
    input  = "results/modules/annotated/metacom_proportions_primary_wide.tsv"
    )


## start with naive primaries
ith_primaries_naive <- ith_by_sample %>%
    filter(
        study != "cell_lines_gabriella_kinker",
   #     sample_type == "p",
    #    treated == FALSE
        ) %>%
    filter(
        !(study == "brca_bhupinder_pal" & tumor_subtype == "predicted_tumour")
    ) %>%
    select(
        sample_study,
        tumor_type,
        sample_type,
        treated,
        n.prop_metacom_untreated_1,
        n.prop_metacom_untreated_2,
        n.prop_metacom_untreated_3,
        n.prop_metacom_untreated_4,
        n.prop_metacom_untreated_5,
        n.prop_metacom_untreated_6
    ) 

## comm data
ith_primaries_naive_mat <- ith_primaries_naive %>%
    select(sample_study, n.prop_metacom_untreated_1:n.prop_metacom_untreated_6) %>%
    as.data.frame()

rownames(ith_primaries_naive_mat) <- ith_primaries_naive_mat$sample_study
ith_primaries_naive_mat$sample_study <- NULL

ith_primaries_naive_mat <- as.matrix(ith_primaries_naive_mat)
ith_primaries_naive_mat[is.na(ith_primaries_naive_mat)] <- 0

shan <- vegan::diversity(
    x = ith_primaries_naive_mat,
    index = "shannon",
    MARGIN = 1,
 #   groups = ith_primaries_naive$tumor_type,
  #  equalize.groups = TRUE
)
ith_primaries_naive$shan <- shan
ith_primaries_naive$tumor_type <- fct_reorder(ith_primaries_naive$tumor_type, ith_primaries_naive$shan)

## genomic ith
genomic_ith <- read_tsv("results/cna/genomic_ith.tsv")

## add n.cells by sample
tcs <- read_tsv("results/annotation/beyondcell_with_therapeutic_clusters.tsv") %>%
    mutate(
        sample_study = paste0(study, "__", sample)
    ) %>%
    group_by(sample_study) %>%
    mutate(
        n.cells = n()
    ) %>%
    filter(
        !is.na(cell)
    )

genomic_ith <- genomic_ith %>%
    left_join(
        y = tcs[, c("sample_study", "n.cells")],
        by = c("study__sample" = "sample_study")
    ) %>%
    distinct() 

genomic_ith_rates <- genomic_ith %>%
    mutate(
        n.prop = nclones * 1000 / n.cells,
        sample_study = paste0(sample, study, sep = "_")
    )

ith_primaries_naive$shan <- shan

genomic_therapeutic_ith <- ith_primaries_naive %>%
    left_join(
        y = genomic_ith_rates,
        by = "sample_study"
    ) %>%
    filter(
        !is.na(n.prop),
        !is.na(shan)
    )

## test
genomic_therapeutic_ith <- genomic_therapeutic_ith %>%
    filter(
    #    sample_type == "p",
    #       treated == FALSE
          )

test <- ggplot(data = genomic_therapeutic_ith, aes(x = shan, y = n.prop)) +
    geom_point() +
    geom_smooth(method = "lm") +
    ggpubr::stat_cor() +
    facet_wrap(~tumor_type, nrow=8, ncol=8)
