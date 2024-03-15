library(vegan)
library(tidyverse)

ith_by_sample <- data.table::fread(
    input  = "results/modules/annotated/metacom_proportions_primary_wide.tsv"
    )


## start with naive primaries
ith_primaries_naive <- ith_by_sample %>%
    filter(
        study != "cell_lines_gabriella_kinker",
        sample_type == "p",
        treated == FALSE
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
    groups = ith_primaries_naive$tumor_type,
    equalize.groups = TRUE
)
ith_primaries_naive$shan <- shan
ith_primaries_naive$tumor_type <- fct_reorder(ith_primaries_naive$tumor_type, ith_primaries_naive$shan)

## test
ith_shan <- ggplot(data = ith_primaries_naive, aes(y = shan, x = tumor_type)) +
    geom_boxplot() +
    geom_point() +
    scale_y_continuous(limits = c(0,2))
