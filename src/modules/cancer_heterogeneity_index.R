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
    select(
        sample_study,
        tumor_type,
        n.prop_metacom_untreated_1,
        n.prop_metacom_untreated_2,
        n.prop_metacom_untreated_3,
        n.prop_metacom_untreated_4,
        n.prop_metacom_untreated_5,
        n.prop_metacom_untreated_6
    ) %>%
    pivot_longer(
        cols = n.prop_metacom_untreated_1:n.prop_metacom_untreated_6,
        names_to = "metacommunity",
        values_to = "n.prop"
    ) %>%
    group_by(
        
    ) %>%
    mutate(
        
    )


dispersion_plot <- ggplot(data = ith_primaries_naive, aes(x = tumor_type, y = n.prop_metacom_untreated_6)) +
    geom_boxplot() +
    theme_minimal()
