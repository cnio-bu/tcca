library(tidyverse)

metacommunities_mt1 <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
    )

total_moas <- metacommunities_mt1 %>%
    group_by(collapsed.MoAs) %>%
    summarise(
        n.drugs = n()
    )

colnames(total_moas) <- c("collapsed.MoAs", "n.background")

moas_by_metacom <- metacommunities_mt1 %>%
    group_by(collapsed.MoAs, meta_community) %>%
    summarise(
        n.appearances = n()
    ) %>%
    group_by(collapsed.MoAs) %>%
    mutate(
        n.total = sum(n.appearances),
        n.prop = n.appearances / n.total
    ) %>%
    left_join(
        y = total_moas,
        by = c("collapsed.MoAs")
    ) %>%
    pivot_wider(
        id_cols = c("collapsed.MoAs", "n.background"),
        names_from = "meta_community",
        values_from = "n.prop"
    )

write_tsv(
    x = moas_by_metacom, 
    file = "results/modules/annotated/metacommunity_moa_summary.tsv"
    )

