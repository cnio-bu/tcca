library(ggalluvial)
library(tidyverse)


# ggplot(data = gg_data,
#        aes(axis1 = models, axis2 = clusters)) +
#     geom_alluvium(aes(fill = models), width = 1/12) +
#     geom_stratum(width = 1/12, fill = "black", color = "grey") +
#     geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#     scale_x_discrete(limits = c("models", "clusters"), expand = c(.05, .05)) +
#     scale_fill_brewer(type = "qual", palette = "Set1") +
#     ggtitle("Distribución de models en los clusters")


metagroup_drugset_1 <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
    ) %>%
    group_by(signature) %>%
    slice_max(order_by = n.appearances, n = 1, with_ties = FALSE) %>%
    mutate(
        metagroup = "MT1"
    ) %>%
    select(
        meta_community, signature, metagroup
    )

metagroup_drugset_2 <- read.table(
    "results/modules/annotated/metagroup_patients_treated_consensus_drugs.tsv"
    ) %>%
    group_by(signature) %>%
    slice_max(order_by = n.appearances, n = 1, with_ties = FALSE) %>%
    mutate(
        metagroup = "MT2"
    ) %>%
    select(
        meta_community, signature, metagroup
    )



metagroups <- bind_rows(metagroup_drugset_1, metagroup_drugset_2)  %>%
    mutate(
        meta_community = as.character(meta_community)
    ) %>%
    pivot_wider(
        id_cols = "signature",
        names_from = "metagroup",
        values_from = "meta_community"
    ) %>%
    replace_na(list(MT1 = "mt2_unique", MT2 = "mt1_unique")) 
    mutate(
        MT1 = paste0("1_", MT1),
        MT2 = paste0("2_", MT2)
    )

sankey <- ggplot(data = metagroups,
                 aes(axis1 = MT1, axis2 = MT2)) +
    geom_alluvium(aes(fill = MT1), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("MT1", "MT2"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme_minimal()
    ggtitle("Distribución de models en los clusters")


