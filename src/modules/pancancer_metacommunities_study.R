library(tidyverse)

pancancer_therapeutic_modules <- read.table(
    file = "results/modules/annotated_hq/patient_primary_untreated_pancancer_communities.tsv",
    sep = "\t",
    header = TRUE
    )

moa_reference <- read.table(
    "reference/final_moas - Collapsed.tsv",
    sep = "\t",
    header = TRUE
    ) %>%
    distinct(IDs, collapsed.MoAs, .keep_all = TRUE)


tm_drugs <- pancancer_therapeutic_modules %>%
    select(signature) %>%
    distinct() %>%
    mutate(
        has_three_underscores = case_when(
            stringr::str_count(string = tm_drugs$signature, pattern = "_") >= 3 ~ TRUE,
            TRUE ~ FALSE
        ),
        signature_fixed = case_when(
            has_three_underscores == TRUE ~ stringr::str_replace(string = tm_drugs$signature, pattern = "_", replacement = "-"),
            TRUE ~ signature
        )
    ) %>%
    separate(
        col = signature_fixed,
        into = c("drug_name", "study", "id"),
        sep = "_"
    )

tm_drugs_moa <- tm_drugs %>%
    left_join(
        y = moa_reference[, c("original.IDs", "studies", "collapsed.MoAs")],
        by = c("id" = "original.IDs", "study" = "studies")
        ) %>%
    select(signature, collapsed.MoAs) %>%
    deframe()

pancancer_therapeutic_modules$moa <- tm_drugs_moa[pancancer_therapeutic_modules$signature]

write.table(
    x = pancancer_therapeutic_modules,
    file = "results/modules/annotated_hq/pancancer_communities_annotated_moa.tsv"
    )

pancancer_drug_composition <- pancancer_therapeutic_modules %>%
    group_by(community, signature) %>%
    mutate(
        n.appearances = n()
    ) %>%
    filter(n.appearances >= 5) %>%
    select(community, signature, moa, n.appearances)


## metacom specificity
pancancer_drug_comp_by_community <- as.data.frame(table(
    pancancer_drug_composition$signature,
    pancancer_drug_composition$community
    )
)
colnames(pancancer_drug_comp_by_community) <- c("drug", "community", "count")

pancancer_drug_comp_by_community <- pancancer_drug_comp_by_community %>%
    pivot_wider(
    id_cols = drug,
    names_from = community,
    values_from = count
)

colnames(pancancer_drug_comp_by_community) <- c("drug", "TM1", "TM2", "TM3")

pancancer_drug_comp_by_community <- pancancer_drug_comp_by_community %>%
    rowwise() %>%
    mutate(
        total_appearances = TM1 + TM2 + TM3
    )

pancancer_drug_comp_by_community$moa <- tm_drugs_moa[as.character(pancancer_drug_comp_by_community$drug)]

write.table(
    x = pancancer_drug_comp_by_community,
    file = "results/modules/annotated_hq/pancancer_drug_composition_by_module.tsv",
    sep = "\t"
    )
