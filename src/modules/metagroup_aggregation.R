library(igraph)
library(tidyverse)

all_modules_mt1 <- list.files(
    path = "results/modules/annotated",
    pattern = "communities.tsv",
    full.names = TRUE
    )


treated_modules <- grep(
    pattern = "_treated_",
    x = all_modules_mt1,
    value = FALSE
    )

all_modules_treated <- all_modules_mt1[treated_modules]
all_modules_non_treated <- all_modules_mt1[-treated_modules]

all_modules_non_treated <- all_modules_non_treated %>% 
    map(read_tsv, id = "module_source")

all_modules_non_treated <- all_modules_non_treated[
    sapply(all_modules_non_treated, nrow) > 0
    ]

mt1_modules_by_cancer <- all_modules_non_treated %>%
    bind_rows() 

communities_tumor_type <- mt1_modules_by_cancer %>%
    group_by(tumor_type, community, signature) %>%
    mutate(
        n.appearances = n()
    ) %>%
    ungroup() %>%
    group_by(tumor_type, community) %>%
    mutate(
        n.samples = length(unique(sample))
    )



## Collapse to constellation plot

