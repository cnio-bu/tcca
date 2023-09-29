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

## Collapse to constellation plot
extract_modules <- function(dt){
    dt <- dt %>%
        group_by(community) %>%
        mutate(
            n.samples = length(unique(sample)),
            n.edges = length(unique(edge))
        ) %>%
        group_by(community, signature) %>%
        reframe(
            n.samples = n.samples,
            n.edges = n.edges,
            n.appearances = n (),
            collapsed.MoAs = collapsed.MoAs
        ) %>%
        filter(
            n.appearances >= round(0.5 * n.samples, digits = 0)
        ) %>%
        distinct() 
    
    return(dt)
}

modules_fixed <- lapply(all_modules_non_treated, FUN = extract_modules)

c <- modules_fixed[[3]]

tes2 <- tes %>%
    group_by(community) %>%
    mutate(
        n.samples = length(unique(sample)),
        n.edges = length(unique(edge))
    ) %>%
    group_by(community, signature) %>%
    reframe(
        n.samples = n.samples,
        n.edges = n.edges,
        n.appearances = n (),
        collapsed.MoAs = collapsed.MoAs
    ) %>%
    filter(
        n.appearances >= round(0.5 * n.samples, digits = 0)
    ) %>%
    distinct() 

tes3 <- tes2 %>%
    group_by(community, collapsed.MoAs) %>%
    summarise(
        n.appearances = n()
    ) %>%
    arrange(community, desc(n.appearances))
