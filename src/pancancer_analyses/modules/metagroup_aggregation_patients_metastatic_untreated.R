library(igraph)
library(scales)
library(tidyverse)

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
            n.appearances = n(),
            collapsed.MoAs = collapsed.MoAs,
            tumor_type = tumor_type
        ) %>%
        filter(
            n.appearances >= round(0.2 * n.edges, digits = 0)
        ) %>%
        distinct()
    
    return(dt)
}


## Patients primary treated
all_modules_treated <- list.files(
    path = "results/modules/annotated",
    pattern = "patient_metastatic_untreated_.*_communities\\.tsv",
    full.names = TRUE
)


moas <- readr::read_tsv(file = "reference/final_moas - Collapsed.tsv") %>%
    distinct(IDs, collapsed.MoAs, .keep_all = TRUE)

all_modules_treated <- all_modules_treated %>%
    map(read_tsv, id = "module_source")

## Check if no com
all_modules_treated <- all_modules_treated[
    sapply(all_modules_treated, nrow) > 0
]

mt1_modules_by_cancer <- all_modules_treated %>%
    bind_rows()

## Prepare communities for aggregation
modules_fixed <- lapply(
    all_modules_treated,
    FUN = extract_modules
)

modules_agg <- modules_fixed %>%
    bind_rows() %>%
    mutate(
        cancer_edge = paste0(tumor_type, "_", community)
    ) %>%
    as.data.frame()

## Generate edge data
module_edges <- modules_agg %>%
    select(cancer_edge, tumor_type) %>%
    distinct(.keep_all = FALSE) %>%
    as.data.frame()

combinations <- combn(module_edges$cancer_edge, m = 2, simplify = TRUE)
combinations <- as.data.frame(t(combinations))

## Generate vertex data
relations <- data.frame(
    from = combinations$V1,
    to = combinations$V2
) %>%
    separate(col = "from", into = c("tumor_1", "community_1"), sep = "_", remove = FALSE) %>%
    separate(col = "to", into = c("tumor_2", "community_2"), sep = "_", remove = FALSE) %>%
    filter(tumor_1 != tumor_2) %>%
    rowwise() %>%
    mutate(
        positive_weight_intersect = length(
            intersect(
                modules_agg[modules_agg$cancer_edge == from,
                            "signature"
                ],
                modules_agg[
                    modules_agg$cancer_edge == to,
                    "signature"
                ]
            )
        ),
        positive_weight_union = length(
            dplyr::union(
                modules_agg[modules_agg$cancer_edge == from,
                            "signature"
                ],
                modules_agg[
                    modules_agg$cancer_edge == to,
                    "signature"
                ]
            )
            
        )
    )


relations_filtered <- relations %>%
    mutate(
        weight = positive_weight_intersect / positive_weight_union
    ) %>%
    filter(weight > 0.1) %>% ## prune edges with 0 few conn.
    select(from, to, weight)

g <- graph_from_data_frame(
    relations_filtered,
    directed = FALSE,
    vertices=module_edges
)

fc <- fastgreedy.community(
    graph = (g),
    weights = relations_filtered$weight
)


piti = c("#6cca8e","#8398dc","#ea95ae","#1dade6", "#ff5f76", "#ffb6b6","#fff154","#ba7fff","#ffdd56", "#4b71e5",# "#cccccc",
                  "#ff6600","#add82f","#ff3333","#0dba3c", "#ff864c", "#c4ea94","#666699","#888888","#b8c0ba", "#d58aca","#6da753","#ca9a8c","#ff4430","#e06d23")
                  
names(piti) <- c(1:24)

V(g)$color <- piti[membership(fc)]
g = simplify(g)
Isolated = which(degree(g) == 0)
G2 = delete.vertices(g, Isolated)
plot(G2, vertex.size = 4)

plot(g, vertex.size = 4)

comms <- data.frame(edge = fc$names, meta_community = fc$membership) %>%
    group_by(meta_community) %>%
    mutate(
        n.cancer_edges = n()
    ) %>%
    filter(
        n.cancer_edges >= 5
    ) %>%
    left_join(
        y = modules_agg,
        by = c("edge" = "cancer_edge")
    )  %>%
    left_join(
        y = moas[, c("IDs", "preferred.drug.names")],
        by = c("signature" = "IDs")
    ) %>%
    select(
        edge,
        community,
        meta_community,
        tumor_type,
        signature,
        n.cancer_edges,
        collapsed.MoAs,
        preferred.drug.names
    )

write.table(
    x = comms,
    file = paste0("results/modules/annotated/patients_metastatic_untreated_meta_groups.tsv"),
    sep = "\t",
    row.names = FALSE
)


comms_agg <- comms %>%
    group_by(meta_community, collapsed.MoAs) %>%
    mutate(n_sigs_in_moa = n()) %>%
    ungroup() %>%
    group_by(meta_community) %>%
    mutate(
        meta_size = n(),
        n_sigs_moa_by_com = round(n_sigs_in_moa / meta_size, digits = 2)
    ) %>%
    select(meta_community, n_sigs_in_moa, meta_size, collapsed.MoAs) %>%
    distinct() %>%
    arrange(meta_community, desc(n_sigs_in_moa))

drugs_com <- comms %>%
    group_by(meta_community, signature) %>%
    summarise(
        n.appearances = n()
    ) %>%
    filter(n.appearances >= 3) %>%
    arrange(meta_community, desc(n.appearances)) %>%
    left_join(y = moas[,c("IDs", "preferred.drug.names", "collapsed.MoAs")],
              by = c("signature" = "IDs")
    )

## Check the community to decide upon cutoff

dplot <- ggplot(data = drugs_com,
                aes(x  = n.appearances,
                    colour = as.factor(meta_community)
                )
) +
    geom_line(aes(y = 1 - ..y..), stat = "ecdf") +
    scale_x_continuous(
        name = "N. appearances for drug",
        limits = c(1, max(drugs_com$n.appearances)),
        n.breaks = 10
    ) +
    scale_y_continuous(
        name = "Proportion of drugs above threshold",
        labels = scales::percent_format()
    ) +
    scale_colour_discrete(name = "Meta community") +
    theme_bw()

ggsave(
    plot = dplot,
    filename = "results/figures/proportion_drugset_dist_metastatic_untreated.png",
    dpi = 100,
    height = 10,
    width = 10
)

write.table(
    x = drugs_com,
    file = "results/modules/annotated/metagroup_patients_metastatic_untreated_consensus_drugs.tsv"
)
