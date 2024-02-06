## Start again but now perform the analysis comm-wise instead of metacom-wise
## TADO: move this to their own script. This is awful
mmieloma_comms <- read.table(
    "results/modules/annotated/patient_primary_treated_MM_communities.tsv",
    header = TRUE
)

drugs_com <- mmieloma_comms %>%
    group_by(community, signature) %>%
    summarise(
        n.appearances = n()
    ) 

## Check the community to decide upon cutoff
dplot <- ggplot(data = drugs_com,
                aes(x  = n.appearances,
                    colour = as.factor(community)
                )
) +
    geom_line(aes(y = 1 - ..y..), stat = "ecdf") +
    scale_x_continuous(
        name = "",
        limits = c(1, max(drugs_com$n.appearances)),
        n.breaks = max(drugs_com$n.appearances)
    ) +
    scale_y_continuous(
        name = "Proportion of drugs above threshold",
        labels = scales::percent_format()
    ) +
    scale_colour_discrete(name = "Community") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
    )

ggsave(
    filename = "results/figures/mmieloma_communities_by_drug.png",
    plot = dplot,
    dpi = 100
)

drugs_com_filtered <- drugs_com %>%
    filter(n.appearances >= 4) %>%
    arrange(community, desc(n.appearances)) %>%
    left_join(y = drugs[,c("IDs", "preferred.drug.names", "collapsed.MoAs")],
              by = c("signature" = "IDs")
    )

coms_set <- split(
    drugs_com_filtered$signature,
    drugs_com_filtered$community
)

bc <- readRDS("results/")
bc <- AddModuleScore(
    bc,
    features = coms_set,
    slot = "data",
    seed = 120394,
    ctrl = 10,
    name = "com_"
)

module_mat <- bc@meta.data[colnames(bc@assays$sketch_10k_new$counts), ] %>%
    rownames_to_column("cell_barcode") %>%
    select(cell_barcode, com_1:com_16, treatment_group) %>%
    as.data.frame()

write.table(
    x = module_mat,
    file = "results/mmieloma/community_mat.tsv",
    sep = "\t"
)