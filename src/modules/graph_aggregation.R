library(igraph)
library(tidyverse)

moas <- readr::read_tsv(file = "reference/final_moas - Collapsed.tsv") %>%
    distinct(IDs, collapsed.MoAs, .keep_all = TRUE)

clinical_database <- data.table::fread(
    "results/annotation/clinical_metadata_v2_clean.tsv"
) %>%
    filter(
        tumor_subtype != "predicted_tumour" ## Get rid of non clear tumoral stuff
    )

module_list <- list.files(
    "results/modules",
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*_clusters.tsv"
)

## Load all modules
modules <- module_list %>%
    map(data.table::fread) %>%
    bind_rows(.id = "sample_study")

modules$sample_study <- module_list[as.numeric(modules$sample_study)]

modules_annotated <- modules %>%
    mutate(
        sample_study = stringr::str_remove(
            pattern = "results/modules/", string = sample_study)
    ) %>%
    separate(
        col = sample_study,
        into = c("study", "sample"),
        sep = "/"
    ) %>%
    mutate(
        sample = stringr::str_remove(string = sample, pattern = "_clusters.tsv")
    )

## Use right join to get rid of duplicate samples
modules_annotated_clinical <- modules_annotated %>%
    right_join(
        y = clinical_database,
        by = c("sample" = "sample", "study" = "study")
    ) %>%
    mutate(
        metagroup = case_when(
            treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_untreated",
            treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_untreated",
            treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_treated",
            treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_treated",
            study == "cell_lines_gabriella_kinker" ~ "cell_line",
            TRUE ~ "other"
        )
    )

## Perform aggregation of metagroup 1
modules_mt1 <- modules_annotated_clinical %>%
    filter(metagroup == "patient_primary_untreated")

cancers_to_consider <- modules_mt1 %>%
    group_by(tumor_type) %>%
    select(sample) %>%
    summarise(
        n.samples = length(unique(sample))
    ) %>%
    filter(
        n.samples >= 5
    ) %>%
    pull(tumor_type)

for(cancer in cancers_to_consider){
    print(paste0("Analyzing ", cancer))
    
    modules_mt_by_cancer <- modules_mt1 %>%
        filter(
            tumor_type == cancer
        ) %>%
        mutate(
            sample_bicluster = paste0(sample, "_", bicluster)
        )
    
    print("Calculating cutoff")
    dist_contribution <- modules_mt_by_cancer %>%
        group_by(study, sample, bicluster) %>%
        summarise(
            information_content = unique(information_content),
            n.drugs = n()
        ) %>%
        arrange(bicluster) %>%
        mutate(
            bicluster = as_factor(bicluster)
        ) %>%
        group_by(study, sample) %>%
        mutate(
            information_content = information_content / sum(information_content),
            information_content = round(information_content, digits = 2),
            cumulative_info = cumsum(information_content),
            sample_bicluster = paste0(sample, "_", bicluster)
        ) 
    
    print("Plotting distribution and cutoff")
    dist_plot_2 <- ggplot(data = dist_contribution) +
        geom_boxplot(aes(x = bicluster, y= cumulative_info)) +
        geom_smooth(aes(x = sort(as.numeric(bicluster)), y = cumulative_info)) +
        geom_hline(yintercept = 0.75) +
        scale_y_continuous(labels = scales::percent, n.breaks = 10) +
        ylab(label = "Relative contribution to total information content") +
        theme_classic() 
    
    ggsave(
        plot = dist_plot_2,
        filename = paste0("results/modules/annotated/", cancer, "_info_content.png"),
        height = 10,
        width = 10,
        dpi = 300
    )
    
    dist_contribution <- dist_contribution %>%
        filter(
            cumulative_info <= 0.75
        ) %>%
        pull(sample_bicluster)
    
    print("Generating edge data")
    module_edges <- modules_mt_by_cancer %>%
        filter(sample_bicluster %in% dist_contribution) %>%
        select(sample_bicluster, tumor_type) %>%
        distinct(.keep_all = FALSE) %>%
        as.data.frame()
    
    combinations <- combn(module_edges$sample_bicluster, m = 2, simplify = TRUE)
    combinations <- as.data.frame(t(combinations))
    
    print("Generating vertex data")
    relations <- data.frame(
        from = combinations$V1,
        to = combinations$V2
    ) %>%
        separate(col = "from", into = c("sample_1", "bicluster_1"), sep = "_") %>%
        separate(col = "to", into = c("sample_2", "bicluster_2"), sep = "_") %>%
        filter(sample_1 != sample_2) %>%
        rowwise() %>%
        mutate(
            positive_weight_intersect = length(
                intersect(
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_1 &
                            modules_mt_by_cancer$bicluster == bicluster_1 &
                            sign(modules_mt_by_cancer$cluster_contribution) == 1),
                        "signature"
                    ],
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_2 &
                            modules_mt_by_cancer$bicluster == bicluster_2 &
                            sign(modules_mt_by_cancer$cluster_contribution) == 1),
                        "signature"
                    ]
                )
                
            ),
            negative_weight_intersect = length(
                intersect(
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_1 &
                            modules_mt_by_cancer$bicluster == bicluster_1 &
                            sign(modules_mt_by_cancer$cluster_contribution) == -1),
                        "signature"
                    ],
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_2 &
                            modules_mt_by_cancer$bicluster == bicluster_2 &
                            sign(modules_mt_by_cancer$cluster_contribution) == -1),
                        "signature"
                    ]
                )
            ),
            positive_weight_union = length(
                dplyr::union(
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_1 &
                            modules_mt_by_cancer$bicluster == bicluster_1 &
                            sign(modules_mt_by_cancer$cluster_contribution) == 1),
                        "signature"
                    ],
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_2 &
                            modules_mt_by_cancer$bicluster == bicluster_2 &
                            sign(modules_mt_by_cancer$cluster_contribution) == 1),
                        "signature"
                    ]
                )
                
            ),
            negative_weight_union = length(
                dplyr::union(
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_1 &
                            modules_mt_by_cancer$bicluster == bicluster_1 &
                            sign(modules_mt_by_cancer$cluster_contribution) == -1),
                        "signature"
                    ],
                    modules_mt_by_cancer[(
                        modules_mt_by_cancer$sample == sample_2 &
                            modules_mt_by_cancer$bicluster == bicluster_2 &
                            sign(modules_mt_by_cancer$cluster_contribution) == -1),
                        "signature"
                    ]
                )
            )
            
            
        )
    
        print("Generating graph data")
        relations_filtered <- relations %>%
            mutate(
                positive_weight = positive_weight_intersect / positive_weight_union,
                negative_weight = negative_weight_intersect / negative_weight_union,
                weight = positive_weight + negative_weight
                ) %>%
            filter(weight > 0.25) %>% ## prune edges with 0 few conn.
            mutate(
                from = paste0(sample_1, "_", bicluster_1),
                to = paste0(sample_2, "_", bicluster_2)
            ) %>%
            select(from, to, weight)
        
        g <- graph_from_data_frame(
            relations_filtered,
            directed = FALSE,
            vertices=module_edges
        )
        
        print("Fast greedy classif")
        fc <- fastgreedy.community(
            graph = (g),
            weights = relations_filtered$weight
            )
        
        comms <- data.frame(edge=fc$names, community = fc$membership) %>%
            group_by(community) %>%
            mutate(
                n.members = n()
            ) %>%
            filter(
                n.members >= 5
            ) %>%
            left_join(
                y = modules_mt_by_cancer,
                by = c("edge" = "sample_bicluster")
            ) %>%
            left_join(
                y = moas[, c("IDs", "collapsed.MoAs")],
                by = c("signature" = "IDs")
            ) %>%
            select(
                edge,
                community,
                sample,
                bicluster,
                tumor_type,
                tumor_subtype, 
                signature,
                n.members,
                collapsed.MoAs
            )
        
        write.table(
            x = comms,
            file = paste0("results/modules/annotated/", cancer, "_communities.tsv"),
            sep = "\t",
            row.names = FALSE
            )
}
