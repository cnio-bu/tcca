library(data.table)
library(igraph)
library(tidyverse)
library(multidplyr)

## setup parallel stuff
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")

moas <- readr::read_tsv(file = "reference/final_moas - Collapsed.tsv") %>%
    distinct(IDs, collapsed.MoAs, .keep_all = TRUE)

modules_annotated_clinical <- read_tsv(
    file = "results/annotation/modules_annotated_clinical.tsv"
    )

extract_communities <- function(group, n_samples){
    
    is_absent = !(group %in% unique(modules_annotated_clinical$metagroup))
    
    if (is_absent) {
        stop("Unknown metagroup")
    }
    
    ## Perform aggregation of metagroup 
    modules_mt1 <- modules_annotated_clinical %>%
        filter(metagroup == group)
    
    cancers_to_consider <- modules_mt1 %>%
        group_by(tumor_type) %>%
        select(sample) %>%
        summarise(
            n.samples = length(unique(sample))
        ) %>%
        filter(
            n.samples >= n_samples
        ) %>%
        pull(tumor_type)
    
    for (cancer in cancers_to_consider) {
        print(paste0("Analyzing ", cancer, " for ", group))
        
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
            geom_boxplot(aes(x = bicluster, y = cumulative_info)) +
            geom_smooth(aes(x = bicluster, y = cumulative_info)) +
           # geom_hline(yintercept = 0.75) +
            scale_y_continuous(labels = scales::percent, n.breaks = 10) +
            ylab(label = "Relative contribution to total information content") +
            theme_classic() 
        
        ggsave(
            plot = dist_plot_2,
            filename = paste0("results/modules/annotated/", group, "_", cancer, "_info_content.png"),
            height = 10,
            width = 10,
            dpi = 300
        )
        
        dist_contribution <- dist_contribution %>%
            filter(
                cumulative_info <= 1 # Do not perform this filtering anymore
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
        
        print("Copying module info to threads")
        cluster_copy(cluster, "modules_mt_by_cancer")
        
        print("Generating vertex data")
        relations <- data.table(
            from = combinations$V1,
            to = combinations$V2
        ) %>%
            separate(col = "from", into = c("sample_1", "bicluster_1"), sep = "_(?=[^_]+$)") %>%
            separate(col = "to", into = c("sample_2", "bicluster_2"), sep = "_(?=[^_]+$)") %>%
            filter(sample_1 != sample_2) %>%
            rowwise() %>%
            partition(cluster) %>%
            mutate(
                weight_intersect = nrow(
                    intersect(
                        modules_mt_by_cancer[(
                            modules_mt_by_cancer$sample == sample_1 &
                                modules_mt_by_cancer$bicluster == bicluster_1
                        ),
                             #   sign(modules_mt_by_cancer$cluster_contribution) == 1),
                            "signature"
                        ],
                        modules_mt_by_cancer[(
                            modules_mt_by_cancer$sample == sample_2 &
                                modules_mt_by_cancer$bicluster == bicluster_2
                        ),
                          #      sign(modules_mt_by_cancer$cluster_contribution) == 1),
                            "signature"
                        ]
                    )
                    
                ),
                weight_union = nrow(
                    dplyr::union(
                        modules_mt_by_cancer[(
                            modules_mt_by_cancer$sample == sample_1 &
                                modules_mt_by_cancer$bicluster == bicluster_1
                        ),
                                #sign(modules_mt_by_cancer$cluster_contribution) == 1),
                            "signature"
                        ],
                        modules_mt_by_cancer[(
                            modules_mt_by_cancer$sample == sample_2 &
                                modules_mt_by_cancer$bicluster == bicluster_2
                        ),
                                #sign(modules_mt_by_cancer$cluster_contribution) == 1),
                            "signature"
                        ]
                    )
                    
                )
                
            ) %>%
            collect()
        
        print("Generating graph data")
        relations_filtered <- relations %>%
            mutate(
                weight = weight_intersect /weight_union
            ) %>%
            filter(weight >= 0.25) %>% ## prune edges with 0 few conn.
            mutate(
                from = paste0(sample_1, "_", bicluster_1),
                to = paste0(sample_2, "_", bicluster_2)
            ) %>%
            select(from, to, weight)
        
    g <- graph_from_data_frame(
        relations_filtered,
        directed = FALSE,
        vertices = module_edges
    )
    
    print("Fast greedy classif")
    fc <- fastgreedy.community(
        graph = (g),
        weights = relations_filtered$weight
    )
    g = simplify(g)
    plot(g, vertex.size = 4, vertex.label=NA)
    
    png(
        filename = paste0(
            "results/modules/annotated/",
            group, "_",
            cancer,
            "_graph.png"
            ),
        height = 7,
        width = 7,
        units = "in",
        res = 100
        )
    
    plot(g, vertex.size = 4, vertex.label=NA)
    dev.off()

   
        comms <- data.frame(edge = fc$names, community = fc$membership) %>%
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
            x = relations_filtered,
            file = paste0("results/modules/annotated/", group, "_", cancer, "_relationships.tsv"),
            sep = "\t",
            row.names = FALSE
        )
        
        write.table(
            x = module_edges,
            file = paste0("results/modules/annotated/", group, "_", cancer, "_edges.tsv"),
            sep = "\t",
            row.names = FALSE
        )
        
        write.table(
            x = comms,
            file = paste0("results/modules/annotated/", group, "_", cancer, "_communities.tsv"),
            sep = "\t",
            row.names = FALSE
        )
        
    }
    
    
}

## For each metagroup, call the FUN
extract_communities(group = "patient_primary_untreated", n_samples = 5)
extract_communities(group = "patient_primary_treated", n_samples = 3)
extract_communities(group = "patient_metastatic_untreated", n_samples = 3)
extract_communities(group = "patient_metastatic_treated", n_samples = 3)
