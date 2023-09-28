library(igraph)

## test
modules_mt1_breast <- modules_mt1 %>%
    filter(tumor_type == "PAAD" & bicluster <= 20)

## let's see
modules_edges <- modules_mt1_breast %>%
    select(sample, bicluster, tumor_subtype) %>%
    mutate(sample = paste0(sample, "_", bicluster)) %>%
    distinct(.keep_all = FALSE) %>%
    as.data.frame()

combinations <- combn(modules_edges$sample, m = 2, simplify = TRUE)
combinations <- as.data.frame(t(combinations))


relations <- data.frame(
    from = combinations$V1,
    to = combinations$V2
    ) %>%
    separate(col = "from", into = c("sample_1", "bicluster_1"), sep = "_") %>%
    separate(col = "to", into = c("sample_2", "bicluster_2"), sep = "_") %>%
    filter(sample_1 != sample_2) %>%
    rowwise() %>%
    mutate(
        weight_intersect = length(intersect(
            modules_mt1_breast[(modules_mt1_breast$sample == sample_1 & modules_mt1_breast$bicluster == bicluster_1), "signature"],
            modules_mt1_breast[(modules_mt1_breast$sample == sample_2 & modules_mt1_breast$bicluster == bicluster_2), "signature"]
        ))
    )

## filter low strength relationships
summary(relations$weight_intersect)

relations_filtered <- relations %>%
    filter(weight_intersect > 3) %>%
    mutate(
        from = paste0(sample_1, "_", bicluster_1),
        to = paste0(sample_2, "_", bicluster_2)
    ) %>%
    select(from, to, weight_intersect) %>%
    rename(
        "weight" = "weight_intersect"
    )
    


g <- graph_from_data_frame(
    relations_filtered,
    directed=FALSE,
    vertices=modules_edges
    ) 

#only works for undirected graphs, which this example is fine since symetric
fc <- fastgreedy.community((g), weights = relations_filtered$weight)

#make colors for different communities
colors_clusters <- membership(fc)

colors_dict <- c("red", "blue", "green", "pink", "orange", "brown")
names(colors_dict) <- c(1,2,3,4,5,6)
    
V(g)$color <- colors_dict[membership(fc)]
g = simplify(g)
plot(g, vertex.size = 4, vertex.label=NA)

com1 <- data.frame(edge=fc$names, community = fc$membership) %>%
    filter(
        community == 4
    )

comm1_drugs <- modules_mt1_breast %>%
    mutate(
        sample_cluster = paste0(sample, "_", bicluster)
    ) %>%
    filter(
        sample_cluster %in% com1$edge
    )


drugs_to_keep <- table(comm1_drugs$signature)
nsamples <- length(unique(modules_mt1_breast$sample))
drug_set <- names(drugs_to_keep[drugs_to_keep > 5])

## drug db
moas <- readr::read_tsv(file = "reference/final_moas - Collapsed.tsv") %>%
    filter(
        IDs %in% drug_set
    ) %>%
    group_by(collapsed.MoAs) %>%
    summarise(
        n.elements = n()
    )

