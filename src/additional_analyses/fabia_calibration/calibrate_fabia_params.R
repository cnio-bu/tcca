library(beyondcell)
library(ComplexHeatmap)
library(fabia)
library(NbClust)
library(tidyverse)

pleural <- readRDS("results/pleural_rui_dong.rds")
pleural <- pleural[[1]]

norm_mat <- pleural@normalized

norm_mat_cor <- cor(t(norm_mat), method = "pearson")

hc <- hclust(as.dist(1 - norm_mat_cor), method = "complete")

## test blind clustering
png(
    filename = "results/figures/pleural_raw_similarity_scores.png",
    width = 14,
    height = 14,
    units = "in",
    res = 100
    )

blind_correlation_heat <- ComplexHeatmap::Heatmap(
    name = "Similarity score (pearson)",
    norm_mat_cor,
    cluster_rows = hc,
    cluster_columns = hc,
    show_column_names = FALSE,
    show_row_names = FALSE
    )

draw(blind_correlation_heat)
dev.off()

## Default fabia biclusters as done by bc-meta 2023
res <- fabia::fabias(
    X = norm_mat,
    p = 50, ## Hidden factors = biclusters. Max. to 50 (581 / 50 ~ 11 drugs)
    cyc = 500, ## iterations, keep it at 500
    spz = 0.5, ## minimum sparseness, Laplace.
    non_negative = 0, ## Allow negative factors a.k.a. negative vectors
    random = 1.0, ## allow random initialization of loadings,
    center = 2, ## median centering
    lap = 1, ## minimal value of the variational param
    nL = 1, ## do not allow drugs to be in > 1 biclust, it will report mirrowed biclusts otherwise,
    lL = 100, ## do not allow biclusters > 100 drugs. It tends to aggregate spurious clusters
)

## Judging from the heatmap, 50 biclusters is way too high
biclusters <- extractBic(fact = res, thresZ = 0.5)

## Get the most informative cluster
this_biclust <- biclusters$bic[1, ]
## THE FUCK 2 signatures? 



## use default bc-meta table of biclust generator
all_biclusters <- list()
all_cells <- list()
for(i in c(1:50)){
    this_biclust <- biclusters$bic[i, ]
    ## check length of the bicluster rowwise (aka drugs)
    if(length(this_biclust$bixn) <= 5) {
        next
    }else {
        named_list_sigs <- this_biclust$bixn
        sig_contributions <- this_biclust$bixv
        named_list_cells <- this_biclust$biypn
        cell_contributions <- this_biclust$biypv
        names(sig_contributions) <- named_list_sigs
        names(cell_contributions) <- named_list_cells
        all_biclusters[[i]] <- sig_contributions
        all_cells[[i]] <- cell_contributions
    }
}

bicluster_table <- enframe(all_biclusters) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "signature" = value_id
    )

table(bicluster_table$bicluster)

## Most biclusters are small... that is not what the hc suggest

## let's draw them
png("results/figures/fabia_aware_heat.png", width = 14, height = 14, units = "in", res = 100)
fabia_aware_correlation_heat <- ComplexHeatmap::Heatmap(
    name = "Similarity score (pearson)",
    norm_mat_cor[bicluster_table[order(bicluster_table$bicluster), ]$signature, 
                 bicluster_table[order(bicluster_table$bicluster), ]$signature
                 ],
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    column_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    row_split = bicluster_table$bicluster,
    column_split = bicluster_table$bicluster
)
draw(fabia_aware_correlation_heat)
dev.off()

## REDO fabia params
## Default fabia biclusters as done by bc-meta 2023
res_big <- fabia::fabias(
    X = norm_mat,
    p = 13, ## Hidden factors = biclusters. Max. to 50 (581 / 50 ~ 11 drugs)
    cyc = 50, ## iterations, keep it at 500
    spz = 0.5, ## minimum sparseness, Laplace.
    non_negative = 0, ## Allow negative factors a.k.a. negative vectors
    random = 1.0, ## allow random initialization of loadings,
    center = 0, ## median centering
    norm = 0,
    lap = 1, ## minimal value of the variational param
    nL = 1, ## do not allow drugs to be in > 1 biclust, it will report mirrowed biclusts otherwise,
    lL = 0, ## do not allow biclusters > 100 drugs. It tends to aggregate spurious clusters
)

## Judging from the heatmap, 50 biclusters is way too high
biclusters_big <- extractBic(fact = res_big, thresZ = 0.5)

## Get the most informative cluster
this_biclust <- biclusters_big$bic[1, ]
## 18 drugs now

## use default bc-meta table of biclust generator
all_biclusters <- list()
all_cells <- list()
for(i in c(1:13)){
    this_biclust <- biclusters_big$bic[i, ]
        named_list_sigs <- this_biclust$bixn
        sig_contributions <- this_biclust$bixv
        named_list_cells <- this_biclust$biypn
        cell_contributions <- this_biclust$biypv
        names(sig_contributions) <- named_list_sigs
        names(cell_contributions) <- named_list_cells
        all_biclusters[[i]] <- sig_contributions
        all_cells[[i]] <- cell_contributions
    
}

bicluster_table_big <- enframe(all_biclusters) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "signature" = value_id
    )

table(bicluster_table_big$bicluster)


## Biclusters are from 6 to 90 drugs
## let's draw them

## Try to remove drugs whose sign is opposite from the bicluster consensus
bicluster_table_big <- bicluster_table_big %>%
    mutate(
        sig_dir = sign(cluster_contribution)
    ) %>%
    group_by(bicluster) %>%
    mutate(
       dom_sig =  case_when(
           sum(sig_dir) >= 0 ~ 1,
           TRUE ~ -1
       )
    ) %>%
    filter(
        dom_sig == sig_dir
    )


png("results/figures/fabia_aware_heat_big.png", width = 14, height = 14, units = "in", res = 100)
fabia_aware_correlation_heat <- ComplexHeatmap::Heatmap(
    name = "Similarity score (pearson)",
    norm_mat_cor[bicluster_table_big[order(bicluster_table_big$bicluster), ]$signature, 
                 bicluster_table_big[order(bicluster_table_big$bicluster), ]$signature
                 ],
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = bicluster_table_big[order(bicluster_table_big$bicluster), ]$signature,
    column_order = bicluster_table_big[order(bicluster_table_big$bicluster), ]$signature,
    row_split = bicluster_table_big$bicluster,
    column_split = bicluster_table_big$bicluster
)
draw(fabia_aware_correlation_heat)
dev.off()



## blind double clustered bc
png("results/figures/blind_heat_bc.png", width = 14, height = 14, units = "in", res = 100)
blind_bc_clusts <- ComplexHeatmap::Heatmap(
    norm_mat,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cluster_row_slices = TRUE,
    cluster_column_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    show_column_names = FALSE,
    show_column_dend = FALSE
)
draw(blind_bc_clusts)
dev.off()

## bc aware of fabia clusts now
png("results/figures/arranged_heat_bc.png", width = 14, height = 14, units = "in", res = 100)
blind_bc_clusts <- ComplexHeatmap::Heatmap(
    norm_mat[bicluster_table_big$signature, ],
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    cluster_row_slices = TRUE,
    cluster_column_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    show_column_names = FALSE,
    show_column_dend = FALSE,
    row_split  = bicluster_table_big$bicluster
)
draw(blind_bc_clusts)
dev.off()
