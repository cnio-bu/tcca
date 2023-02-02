library("beyondcell")
library("patchwork")
library("tidyverse")

setwd("/local/sagarcia/bc-meta/seurat")
bc <- readRDS("dalia_merged_bc.rds")

bc@data[is.nan(bc@data)] <- 0
bc@normalized[is.nan(bc@normalized)] <- 0
bc@scaled[is.nan(bc@scaled)] <- 0

bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))
bc <- bcUMAP(bc = bc, pc = 10, k.neighbors = 20, res = 0.05)


saveRDS(object = bc, file = "dalia_merged_bc_rgrss.rds")

test <- Seurat::ScaleData(
    object = mat,
    features = rownames(mat),
    do.scale = TRUE,
    do.center = TRUE,
)

sv <- svd(t(test))

U <- sv$u
V <- sv$v
D <- sv$d

## U is an un-scaled PC, Z is scaled PC
Z <- t(mat) %*% V

## let's see our clusters
clusters <- bcClusters(bc = bc, idents = "bc_clusters_res.0.1")

## pca of tcs
bcpca <- prcomp(x = bc@normalized)

bcmeta <- bc@meta.data
bcmeta$merged_codes <- rownames(bcmeta)

pcarot <- bcpca$rotation %>%
    as_tibble(rownames = "cell") %>%
    select(cell, PC1, PC2, PC3, PC4) %>%
    left_join(
        y = bcmeta[, c("merged_codes", "bc_clusters_res.0.1", "cancer")],
        by = c("cell" = "merged_codes")
    )


bc_pca <- ggplot(data = pcarot, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = bc_clusters_res.0.1)) +
    ggtitle("Top row: therapeutic clusters") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

bc_pca_coop <- ggplot(data = pcarot, aes(x = PC2, y = PC3)) +
    geom_point(aes(color = bc_clusters_res.0.1)) + 
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

bc_pca_can <- ggplot(data = pcarot, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = cancer)) +
    ggtitle("Bottom row: cancer types") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

bc_pca_coop_can <- ggplot(data = pcarot, aes(x = PC2, y = PC3)) +
    geom_point(aes(color = cancer)) + 
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

third_dim_pca <- bc_pca + bc_pca_coop + bc_pca_can + bc_pca_coop_can +
    plot_layout(guides = "collect", ncol = 2) +
    plot_annotation(tag_levels = "A") 



first_axis_loadings <- V[, 1]
names(first_axis_loadings) <- rownames(mat)

first_axis_loadings_tib <- tibble::enframe(
    first_axis_loadings,
    name = "drug", 
    value = "loading"
) %>%
    arrange(desc(abs(loading)))


axis_distribution <- first_axis_loadings_tib %>%
    ggplot(aes(x = loading)) +
    geom_histogram(col = "white", bins = 100) + 
    geom_density()  +
    ggtitle("Distribution of the loadings for PC1") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

cutoff <- 0.05

## Get rid of the abs stuff, focus on the pos. tail for cluster id.
top_axis_candidates <- first_axis_loadings_tib %>%
    filter(loading > cutoff)

## Annotate the candidates
## get top level moa annot
moas <- drugInfo$MoAs


top_axis_candidates_annotated <- top_axis_candidates %>%
    left_join(y = moas[, c("IDs", "main.MoAs")], by = c("drug" = "IDs")) %>%
    left_join(
        y = drugInfo$IDs[, c("IDs", "preferred.drug.names", "studies")],
        by = c("drug" = "IDs")) %>%
    distinct(.keep_all = TRUE)

top_moas <- top_axis_candidates_annotated %>%
    group_by(main.MoAs) %>%
    summarise(appearances = n())

## let's assume that the main MOA is MAPK inhibitor driven :)

## plot the top drug
top_drug <- top_axis_candidates_annotated %>%
    arrange(desc(abs(loading))) %>%
    head(n = 1) %>%
    pull(drug)


## show me the tcs and cancers
tcs_cancers <- bcClusters(bc = bc, idents = "cancer", UMAP = "beyondcell")


## let's see some diff. enrichment in the module
module_drugs <- top_axis_candidates_annotated$drug

norm_bcs <- t(bc@normalized[module_drugs, ]) %>%
    as_tibble(rownames = "cell") %>%
    pivot_longer(
        cols = all_of(module_drugs),
        names_to = "signature",
        values_to = "score"
    ) %>%
    left_join(
        y = bcmeta[, c("merged_codes", "bc_clusters_res.0.1", "cancer")],
        by = c("cell" = "merged_codes")
    ) %>%
    group_by(signature, bc_clusters_res.0.1) %>%
    summarise(
        average_cell_enrichment = median(score)
    ) 


module_box <- ggplot(data = norm_bcs, aes(x = bc_clusters_res.0.1, y = average_cell_enrichment)) +
    geom_point(alpha = 0.3) +
    geom_boxplot(aes(fill = bc_clusters_res.0.1)) +
    labs(x = "", y = "Average cell-wise enrichment by drug") +
    scale_fill_discrete(name = "Beyondcell therapeutic clusters") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    )

ggsave(plot = module_box, filename = "new_approach/boxplot_bc_by_module1.png")
ggsave(
    plot = third_dim_pca,
    filename = "new_approach/multipca.png",
    height = 14,
    width = 14
)

ggsave(plot = clusters, filename = "new_approach/bc_clusters_tcs.png")
ggsave(plot = tcs_cancers, filename = "new_approach/bc_tcs_cancers.png")
ggsave(plot = axis_distribution, filename = "new_approach/axis_distribution_pc1.png")

saveRDS(object = bc, "new_approach/bc_dalia_rgrss.rds")
saveRDS(
    object = top_axis_candidates_annotated,
    file = "new_approach/top_axis1_candidates_annotated.rds"
)

saveRDS(object = V, file = "new_approach/V_mat.rds")

## now, go for a simple  NMF
library("RcppML")

## non negative mat from bcnorm
non_negative_mat <- bc@normalized
non_negative_mat[non_negative_mat < 0] <- 0

model <- RcppML::nmf(non_negative_mat, 10, verbose = T, seed = 1234)

## again, decompose
w <- model$w
d <- model$d
h <- model$h

## amplitude matrix is now W
dim(w)

rownames(w) <- rownames(non_negative_mat)
colnames(w) <- paste0("component", 1:10)

rownames(h) <- paste0("component", 1:10)
colnames(h) <- colnames(non_negative_mat)

## same as BC clusters
kmeans_NMF_res <- kmeans(t(h), centers = 13)
kmeans_NMF_clusters <- kmeans_NMF_res$cluster

bc@meta.data$kmeans <- kmeans_NMF_clusters[rownames(bc@meta.data)]

bcClusters(bc = bc, idents = "kmeans")

## component vs clusters
nmf_df <- t(h) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cell") %>%
    left_join(bc@meta.data, by = c("cell" = "merged_barcodes"))

kmeans_components <- ggplot(
    data = nmf_df, aes(x = bc_clusters_res.0.1, y = component10)) +
    geom_point(aes(color = bc_clusters_res.0.1)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    xlab("") + 
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()
    ) 

## get top drugs by nmf
drugs_nmf <- w %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "drug") %>%
    arrange(desc(component10)) 


library(pheatmap)
pheatmap::pheatmap(
    model$w,
    show_rownames = FALSE,
    show_colnames = FALSE
)

library("fastICA")

norm_mat <- bc@normalized

ica_res <- fastICA(X = t(norm_mat), n.comp = 10, alg.typ = "deflation")

## Pattern matrix is now S
dim(ica_res$S)


##
## The data matrix X is considered to be a linear combination of non-Gaussian (independent) components
## i.e. X = SA where columns of S contain the independent components and A is a linear mixing matrix.
## In short ICA attempts to ‘un-mix’ the data by estimating an un-mixing matrix W where XW = S.
##
##

ica_df <- ica_res$S %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cell") %>%
    left_join(
        y = bcmeta[,c("bc_clusters_res.0.1", "merged_codes", "cancer")],
        by = c("cell" = "merged_codes")
    ) 

## I tricked and plot many of them, and found component 10 is associated with B cells
ica_distribution <- ggplot(ica_df, aes(x = bc_clusters_res.0.1, y = V2)) +
    geom_boxplot(aes(fill = bc_clusters_res.0.1)) +
    theme_classic(base_size = 14) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## ICA cell 2 -> let's see weights

ica_v2_drugs <- tibble(
    drug = rownames(bc@normalized),
    weights = ica_res$A[2, ]
) %>%
    arrange(desc(weights))


## ICA V2, cluster 3
pheatmap::pheatmap(ica_res$A, name = "Amplitude matrix ICA focusing on V2, (TC3)")
