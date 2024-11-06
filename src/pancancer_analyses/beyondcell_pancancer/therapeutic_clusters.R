library(tidyverse)
library(BPCells)
library(Seurat)
library(clustree)


## Set options
options(future.globals.maxSize = 20 * 1024^3)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell_immuno")

# #Load full beyondcell mat
mat <- open_matrix_dir(dir = "./full_mat_beyondcell")

# Load full database
meta.data_full_clinical <- read_tsv(
    file = "./beyondcell_metadata_with_clinical.tsv"
    )

meta.data_full_clinical <- meta.data_full_clinical %>%
    as.data.frame()

colnames(mat) <- paste0("c", c(1:ncol(mat)))

meta.data_full_clinical$new_cell_id <-  paste0("c", c(1:nrow(meta.data_full_clinical)))
rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

mat <- mat[1:589, ]

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "beyondcell_pancancer",
    meta.data = meta.data_full_clinical
)

# subset bc to remove bhupinder pal samples that are predicted tumours (not in lvl2 nor lvl1)
seu <- subset(seu,
              subset = study != "brca_bhupinder_pal" | tumor_subtype != "predicted_tumour"
)

# redo the index
seu$new_cell_id <- paste0("c", c(1:ncol(seu)))

## redo colnames
colnames(seu) <- seu$new_cell_id
rownames(seu@meta.data) <- seu$new_cell_id

seu[["RNA"]]$data <- seu[["RNA"]]$counts

# Check for the number of zeros for each cell
mat <- as.matrix(seu[["RNA"]]$counts)
seu@meta.data$zero_count <-  colSums(mat == 0)

# Remove cells with 0 for more than 10 % of the signatures
seu <- subset(seu, subset = zero_count < 0.1 * nrow(seu))

## Use the full dataset as variable feat. We are talking 581. Let PCA summarise.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu)

## Go for sketch, it's faster
seu <- SketchData(
    object = seu,
    ncells = 50000,
    method = "LeverageScore",
    over.write = TRUE,
    sketched.assay = "sketch_50k",
    seed = 120394
)

DefaultAssay(seu) <- "sketch_50k"

saveRDS(seu, "beyondcell_pancancer_sketch.rds")

seu <- readRDS("beyondcell_pancancer_sketch.rds")

seu <- FindVariableFeatures(
    object = seu,
    selection.method = "disp",
    nfeatures = 300
    )
seu <- ScaleData(seu, do.scale = TRUE, do.center = TRUE, scale.max = 20)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 100)

## QC
vimdizload <- VizDimLoadings(seu, dims = 1:2, reduction = "pca")

ggsave(
    filename = "results/plots_qc/vim_loadings_pca.png",
    plot = vimdizload,
    height = 9,
    width = 16,
    dpi = 100
    )

seu_pca <- DimPlot(seu, reduction = "pca")

ggsave(
    filename = "results/plots_qc/seu_pca_sketched_50.png",
    plot = seu_pca,
    height = 9,
    width = 16,
    dpi = 100
)

heats <- DimHeatmap(
    seu,
    dims = 1:33, 
    cells = 500,
    balanced = TRUE,
    slot = "scale.data",
    ncol = 3,
    nfeatures = 10,
    fast = FALSE,
    combine = TRUE
    )

ggsave(
    filename = "results/plots_qc/loadings_lose_heats_33.png",
    heats,
    width = 29,
    height = 16,
    dpi = 300
    )

heats <- DimHeatmap(
    seu,
    dims = 1:10, 
    cells = 500,
    balanced = TRUE,
    slot = "scale.data",
    ncol = 2,
    nfeatures = 10,
    fast = FALSE,
    combine = TRUE
)

ggsave(
    filename = "results/plots_qc/loadings_lose_heats_10.png",
    heats,
    width = 16,
    height = 9,
    dpi = 300
)

elbow <- ElbowPlot(seu, ndims = 148)
ggsave(
    filename = "results/plots_qc/elbow_plot_seu.png",
    elbow,
    width = 16,
    height = 9,
    dpi = 100
)

## Try top 10 clustering parameters based on Davies Boulding Index.
# metrics_tc <- readRDS("clustering_metrics.rds")
# names(metrics_tc) <- paste0("k.", c(20, 30, 50, 75, 100, 150, 200, 300))
# metrics_tc <- lapply(metrics_tc, function(X) {
#   names(X) <- paste0("res.", seq(.1,1,0.1))
#   return(X)
#   }
# )

# metrics_tc <- do.call(c, metrics_tc)
# metrics_tc <- lapply(metrics_tc, unlist)
# metrics_tc <- do.call(rbind, metrics_tc)
# params <- as.data.frame(metrics_tc) %>%
# arrange(dbi, desc(sil), desc(purity)) %>%
# slice_head(n = 10) %>% 
# row.names()

# for (x in params){
#     k.param <- unlist(strsplit(x, split = ".", fixed = TRUE))[2]
#     print(k.param)
#     res <- paste0(unlist(strsplit(x, split = ".", fixed = TRUE))[4:5], collapse = ".")
#     print(res)
#     seu <- FindNeighbors(
#                         object = seu,
#                         dims = 1:30,
#                         reduction = "pca",
#                         k.param = as.integer(k.param)
#     )
#     seu <- FindClusters(
#                         seu,
#                         resolution = as.numeric(res), 
#                         # method = "igraph",
#                         # algorithm = 2,
#                         random.seed = 120394,
#                         group.singletons = TRUE,
#                         cluster.name = paste0("therapeutic_clusters_sketched_k.", k.param, "_res.", res)
#     )
#     print(paste0("therapeutic_clusters_sketched_k.", k.param, "_res.", res))
# }

seu <- FindNeighbors(
    object = seu,
    dims = 1:30,
    reduction = "pca",
    k.param = 300
)
seu <- FindClusters(
    seu,
    resolution = 0.5, 
    # method = "igraph",
    # algorithm = 2,
    random.seed = 120394,
    group.singletons = TRUE,
    cluster.name = "therapeutic_clusters_sketched_k.300_res.0.5"
)


# # Evaluate different clustering resolutions with Silhoutte method.
# sampleDist <- stats::dist(t(seu[["sketch_50k"]]$counts), method = "euclidean")


# metadata_sketch <- seu@meta.data[!is.na(seu@meta.data$sketch_50k_snn_res.0.2),]
# clustree_tc <- clustree(metadata_sketch[, grep("sketch_50k_snn_res", colnames(metadata_sketch))],
#                         prefix = "sketch_50k_snn_res.")

# ggsave(
#     filename = "results/plots_qc/clustree_tcs.png",
#     clustree_tc,
#     width = 8,
#     height = 12,
#     dpi = 300
# )


seu <- RunUMAP(
    seu,
    dims = 1:30,
    return.model = T,
    umap.method = "uwot",
    #n.neighbors = 50,
    #metric = "correlation",
    #local.connectivity = 4
    )

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# Idents(seu) <- seu@meta.data$therapeutic_clusters_sketched_k.200_res.0.4
# bc.markers <- FindAllMarkers(seu, only.pos = TRUE)
# bc.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 1)

# bc.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 1) %>%
#     slice_head(n = 10) %>%
#     ungroup() -> top10

# top10_sigs <- DoHeatmap(seu, features = top10$gene) + NoLegend()
# ggsave(
#     filename = "results/plots_qc/DoHeatmap_top10_drugs.png",
#     top10_sigs,
#     width = 16,
#     height = 9,
#     dpi = 100
# )

print(colnames(seu@meta.data))
seu <- Seurat::ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch_50k",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:30,
    refdata = list(therapeutic_clusters_k.300.res.0.5 = "therapeutic_clusters_sketched_k.300_res.0.5"),
)

write_tsv(seu@meta.data, "results/tcs.tsv")

seu <- Seurat::SketchData(
    object = seu,
    assay = "RNA",
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch_5k"
    )

saveRDS(object = seu, file = "results/beyondcell_pancancer_final.Rds")

## export  sketches
sketched_mat <- seu[["sketch_50k"]]$data
sketched_mat_5k <- seu[["sketch_5k"]]$data
write_matrix_dir(mat = sketched_mat, dir = "results/sketch_mat_beyondcell", overwrite = TRUE)
write_matrix_dir(mat = sketched_mat_5k, dir = "results/sketch_mat_beyondcell_5k", overwrite = TRUE)


## plot UMAP
tcs_umap <- DimPlot(
    object = seu,
    group.by = "therapeutic_clusters_k.20.res.0.1",
    reduction = "umap"
    )

saveRDS(object = tcs_umap, file = "results/tcs_umap.rds")


# Check for the number of zeros for each cell
# mat <- as.matrix(seu[["RNA"]]$counts)
# seu@meta.data$zero_count <-  colSums(mat == 0)

# DefaultAssay(seu) <- "RNA"
# tcs_zero <- FeaturePlot(
#     object = seu,
#     features = "zero_count",
#     reduction = "full.umap"
#     ) +
#     ggtitle("") +
#     scale_shape_manual() +
#     xlab("UMAP1") +
#     ylab("UMAP2") +
#     theme(
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#     )

# tcs_zero$layers[[1]]$aes_params$size <- 0.1
# tcs_zero$layers[[1]]$aes_params$alpha <- 0.7

# ggsave(
#     plot = tcs_zero,
#     filename = "results/figures/tcs_zero_umap_full.png",
#     dpi = 300,
#     height = 7,
#     width = 7
#     )

# seu@meta.data$zero_count_true <- seu@meta.data$zero_count > 10

# zero_barplot <- ggplot(seu@meta.data, aes(x = therapeutic_clusters_k.200_res.0.4, fill = zero_count_true)) +
# geom_bar(position = "fill") +
# guides(fill = guide_legend(ncol = 1)) +
# theme_bw()

# ggsave("results/figures/zero10_barplot_tc.png", zero_barplot, width = 8, height = 8, dpi = 500)

## test
# scaled.matrix <- t(apply(mat, 1, scales::rescale, to = c(0, 1)))
# seu@assays$sketch_50k$scale.data <- scaled.matrix

# seu <- RunPCA(object = seu, npcs = 90)

# heats <- DimHeatmap(
#     seu,
#     dims = 1:10, 
#     cells = 500,
#     balanced = TRUE,
#     slot = "scale.data",
#     ncol = 2,
#     nfeatures = 10,
#     fast = FALSE,
#     combine = TRUE
# )

# seu <- FindNeighbors(seu)
# seu <- FindClusters(seu)
# seu <- RunUMAP(seu, dims = 1:10, umap.method = "uwot", n.components = 2)
