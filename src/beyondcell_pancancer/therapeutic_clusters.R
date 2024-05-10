library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

## Load full beyondcell mat
mat <- open_matrix_dir(dir = "results/beyondcell_bp_hq/full_mat_beyondcell")
## Load full database
meta.data_full_clinical <- read_tsv(
    file = "results/annotation/beyondcell_hq_metadata_with_clinical.tsv"
    )

meta.data_full_clinical <- meta.data_full_clinical %>%
    as.data.frame()

colnames(mat) <- paste0("c", c(1:ncol(mat)))

meta.data_full_clinical$new_cell_id <-  paste0("c", c(1:nrow(meta.data_full_clinical)))
rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "beyondcell_pancancer",
    meta.data = meta.data_full_clinical
)

## subset bc to remove bhupinder pal samples that are predicted tumours
seu <- subset(seu,
              subset = study != "brca_bhupinder_pal" | tumor_subtype != "predicted_tumour"
              )

## subset to remove T10 sample, as it has FEWER cells with bcscore than expr.
seu <- subset(seu, subset = sample != "T10")

## redo the index
seu$new_cell_id <- paste0("c", c(1:ncol(seu)))

## redo colnames
colnames(seu) <- seu$new_cell_id

seu[["RNA"]]$data <- seu[["RNA"]]$counts

## Use the full dataset as variable feat. We are talking 581. Let PCA summarise.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(mat)

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

seu <- FindVariableFeatures(
    object = seu,
    selection.method = "disp",
    nfeatures = 90
    )
seu <- ScaleData(seu, do.scale = TRUE, do.center = TRUE)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 90)

## QC
vimdizload <- VizDimLoadings(seu, dims = 1:2, reduction = "pca")

ggsave(
    filename = "results/tcca/clustering_expr/vim_loadings_pca.png",
    plot = vimdizload,
    height = 9,
    width = 16,
    dpi = 100
    )

seu_pca <- DimPlot(seu, reduction = "pca")

ggsave(
    filename = "results/tcca/clustering_expr/seu_pca_sketched_50.png",
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
    filename = "results/tcca/clustering_expr/loadings_lose_heats_33.png",
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
    filename = "results/tcca/clustering_expr/loadings_lose_heats_10.png",
    heats,
    width = 16,
    height = 9,
    dpi = 300
)

elbow <- ElbowPlot(seu, ndims = 90)
ggsave(
    filename = "results/tcca/clustering_expr/elbow_plot_seu.png",
    elbow,
    width = 16,
    height = 9,
    dpi = 100
)


seu <- FindNeighbors(
    object = seu,
    dims = 1:30,
    reduction = "pca"
    )

seu <- FindClusters(
    seu,
    resolution = 0.2,
    #method = "igraph",
    #algorithm = 2,
    random.seed = 120394,
    group.singletons = TRUE,
    cluster.name = "therapeutic_clusters_sketched_0.2"
    )

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
bc.markers <- FindAllMarkers(seu, only.pos = TRUE)
bc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

bc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

DoHeatmap(seu, features = top10$gene) + NoLegend()



seu <- Seurat::ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch_50k",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:50,
    refdata = list(therapeutic_clusters_0.2 = "therapeutic_clusters_sketched_0.2")
)

write_tsv(seu@meta.data, "results/annotation/tcs.tsv")

seu <- Seurat::SketchData(
    object = seu,
    assay = "RNA",
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch_5k"
    )

## export  sketches
sketched_mat <- seu[["sketch_50k"]]$data
sketched_mat_5k <- seu[["sketch_5k"]]$data
write_matrix_dir(mat = sketched_mat, dir = "results/beyondcell_bp/sketch_mat_beyondcell", overwrite = TRUE)
write_matrix_dir(mat = sketched_mat_5k, dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k", overwrite = TRUE)

## Calculate module scores for all samples for metacoms 1/2
DefaultAssay(seu) <- "RNA"

drugs_metacommunities_untreated <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
)

drugs_metacommunities_treated <- read.table(
    "results/modules/annotated/metagroup_patients_treated_consensus_drugs.tsv"
)

meta_coms_set_untreated <- split(
    drugs_metacommunities_untreated$signature,
    drugs_metacommunities_untreated$meta_community
)

meta_coms_set_treated <- split(
    drugs_metacommunities_treated$signature,
    drugs_metacommunities_treated$meta_community
)

seu <- AddModuleScore(
    object = seu,
    features = meta_coms_set_untreated,
    seed = 120394,
    slot = "data",
    name = "metacom_untreated_",
    ctrl = 20
)

seu <- AddModuleScore(
    object = seu,
    features = meta_coms_set_treated,
    seed = 120394,
    slot = "data",
    name = "metacom_treated_",
    ctrl = 20
)

## export object
saveRDS(object = seu, file = "results/beyondcell_bp_hq/beyondcell_pancancer.Rds")


## plot UMAP
tcs_umap <- DimPlot(
    object = seu,
    group.by = "therapeutic_clusters_0.2",
    reduction = "umap"
    )

saveRDS(object = tcs_umap, file = "results/beyondcell_bp/tcs_raw.rds")



## test
scaled.matrix <- t(apply(mat, 1, scales::rescale, to = c(0, 1)))
seu@assays$sketch_50k$scale.data <- scaled.matrix

seu <- RunPCA(object = seu, npcs = 90)

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

seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:10, umap.method = "uwot", n.components = 2)
