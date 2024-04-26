library(tidyverse)
library(BPCells)
library(Seurat)

## Set options
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

## Load full beyondcell mat
mat <- open_matrix_dir(dir = "results/beyondcell_bp/full_mat_beyondcell")
## Load full database
meta.data_full_clinical <- read_tsv(
    file = "results/beyondcell_bp/beyondcell_metadata_with_clinical.tsv"
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

seu[["RNA"]]$data <- mat

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
    nfeatures = 300
    )

seu <- ScaleData(seu, do.scale = TRUE, do.center = TRUE, scale.max = 20)
seu <- RunPCA(seu, npcs = 100, reduction.name = "pca")

seu <- FindNeighbors(
    object = seu,
    dims = 1:50,
    reduction = "pca"
    )

seu <- FindClusters(
    seu,
    resolution = 0.2,
    method = "igraph",
    algorithm = 2,
    random.seed = 120394,
    group.singletons = TRUE,
    cluster.name = "therapeutic_clusters_sketched_0.2"
    )

seu <- RunUMAP(seu, dims = 1:50, return.model = T, umap.method = "uwot")

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
write_matrix_dir(mat = sketched_mat, dir = "results/beyondcell_bp/sketch_mat_beyondcell")
write_matrix_dir(mat = sketched_mat_5k, dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k")

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
saveRDS(object = seu, file = "results/beyondcell_bp/beyondcell_pancancer.Rds")


## plot UMAP
tcs_umap <- DimPlot(
    object = seu,
    group.by = "therapeutic_clusters_0.2",
    reduction = "umap"
    )

saveRDS(object = tcs_umap, file = "results/beyondcell_bp/tcs_raw.rds")
