library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)

memory_mat <- readRDS("results/paad/bc_normalized_mat_rgrss.rds")
meta.data <- read.table("results/paad/seu_metadata_malignants.tsv")

seu <- CreateSeuratObject(
    counts = memory_mat,
    assay = "RNA",
    meta.data = meta.data
    )

seu@assays$RNA$data <- seu@assays$RNA$counts
seu <- ScaleData(object = seu, do.scale = TRUE, do.center = TRUE)

## Calculate metacom enrichment
drugs_metacommunities_untreated <- read.table(
    "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
)

meta_coms_set_untreated <- split(
    drugs_metacommunities_untreated$signature,
    drugs_metacommunities_untreated$meta_community
)

seu <- AddModuleScore(
    object = seu,
    features = meta_coms_set_untreated,
    seed = 120394,
    slot = "data",
    name = "metacom_untreated_",
    ctrl = 20
)

## Setting up flag for variable features
seu <- FindVariableFeatures(seu)
VariableFeatures(seu) <- rownames(seu)

## All patients sketch by metacom variability
seu <- SketchData(
    object = seu,
    assay = "RNA",
    ncells = 2000,
    sketched.assay = "sketch_200",
    method = "LeverageScore",
    var.name = "leverage.score",
    seed = 120394,
    verbose = TRUE,
    over.write = TRUE
)

DefaultAssay(seu) <- "sketch_200"

## Export sketched mats for each treatment group to draw a heatmap
module_mat_sketch <- seu@meta.data[colnames(seu@assays$sketch_200$data), ] %>%
    select(metacom_untreated_1:metacom_untreated_6) %>%
    as.data.frame()

write.table(module_mat_sketch, file = "results/paad/metacom_mat_2k_sketch.tsv")

## Export full mat for each treatment group for other analyses
module_mat <- seu@meta.data[colnames(seu@assays$RNA$data), ] %>%
    select(metacom_untreated_1:metacom_untreated_6) %>%
    as.data.frame()

write.table(x = module_mat, file = "results/paad/metacom_mat_full.tsv")

## switch back to full dataset for clustering. No need for sketchs with 16k cells
DefaultAssay(seu) <- "RNA"

seu <- RunPCA(object = seu, npcs = 50)
seu <- FindNeighbors(object = seu)
seu <- FindClusters(seu, resolution = 0.2, cluster.name = "therapeutic_clusters")
seu <- RunUMAP(object = seu, dims = 1:50)

saveRDS(object = seu, file = "results/paad/seu_beyondcell_with_tcs.rds")
