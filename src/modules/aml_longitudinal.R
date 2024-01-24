library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)


mat <- open_matrix_dir(dir = "results/aml/aml_beyondcell_mat")
meta.data <- read.table("results/aml/bc_meta.tsv")

bc <- CreateSeuratObject(counts = mat, assay = "RNA", meta.data = meta.data)

## Add expression related features
seu <- readRDS("results/aml/seu.Rds")
study_metadata <- seu@meta.data

study_metadata <- study_metadata %>%
    rownames_to_column("cell_id")

rm(seu)
gc()

meta.data_annotated <- meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(y = study_metadata, by = "cell_id")

rownames(meta.data_annotated) <- meta.data_annotated$cell_id
bc@meta.data <- meta.data_annotated

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

bc@assays$RNA$data <- bc@assays$RNA$counts
bc <- ScaleData(object = bc, do.scale = TRUE, do.center = TRUE)

bc <- AddModuleScore(
    object = bc,
    features = meta_coms_set_untreated,
    seed = 120394,
    slot = "data",
    name = "metacom_untreated_",
    ctrl = 20
    )

bc <- AddModuleScore(
    object = bc,
    features = meta_coms_set_treated,
    seed = 120394,
    slot = "data",
    name = "metacom_treated_",
    ctrl = 20
)

## load drug data
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    distinct() %>%
    as.data.frame()

rownames(drugs) <- drugs$IDs

## Setting up flag for variable features
bc <- FindVariableFeatures(bc)
VariableFeatures(bc) <- rownames(bc)

## All patients sketch by metacom variability
bc <- SketchData(
    object = bc,
    assay = "RNA",
    ncells = 20000,
    sketched.assay = "sketch_10k",
    method = "LeverageScore",
    var.name = "leverage.score",
    seed = 120394,
    verbose = TRUE,
    over.write = TRUE
)

DefaultAssay(bc) <- "sketch_10k"

## Export sketched mats for each treatment group to draw a heatmap
module_mat_sketch <- bc@meta.data[colnames(bc@assays$sketch_10k$data), ] %>%
    select(metacom_untreated_1:metacom_treated_6) %>%
    as.data.frame()


write.table(module_mat_sketch, file = "results/aml/metacom_mat_20k_sketch.tsv")

## Export full mat for each treatment group for other analyses
module_mat <- bc@meta.data[colnames(bc@assays$RNA$data), ] %>%
    select(metacom_untreated_1:metacom_treated_6) %>%
    as.data.frame()

write.table(x = module_mat, file = "results/aml/metacom_mat_full.tsv")
