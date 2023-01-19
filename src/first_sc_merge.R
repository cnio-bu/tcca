library("Seurat")

setwd("/local/sagarcia/bc-meta/single_cell/obj/")

brain_mets <- readRDS("brmets_hugo_gonzalez/all_malignant.rds")
full_brain <- merge(
    x = brain_mets[[1]],
    y = brain_mets[2:length(brain_mets)]
)

rm(brain_mets)
pancancer_dahlia <- readRDS("pancancer_dalia_barkley/all_malignant.rds")

merged_pancancer <- merge(x = full_brain, y = pancancer_dahlia)
rm(pancancer_dahlia)
gc()
saveRDS(object = merged_pancancer, file = "/local/sagarcia/test_merge.rds")

pancancer_sunny_wu <- readRDS(file = "pancancer_sunny_wu/all_malignant.rds")
merged_pancancer_v2 <- merge(x = merged_pancancer, y = pancancer_sunny_wu)

## test
merged_pancancer_v2 <- Seurat::FindVariableFeatures(
    merged_pancancer_v2,
    selection.method = "vst",
    nfeatures = 2000
    )

merged_pancancer_v2 <- Seurat::ScaleData(merged_pancancer_v2)
merged_pancancer_v2 <- RunPCA(
    merged_pancancer_v2,
    features = VariableFeatures(object = merged_pancancer_v2)
    )

my_elbow <- ElbowPlot(merged_pancancer_v2)


merged_pancancer_v2 <- FindNeighbors(merged_pancancer_v2, dims = 1:10)
merged_pancancer_v2 <- FindClusters(merged_pancancer_v2, resolution = 0.5)
merged_pancancer_v2 <- RunUMAP(object = merged_pancancer_v2, dims = 1:10)

my_umap <- DimPlot(merged_pancancer_v2, reduction = "umap")

library(ggplot2)

ggsave(plot = my_elbow, filename = "../../../elbow_test.png", dpi = 100)
ggsave(plot = my_umap, filename = "../../../clusters_umap.png", dpi = 100)

my_umap2 <- DimPlot(object = merged_pancancer_v2, group.by = "orig.ident")

ggsave(my_umap2, filename = "../../../studies_umap.png", dpi = 100, height = 10, width = 20)

first_clustering <- table(merged_pancancer_v2$orig.ident, merged_pancancer_v2$RNA_snn_res.0.5)
first_clustering <- as.data.frame(first_clustering)
colnames(first_clustering) <- c("sample", "cluster", "cells")

my_umap3 <- DimPlot(object = merged_pancancer_v2, group.by = "Cell_Type")

## test cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

merged_pancancer_v2 <- CellCycleScoring(
    object = merged_pancancer_v2,
    s.features = s.genes,
    g2m.features = g2m.genes
    )

my_umap4 <- DimPlot(object = merged_pancancer_v2, group.by = "Phase")

ggsave(
    filename = "../../../umap_cell_cycle.png",
    plot = my_umap4,
    height = 7,
    width = 7,
    dpi = 100
    )


merged_pancancer_v2 <- RunPCA(merged_pancancer_v2, features = c(s.genes, g2m.genes))
phase_pca_separation <- DimPlot(merged_pancancer_v2, reduction = "pca", group.by = "Phase")

ggsave(
    filename = "../../../pca_by_phase.png",
    phase_pca_separation,
    height = 7,
    width = 7,
    dpi = 100
)

merged_pancancer_v2 <- ScaleData(merged_pancancer_v2, 
                                 vars.to.regress = c("S.Score", "G2M.Score"),
                                 features = rownames(merged_pancancer_v2)
                                 )
