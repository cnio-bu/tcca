library(BPCells)
library(Seurat)
library(tidyverse)

## set wd, required for bpcells relative path alloc.
setwd("results")

## Tell Seurat to work with on disk storage
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

seu <- readRDS(
    file = "seurat/v5/full_bc_meta_merged/full_bc_meta_merged.Rds"
    )


## Import the database
clinical_metadata <- data.table::fread(
    "annotation/clinical_metadata_v2_clean.tsv"
    )

seurat_meta <- seu@meta.data %>%
    left_join(
        y = clinical_metadata,
        by = c("sample")
    )


## Reannotate seurat object

## Test normalize
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)


## Näive clustering
seu <- SketchData(
    object = seu,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
    )

DefaultAssay(seu) <- "sketch"

seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu) 
seu <- FindNeighbors(seu, dims=1:50)
seu <- FindClusters(seu, resolution = 2)
seu <- RunUMAP(seu, dims=1:50, return.model = TRUE) 

sketched_umap <- DimPlot(seu, label = T, label.size = 3, reduction = "umap") + 
    NoLegend()

## Extend results back to full
seu <- ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims=1:50,
    refdata = list(cluster_full = "seurat_clusters")
)

DefaultAssay(seu) <- "RNA"

full_umap <- DimPlot(
    object = seu,
    label = T,
    label.size = 3,
    reduction = "full.umap",
    group.by = "cluster_full",
    alpha = 0.1
    ) + NoLegend()

ggsave(
    filename = "seurat/full_umap_sketched_default.png",
    plot = full_umap
)



reduction_umap <- as.data.frame(seu[["full.umap"]]@cell.embeddings)

test_tumor_intersection <- intersect(seu@meta.data$sample, clinical_metadata$sample)

seu <- subset(x = seu, subset = sample %in% test_tumor_intersection)

samples_clinical <- seu@meta.data %>%
    select(-cell, -c(17:36)) %>%
    rownames_to_column("cell") %>%
    left_join(
        y = clinical_metadata,
        by = "sample"
    ) %>%
    mutate(
        refined_tumor_site = case_when(
            refined_tumor_site == "" ~ "Unknown",
            TRUE ~ refined_tumor_site
        )
    ) %>%
    distinct(cell, .keep_all = TRUE)

samples_clinical <- as.data.frame(samples_clinical)
rownames(samples_clinical) <- samples_clinical$cell

seu@meta.data <- samples_clinical

full_umap_histology <- DimPlot(
    object = seu,
    label = T,
    label.size = 3,
    reduction = "full.umap",
    group.by = "refined_tumor_site",
    alpha = 0.1,
    repel = TRUE,
    raster = FALSE,
    pt.size = 2
) 

ggsave(
    plot = full_umap_histology,
    filename = "seurat/full_umap_refined_tumor_site.png",
    dpi = 300,
    height = 14,
    width = 28
)


## Go back to sketch
DefaultAssay(seu) <- "sketch"

sketch_umap_histology <- DimPlot(
    object = seu,
    label = T,
    label.size = 3,
    reduction = "umap",
    group.by = "refined_tumor_site",
    alpha = 0.1,
    repel = TRUE,
    raster = FALSE,
) + NoLegend()

ggsave(
    plot = sketch_umap_histology,
    filename = "seurat/sketch_umap_refined_tumor_site.png",
    dpi = 100,
    height = 7,
    width = 7
)

`sketch_umap_primary <- DimPlot(
    object = seu,
    label = T,
    label.size = 3,
    reduction = "umap",
    group.by = "sample_type",
    alpha = 0.1,
    repel = TRUE,
    raster = FALSE,
) 

ggsave(
    plot = sketch_umap_primary,
    filename = "seurat/sketch_umap_sample_type.png",
    dpi = 100,
    height = 7,
    width = 7
)

## is cell_line
seu$is_cell_line <- seu$study == "cell_lines_gabriella_kinker"
seu@meta.data$is_cell_line <- as.factor(seu@meta.data$is_cell_line)
levels(seu@meta.data$is_cell_line) <- c("No", "Yes")

sketch_umap_cell_line <- DimPlot(
    object = seu,
    label = T,
    label.size = 3,
    reduction = "umap",
    group.by = "is_cell_line",
    alpha = 0.1,
    repel = TRUE,
    raster = FALSE,
) 

ggsave(
    plot = sketch_umap_cell_line,
    filename = "seurat/sketch_umap_cell_line.png",
    dpi = 100,
    height = 7,
    width = 7
)

## Integration test
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- FindVariableFeatures(seu, verbose = FALSE)

seu <- SketchData(
    seu,
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch"
    )

