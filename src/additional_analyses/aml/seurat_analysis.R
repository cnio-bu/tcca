library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)

sample_dirs <- list.dirs(
    path = "raw/aml_sander_lambo",
    recursive = FALSE
)

samples <- basename(sample_dirs)

## load metadata
meta.files <- list.files(
    path = "raw/aml_sander_lambo",
    pattern = "^GSM.*.tsv",
    recursive = FALSE,
    full.names = TRUE
)

meta.data <- meta.files %>%
    map(read_tsv) %>%
    bind_rows()


## add geo ids
geos <- data.table::fread(geo_to_samples)
seu_list <- list()  

for (i in c(1:length(samples))){
    mat_dir <- sample_dirs[i]
    sample_name <- samples[i]
    
    mat <- Seurat::Read10X(data.dir = mat_dir, gene.column = 1)
    mat = as.sparse(mat)
    mat_dir = paste0("results/aml/", sample_name)
    write_matrix_dir(mat = mat, dir = mat_dir, overwrite = TRUE)
    mat = open_matrix_dir(dir = mat_dir)
    
    sample_meta <- meta.data %>%
        filter(gsm_id == sample_name) %>%
        as.data.frame()
    
    rownames(sample_meta) <- sample_meta$Cell_Barcode
    
    seu <- CreateSeuratObject(
        counts = mat,
        project = "aml_sander",
        assay = "RNA",
        meta.data = sample_meta
    )
    seu_list <- c(seu_list, seu)
}

seu <- merge(seu_list[[1]], seu_list[2:75])

rm(seu_list)
gc()

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(object = seu)
seu <- ScaleData(seu)

saveRDS(
    object = seu,
    file = "results/aml/seu.Rds",
    destdir = "results/aml", overwrite = TRUE
    )

seu <- subset(seu, subset = GEO_ID != "PAWNPU_Dx_scRNA")

seu <- SketchData(
    seu,
    ncells = 500,
    method = "LeverageScore",
    var.name = "leverage.score",
    sketched.assay = "sketch_500",
    over.write = TRUE,
    seed = 120394
    )

saveRDS(
    object = seu,
    file = "results/aml/seu.Rds",
)

DefaultAssay(seu) <- "sketch_500"

seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

# integrate the datasets
seu <- IntegrateLayers(
    seu,
    method = HarmonyIntegration,
    new.reduction = "integrated.harmony"
)

# cluster the integrated data
seu <- FindNeighbors(seu, reduction = "integrated.harmony")
seu <- FindClusters(seu, resolution = 2, cluster.name = "integrated.clusters")

# check integration
seu <- RunUMAP(
    seu,
    dims = 1:30,
    reduction = "integrated.harmony",
    return.model = T,
    verbose = F
    )

sketch_umap_sample <- DimPlot(
    seu,
    group.by = "gsm_id",
    reduction = "umap",
    label = FALSE
    ) +
    NoLegend()

sketch_umap_major_ctypes <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "Classified_Celltype",
    label = TRUE
    ) +
    NoLegend()

sketch_umap_malignants <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "Malignant"
    )

sketch_umap_treatment <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "Subgroup"
    )



test <- CombinePlots(plots = list(sketch_umap_sample,
                                  sketch_umap_major_ctypes,
                                  sketch_umap_malignants,
                                  sketch_umap_treatment)
                     )

ggsave(
    filename = "results/aml/combined_aml_sketch_umaps_integrated.png",
    plot = test,
    dpi = 300,
    height = 14, 
    width = 20
    )


## Integrate back the full dataset
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
seu <- ProjectIntegration(
    object = seu,
    sketched.assay = "sketch_500",
    assay = "RNA",
    reduction = "integrated.harmony"
    )

seu <- RunUMAP(
    object = seu, 
    reduction = "integrated.harmony.full",
    dims = 1:30,
    reduction.name = "umap.full",
    reduction.key = "UMAP_full_"
    )

saveRDS(object = seu, file = "results/aml/seu.Rds")

DefaultAssay(seu) <- "RNA"

## Merge all layers after integration and clustering
seu <- JoinLayers(object = seu, assay = "RNA")

## Cell cycle
seu <- CellCycleScoring(
    object = seu,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes
)

## Filter and keep malignant cells only
malignants <- subset(seu, subset = Malignant == "Malignant")
saveRDS(object = malignants, file = "results/aml/malignants_seu.rds")

## generate mat and meta for bc
rm(seu)
gc()
bc_metadata <- malignants@meta.data
seu_mat <- as(malignants[["RNA"]]$data, Class = "dgCMatrix")
saveRDS(object = seu_mat, file = "results/aml/normalized_malignants_mat.rds")
write.table(bc_metadata, file = "results/aml/malignants_annotation.tsv")
