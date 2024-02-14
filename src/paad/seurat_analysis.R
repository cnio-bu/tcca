library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)

all_samples <- list.dirs("raw/pdac_shu_zhang",
                         recursive = FALSE,
                         full.names = TRUE
)

for(sample in all_samples){
    sample_name <- basename(sample)
    mat <- Seurat::Read10X(data.dir = sample)
    ## cast to bp cells
    write_matrix_dir(mat = mat, dir = paste0(sample, "_bp"))
}

## generate seurat object
all_samples <- paste0(all_samples, "_bp")

seu_list <- list()
for(sample in all_samples){
    mat <- open_matrix_dir(dir = sample)
    sample_name <- gsub(pattern = "_bp", replacement = "", x = basename(sample))
    
    seu <- CreateSeuratObject(counts = mat, assay = "RNA", project = sample_name)
    seu_list <- c(seu_list, seu)
    
}

seu <- merge(seu_list[[1]], seu_list[2:length(seu_list)])

rm(seu_list)
gc()

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(object = seu)
seu <- ScaleData(seu)

## add non malignant annotation from bc-meta
curated_cell_annot <- data.table::fread(
    "results/annotation/harmonized_metadata_malignant_and_nonmalignant.csv"
    )

curated_pdac <- curated_cell_annot %>%
    filter(study == "pdac_shu_zhang") %>%
    mutate(cell_new_id = gsub(pattern = "_28", replacement = "", x = cell))

curated_pdac <- as.data.frame(curated_pdac)
rownames(curated_pdac) <- curated_pdac$cell_new_id

seu_annot <- seu[, curated_pdac$cell_new_id]
seu_annot <- AddMetaData(seu_annot, metadata = curated_pdac)

## integrate
seu_annot <- FindVariableFeatures(seu_annot)
seu_annot <- ScaleData(seu_annot)
seu_annot <- RunPCA(seu_annot)

seu_annot <- IntegrateLayers(
    object = seu_annot,
    method = HarmonyIntegration,
    new.reduction = "integrated.harmony"
    )

# cluster the integrated data
seu_annot <- FindNeighbors(
    seu_annot,
    reduction = "integrated.harmony"
    )

seu_annot <- FindClusters(
    seu_annot,
    resolution = 2,
    cluster.name = "integrated.clusters"
    )

# check integration
seu_annot <- RunUMAP(
    seu_annot,
    dims = 1:30,
    reduction = "integrated.harmony"
)

expr.clusters <- DimPlot(seu_annot, group.by = "integrated.clusters")
cell_types <- DimPlot(
    seu_annot,
    reduction = "umap",
    group.by = "cell_type_main",
    label = TRUE
) +
    NoLegend()

malignants <- DimPlot(
    seu_annot,
    reduction = "umap",
    group.by = c("malignancy", "tumor_site")
    )

patients <- DimPlot(seu_annot, reduction = "umap", group.by = "patient")
seu_annot <- JoinLayers(object = seu_annot)
seu_malignants <- subset(seu_annot, subset = cell_type_main == "Malignant")

saveRDS(seu_malignants, file = "results/paad/seurat_malignants.rds")

## export mat and meta.data for bc
norm_mat <- seu_malignants@assays$RNA$data
saveRDS(object = as.matrix(norm_mat), file = "results/paad/normalized_seu_malignant_mat.rds")

## export metadata
meta.data <- seu_malignants@meta.data
write.table(meta.data, file = "results/paad/seu_metadata_malignants.tsv")
