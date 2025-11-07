library(DESeq2)
library(Seurat)
library(tidyverse)

bulk <- readRDS("results/functional/pseudo_bulk_metacom_seurat.rds")

## SoupX counts must be rounded down
mat <- bulk@assays$RNA$counts
mat <- as.matrix(mat)
mat <- round(mat, digits = 0)

bulk@assays$RNA$counts <- mat

bulk <- NormalizeData(bulk)

markers.metacoms <- FindAllMarkers(
    object = bulk,
    assay = "RNA",
    slot = "counts",
    test.use = "DESeq2",
    logfc.threshold = 0,
    return.thresh = 0.05
)

markers.metacoms2 <- FindAllMarkers(
    object = bulk,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0,
    only.pos = TRUE,
    return.thresh = 0.05
    )

markers.metacoms <- markers.metacoms %>%
    filter(p_val_adj <= 0.05)

write.csv(x = markers.metacoms, file = "results/metacom_markers.csv")
write.csv(x = markers.metacoms2, file = "results/metacom_markers_wilcox.csv")

## test a classic model
meta.data <- bulk@meta.data
colnames(meta.data) <- "metacom"

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = meta.data,
    design = ~metacom
    )

vsd <- vst(dds)

pseudobulk_metacom_pca <- plotPCA(vsd, intgroup = "metacom")
dds2 <- DESeq(object = dds, test = "Wald", fitType = "parametric")

res <- results(
    object = dds2,
    contrast = c("metacom", "Metacommunity 6", "Metacommunity 5"), 
    tidy = TRUE
    )
