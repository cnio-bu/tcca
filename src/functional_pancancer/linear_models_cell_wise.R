library(limma)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = 'v5')


## Load metaprogram data
mp_mat <- readRDS("results/functional/metaprograms_sc_mat.rds")
mp_meta <- data.table::fread("results/functional/cell_level_mps_metadata.tsv")

## subset to malignant cells
idx_right <- setdiff(rownames(mp_mat), mp_meta$original_rownames)

mp_meta_malignants <- mp_meta %>%
    filter(malignancy == TRUE)

## Focus on primary untreated tumors
mp_meta_malignants <- mp_meta_malignants %>%
    filter(treated == "f" & sample_type == "p")

mp_mat <- mp_mat[mp_meta_malignants$original_rownames, ]

# Load metacom data
cell_level_metacoms <- data.table::fread(
    "results/modules/annotated/malignant_cells_best_metacoms_all_cohort.tsv"
)

cell_level_metacoms <- cell_level_metacoms %>%
    filter(
        treated == FALSE,
        sample_type == "p"
        )

## match cell names
print(paste0("Missmatched cells:", nrow(mp_mat) - nrow(cell_level_metacoms)))

diff_names <- setdiff(mp_meta_malignants$new_cell_id, cell_level_metacoms$new_cell_id)

## pretify names
cell_level_metacoms$metacommunity <- as.factor(cell_level_metacoms$metacommunity)
levels(cell_level_metacoms$metacommunity) <- c(
    "Metacommunity 1",
    "Metacommunity 2",
    "Metacommunity 3",
    "Metacommunity 4",
    "Metacommunity 5",
    "Metacommunity 6"
    
)

## Initial lm fit
design <- model.matrix(~0 + tumor_type + metacommunity, data = cell_level_metacoms)
colnames(design)[30:34] <- paste0("metcom_", c(2:6))

fit <- limma::lmFit(object = t(mp_mat), design = design)

## get all contrasts
all_contrasts <- makeContrasts(
 #   metcom_2 - metcom_1,
 #   metcom_3 - metcom_1,
 #   metcom_4 - metcom_1,
 #   metcom_5 - metcom_1,
 #   metcom_6 - metcom_1,
    metcom_3 - metcom_2,
    metcom_4 - metcom_2,
    metcom_5 - metcom_2,
    metcom_6 - metcom_2,
    metcom_4 - metcom_3,
    metcom_5 - metcom_3,
    metcom_6 - metcom_3,
    metcom_5 - metcom_4,
    metcom_6 - metcom_4,
    metcom_6 - metcom_5,
    levels = design
)

fit <- contrasts.fit(fit, contrasts = all_contrasts)
fit2 <- eBayes(fit = fit, robust = TRUE)

res <- decideTests(fit2)

res2 <- topTable(fit2, coef = 10, number = Inf)

res2 <- res2 %>%
    filter("adj.P.Val" <= 0.05)


## try wilcox test
rownames(mp_meta_malignants) <- mp_meta_malignants$new_cell_id
rownames(mp_mat) <- mp_meta_malignants$new_cell_id

cell_level_metacoms <- as.data.frame(cell_level_metacoms)
rownames(cell_level_metacoms) <- cell_level_metacoms$new_cell_id

seu <- CreateSeuratObject(
    counts = t(mp_mat),
    meta.data = cell_level_metacoms
)

Idents(seu) <- seu$metacommunity

res <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    slot = "counts",
    random.seed = 1,
    only.pos = TRUE,
    logfc.threshold = 0.1,
    return.thresh = 0.05
)

table(res$gene)

## keep top 3
res_top <- res %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 3)

## cute res
real_names <- c(
    "metaprograms10" = "Protein reg. Protein maduration",
    "metaprograms11" = "Protein reg. Translation init",
    "metaprograms14" = "Mesenchymal EMT III",
    "metaprograms15" = "Mesenchymal EMT IV",
    "metaprograms17" = "Interferon and MHC II Interferon and MHC II (I)",
    "metaprograms18" = "Interferon and MHC II Interferon and MHC II (II)",
    "metaprograms19" = "Senescence EpiSen",
    "metaprograms23" = "Secreted, Secreted II",
    "metaprograms24" = "Cilia Cilia",
    "metaprograms25" = "Lineage related Neural astrocytes",
    "metaprograms28" = "Lineage related Neural oligo normal",
    "metaprograms29" = "Lineage related Neural NPC and OPC",
    "metaprograms34" = "Lineage related Haemat. Platelet activation",
    "metaprograms37" = "Lineage related Haemat. Haemato. related II",
    "metaprograms4" = "Cell cycle Chromatin",
    "metaprograms40" = "Unnasigned PDAC related",
    "metaprograms8" = "Protein reg. Proteasomal deg."
    
)

res_top$gene <- real_names[res_top$gene]

## Perform DGE
seu <- readRDS("results/tcca/tcca_seurat_raw.rds")
mat <- BPCells::open_matrix_dir(dir = "results/tcca/raw_matrix_tcca")

seu@assays$RNA$counts <- mat

## keep malignant cells only
seu <- subset(seu, subset = malignancy == TRUE)

## subset sc to remove bhupinder pal samples that are predicted tumours
seu <- subset(
    seu,
    subset = study != "brca_bhupinder_pal" | tumor_subtype != "predicted_tumour"
)

## remove sample w/o metacoms
seu <- subset(seu, subset = sample != "T10")

## Generate cell indexes
seu$new_cell_id <- paste0("c", 1:ncol(seu))
    
## subset primary untreated
seu <- subset(seu, subset =  treated == "f" & sample_type == "p")

## test match
c1 <- seu@meta.data[seu@meta.data$new_cell_id == "c1", ]
c1_meta <- cell_level_metacoms[cell_level_metacoms$new_cell_id == "c1", ]

c350 <- seu@meta.data[seu@meta.data$new_cell_id == "c25000", ]
c350_meta <- cell_level_metacoms[cell_level_metacoms$new_cell_id == "c25000", ]

## match n.cells
print(paste0("Missmatched cells:", ncol(seu) - nrow(cell_level_metacoms)))

cells.diff <- setdiff(seu$new_cell_id, cell_level_metacoms$new_cell_id)
print(paste0("Differring cell sets:", length(cells.diff)))

## Add meta.data
colnames(cell_level_metacoms)

colnames(seu) <- seu$new_cell_id

seu <- AddMetaData(
    object = seu,
    metadata = cell_level_metacoms[, "metacommunity"],
    col.name = "metacommunity"
    )

Idents(seu) <- seu$metacommunity
seu <- NormalizeData(seu)

## Find DGE genes
metacom_dge_genes <- FindMarkers(
    object = seu,
    assay = "RNA",
    slot = "data",
    random.seed = 1,
    ident.1 = "Metacommunity 6",
    ident.2 =  "Metacommunity 1"
)


## Cell level analyses crashes, resort to pseudobulk
seu$sample_study <- paste0(seu$sample, "_", seu$study)

bulk <- AggregateExpression(
    object = seu,
    return.seurat = T,
    slot = "counts",
    assays = "RNA",
    group.by = c("metacommunity", "sample_study")
    )

bulk <- NormalizeData(bulk)

res.markers <- FindAllMarkers(
    object = bulk,
    test.use = "wilcox",
    slot = "counts"
)

saveRDS(object = bulk, file = "results/functional/pseudo_bulk_metacom_seurat.rds")
