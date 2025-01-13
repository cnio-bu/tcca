library(BPCells)
library(Seurat)
library(tidyverse)

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 3e+09)


annotate_clinical_data <- function(sc){
    this_meta <- metadata %>%
        filter(patient == unique(sc$orig.ident)) %>%
        as.data.frame()
    
    common_cells <- intersect(colnames(sc), this_meta$Cell_barcode)
    this_meta <- this_meta[this_meta$Cell_barcode %in% common_cells, ]
    rownames(this_meta) <- this_meta$Cell_barcode
    sc <- sc[, common_cells]
    sc <- AddMetaData(sc, metadata = this_meta)
    return(sc)
    
}


all_samples <- list.files("raw/multiple_myeloma_stephan_tirier", full.names = TRUE)
seu_list <- list()

for (sample in all_samples) {
    sample_name <- basename(sample)
    sample_name <- stringr::str_remove(string = sample_name, pattern = "^[^_]*_")
    sample_name <- stringr::str_remove(string = sample_name, pattern = ".csv")
    
    mat <- data.table::fread(sample) %>%
        as.data.frame()
    
    rownames(mat) <- mat$gene
    mat$gene <- NULL 
    mat = as.sparse(mat)
    mat_dir = paste0("results/mmieloma/", sample_name)
    write_matrix_dir(mat = mat, dir = mat_dir, overwrite = TRUE)
    mat = open_matrix_dir(dir = mat_dir)
    
    seu <- Seurat::CreateSeuratObject(
        counts = mat,
        assay = "RNA",
        project = sample_name
    )
    seu_list <- c(seu_list, seu)
}

names(seu_list) <- basename(all_samples)


metadata <- data.table::fread(input = "raw/GSE161801_K43R_metadata_table.csv") %>%
    as.data.frame()

seu_list = lapply(X = seu_list, FUN = annotate_clinical_data)

seu = merge(x = seu_list[[1]], y = seu_list[2:87])

seu = NormalizeData(seu)
seu = FindVariableFeatures(object = seu)

seu = ScaleData(seu)
seu = RunPCA(seu)

# integrate the datasets
seu <- IntegrateLayers(seu, method = HarmonyIntegration)

# cluster the integrated data
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 2, cluster.name = "harmony_clusters")
seu <- RunUMAP(
    seu,
    reduction = "harmony",
    dims = 1:30,
    return.model = T,
    verbose = T
    )

major_cell_types = DimPlot(object = seu, group.by = "major_celltype")

ggsave(
    plot = major_cell_types,
    filename = "results/mmieloma/integrated_samples_major_types.png"
    )

patients_integrated = DimPlot(
    object = seu,
    group.by = "major_celltype",
    split.by = "timepoint",
    alpha = 0.7,
    )

ggsave(
    plot = patients_integrated,
    filename = "results/mmieloma/integrated_samples_points.png",
    width = 21
    )

seu[["RNA"]] = JoinLayers(seu[["RNA"]])

seu$is_malignant = (seu$major_celltype == "PCs" & seu$celltype_1 == "Myeloma")

patients_mix = DimPlot(
    object = seu,
    group.by = "is_malignant",
    alpha = 0.3
    )

patients_mix = patients_mix +
    scale_color_discrete(name = "Malignant annotation", labels = c("TME", "Myeloma")) +
    ggtitle("")

ggsave(plot = patients_mix, filename = "results/mmieloma/integrated_malignants.png")

## Cell cycle
seu <- CellCycleScoring(
    object = seu,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes
    )

## Filter and keep malignant cells only
malignants = subset(seu, subset = is_malignant == TRUE)
saveRDS(object = malignants, file = "results/mmieloma/malignants_seu.rds")

## generate mat and meta for bc
rm(seu)
gc()
bc_metadata = malignants@meta.data
seu_mat <- as(malignants[["RNA"]]$data, Class = "dgCMatrix")
saveRDS(object = seu_mat, file = "results/mmieloma/normalized_malignants.rds")
write.table(bc_metadata, file = "results/mmieloma/malignants_annotation.tsv")
