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
    ncells = 30000,
    method = "LeverageScore",
    var.name = "leverage.score",
    sketched.assay = "sketch_10k",
    over.write = TRUE,
    seed = 120394
    )

saveRDS(
    object = seu,
    file = "results/aml/seu.Rds",
    destdir = "results/aml", 
    overwrite = TRUE
)
