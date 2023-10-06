library(BPCells)
library(uwot)
library(tidyverse)
library(Seurat)

all_mats <- list.dirs(
    path = "results/beyondcell_bp",
    full.names = TRUE
    )


# yay!
cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

all_mats <- all_mats[2:37]

mats <- map(
    all_mats,
    open_matrix_dir
)

full_mat <- do.call(cbind.fill, mats)
full_mat[is.na(full_mat)] <- 0
full_mat <- as(full_mat, "sparseMatrix")

write_matrix_dir(
    mat = full_mat,
    dir = "results/beyondcell_bp/full_mat_beyondcell"
    )

rm(full_mat, mats)
gc()

mat <- open_matrix_dir(dir = "results/beyondcell_bp/full_mat_beyondcell/")

## load metadata
all.meta <- list.files(
    "results/beyondcell_bp/",
    pattern = "*.tsv",
    full.names = TRUE
    )


meta.data <- all.meta %>%
    map(read.table, row.names = 1) 

meta.data[[27]]$patient <- meta.data[[27]]$orig.ident

for(i in c(1:length(all.meta))){
    print(i)
    this_study <- all.meta[[i]]
    meta.data[[i]]$study <- this_study
}

meta.data_full <-  meta.data  %>%
    map(~.x %>%
             select(
                 nCount_RNA,
                 nFeature_RNA,
                 sample,
                 patient,
                 percent.mt,
                 percent.ribo,
                 S.Score,
                 G2M.Score,
                 Phase,
                 study
                 ) %>%
             mutate(
                 "sample" = as.character(sample),
                 "patient" = as.character(patient)
             )) %>%
    bind_rows() %>%
    mutate(
        study = basename(study),
        study = gsub(pattern="*.tsv", replacement="", x = study)
    )


## get clinical data
clinical <- data.table::fread("results/annotation/clinical_metadata_v2_clean.tsv")
clinical[
    clinical$sample == "T19" &
        clinical$study == "adrenalnb_rui_chong",
    "sample"
] <- "T19_1"

meta.data_full <- meta.data_full %>%
    mutate(
        sample = case_when(
            study == "adrenal_nb_rui_chong" & sample == "T19" ~ "T19_1",
            TRUE ~ sample
        )
    )

## Add clinical metadata
meta.data_full_clinical <- meta.data_full %>%
    rownames_to_column("cell") %>%
    left_join(
        y = clinical,
        by = c("sample" = "sample", "study" = "study")
    ) %>%
    select(-patient.x) %>%
    rename("patient" = patient.y)  %>%
    mutate(
        refined_tumor_site = case_when(
            refined_tumor_site == "" ~ "Unknown",
            TRUE ~ refined_tumor_site
        )
    )

write_tsv(
    x = meta.data_full_clinical, 
    file = "results/annotation/beyondcell_metadata_with_clinical.tsv"
    )

meta.data_full_clinical <- meta.data_full_clinical %>%
    as.data.frame()

rownames(meta.data_full_clinical) <- meta.data_full_clinical$cell
meta.data_full_clinical$cell <- NULL

options(Seurat.object.assay.version = 'v5')
options(future.globals.maxSize = 1e9)

seu <- Seurat::CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    project = "beyondcell_pancancer",
    meta.data = meta.data_full_clinical
    )

seu[["RNA"]]$data <- mat
seu <- FindVariableFeatures(seu, nfeatures = 250)

seu <- SketchData(
    object = seu,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
)

DefaultAssay(seu) <- "sketch"
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 100)
seu <- FindNeighbors(seu, dims = 1:100)
seu <- FindClusters(seu, resolution = 2)

piti = c("#6cca8e","#8398dc","#ea95ae","#1dade6", "#ff5f76", "#ffb6b6","#fff154","#ba7fff","#ffdd56", "#4b71e5",# "#cccccc",
         "#ff6600","#add82f","#ff3333","#0dba3c", "#ff864c", "#c4ea94","#666699","#888888","#b8c0ba", "#d58aca","#6da753","#ca9a8c","#ff4430","#e06d23")



seu <- RunUMAP(
    seu,
    dims = 1:100,
    metric = "correlation",
    return.model = T,
    repulsion.strength = 2,
    min.dist = 0.0001,
    init = "spectral",
    n.neighbors = 50
    )

DimPlot(
    seu,
 #   cols = piti[1:length(unique(seu@meta.data$treated))],
    label = T,
    label.size = 3,
    reduction = "umap",
    pt.size = 0.5,
    group.by = "refined_tumor_site",
    
    )


## Replot to full
seu <- ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
#    dims = 1:100,
#    refdata = list(cluster_full = "seurat_clusters")
)


# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(seu) <- "RNA"

dim(mat)

test <- seu[["sketch"]]$counts
test <- as.matrix(test)
test <- scale(x = test, center = TRUE, scale = TRUE)
test[test > 10] = 5
test[test < 10] = -5

umap <- uwot::umap(
    X = test, n_neighbors = 50,
    metric = "correlation",
    n_components = 3,
    scale = FALSE,
    init = "spectral",
    min_dist = 0.001
    ) 

plot(umap[, 1], umap[, 2])
plot(umap[, 2], umap[, 3])
