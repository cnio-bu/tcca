library(BPCells)
library(ComplexHeatmap)
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

seu <- SketchData(
    object = seu,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
)

DefaultAssay(seu) <- "sketch"


## Extract matrix of leveraged cells for heatmap
sketched_mat <- seu[["sketch"]]$data
sketched_mat <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

## extract meta group annotation
mt1 <- read_tsv(file = "results/modules/annotated/patients_non_treated_meta_groups.tsv") %>%
    select(meta_community, signature, collapsed.MoAs) %>%
    group_by(meta_community, signature) %>%
    mutate(
        n.repeats = n()
    ) %>%
    select(signature, meta_community, collapsed.MoAs) %>%
    as.data.frame()

samples_to_keep <- seu@meta.data %>%
    filter(
        sample_type == "p" & treated == "f"
    ) %>%
    filter(
        study != "cell_lines_gabriella_kinker"
    ) 

cells_to_keep <- rownames(samples_to_keep)

sketched_to_keep <- intersect(cells_to_keep, colnames(sketched_mat))


rownames_annot <- data.frame(
    "signature" = rownames(sketched_mat)
) %>%
    left_join(
        y = mt1
    ) %>%
    arrange(meta_community) %>%
    group_by(signature) %>%
    slice_head() %>%
    arrange(meta_community) %>%
    mutate(
        meta_community = as.character(meta_community),
        meta_community = replace_na(meta_community, "None"),
        meta_community = as_factor(meta_community)
    ) %>%
    as.data.frame()

rownames(rownames_annot) <- rownames_annot$signature
rownames_annot$signature <- NULL

piti = c("#6cca8e","#8398dc","#ea95ae","#1dade6", "#ff5f76", "#ffb6b6","#fff154","#ba7fff","#ca9a8c", "#4b71e5",# "#cccccc",
         "#ff6600","#add82f","#ff3333","#0dba3c", "#ff864c", "#c4ea94","#666699","#888888","#b8c0ba", "#d58aca","#6da753","#ca9a8c","#ff4430","#e06d23")

has_comm <- rownames(rownames_annot[rownames_annot$meta_community != "None", ])
smaller_sketch <- sketched_mat[has_comm, sketched_to_keep]

subset_annot <- rownames_annot[rownames_annot$meta_community != "None", ]

moa_cols <- piti[1:length(unique(subset_annot$collapsed.MoAs))]
names(moa_cols) <- unique(subset_annot$collapsed.MoAs)


meta_comms <- piti[1:length(unique(subset_annot$meta_community))]
names(meta_comms) <- unique(subset_annot$meta_community)

## Prepare heatmap annotations
right_palette <- list(#"Mechanism of action" = moa_cols,
    "Meta-community" = meta_comms
)

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    #   "Mechanism of action" = rownames_annot$collapsed.MoAs,
    "Meta-community" = subset_annot$meta_community,
    which = "row",
    show_annotation_name = FALSE,
    col = right_palette
)

top_annot <- ComplexHeatmap::HeatmapAnnotation(
    "Tumor of origin" = seu@meta.data[colnames(smaller_sketch), "refined_tumor_site"],
    which = "col"
)


png(filename = "results/test.png", width = 24, height = 48, units = "in", res = 100)
ComplexHeatmap::Heatmap(
    matrix = smaller_sketch,
    col = circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"), transparency = 0.3),
    cluster_rows = FALSE,
    cluster_row_slices = TRUE,
    clustering_distance_rows = "pearson",
    row_order = subset_annot[order(subset_annot$meta_community), "signature"],
    row_split = subset_annot$meta_community,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    clustering_distance_columns = "pearson",
    column_split = 10,
    show_row_names = FALSE,
    show_column_names = FALSE,
    right_annotation = right_annotation,
    top_annotation = top_annot 
)
dev.off()


# test
smaller_sketch_dt <- as.data.frame(smaller_sketch) %>%
    as_tibble(rownames = "signature") %>%
    pivot_longer(
        cols = all_of(2:27760),
        names_to = "cell",
        values_to = "enrichment"
    )

cell_level_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, sample, study)


smaller_sketch_dt_annot <- smaller_sketch_dt %>%
    left_join(
        y = cell_level_meta,
        by = "cell"
    ) %>%
    left_join(
        y = mt1,
        by = "signature"
    ) %>%
    filter(
        !is.na(meta_community)
    ) %>%
    group_by(cell, meta_community) %>%
    summarise(
        average_enrichment = median(enrichment)
    )

sketch_remat <- smaller_sketch_dt_annot %>%
    pivot_wider(id_cols = meta_community, names_from = cell, values_from = average_enrichment) %>%
    as.data.frame()

rownames(sketch_remat) <- sketch_remat$meta_community
sketch_remat$meta_community <- NULL
sketch_remat <- as.matrix(sketch_remat)

png(filename = "results/test_average_meta.png", width = 24, height = 48, units = "in", res = 100)
ComplexHeatmap::Heatmap(
    matrix = sketch_remat,
    col = circlize::colorRamp2(c(-2,0,2), c("blue", "white", "red"), transparency = 0.3),
  #  cluster_rows = FALSE,
  #  cluster_row_slices = TRUE,
  #  clustering_distance_rows = "pearson",
  #  row_order = subset_annot[order(subset_annot$meta_community), "signature"],
  #  row_split = subset_annot$meta_community,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
  #  clustering_distance_columns = "pearson",
  #  column_split = 10,
  #  show_row_names = FALSE,
    show_column_names = FALSE,
  #  right_annotation = right_annotation,
  #  top_annotation = top_annot 
)
dev.off()




## CLUSTERING

### SEURAT MODES
options(Seurat.object.assay.version = "v5")

DefaultAssay(seu) <- "sketch"
seu <- FindVariableFeatures(seu, nfeatures = 100, selection.method = "dispersion")

seu <- ScaleData(seu)
seu <- FindNeighbors(seu, features = VariableFeatures(seu))
seu <- FindClusters(
    seu,
    resolution = 2,
    method = "igraph",
    group.singletons = TRUE,
    cluster.name = "sketched.clusters"
    )

seu <- RunUMAP(
    seu,
    dims = NULL,
    features = VariableFeatures(seu),
    return.model = T,
    assay = "sketch"
)
    

nonsense_plot <- DimPlot(
    seu,
    label = T,
    label.size = 3,
    reduction = "umap",
    pt.size = 1
    ) + NoLegend()


ggsave(
    plot = nonsense_plot,
    filename = "results/beyondcell_bp/bc_malignants_umap_default.png",
    height = 7,
    width = 7
    )

seu <- ProjectData(
    object = seu,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:50,
    refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(seu) <- "RNA"

DimPlot(
    seu, label = T,
    label.size = 3,
    reduction = "umap",
    group.by = "cluster_full",
    alpha = 0.1
    ) + NoLegend()





####### TEST WITH SKET
sketched_mat_scaled <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

# calculate the variance for each gene
rv <- matrixStats::rowVars(sketched_mat_scaled)

library(ggplot2)

sketched_rv_dist <- ggplot(data = as.data.frame(rv), aes(x=rv)) +
    geom_density() +
    geom_rug() +
    scale_x_continuous(name = "Drug variance", n.breaks = 10) +
    theme_bw() +
    theme(
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )


ggsave(
    plot = sketched_rv_dist,
    filename = "results/beyondcell_bp/pancancer_drug_variance.png",
    height = 7,
    width = 7,
    dpi = 100
    )

## keep drugs with > 2 var
drugs_to_keep <- rv[rv >= 2]

sketched_mat_best <- sketched_mat_scaled[names(drugs_to_keep), ]

# perform a PCA on the data in assay(x) for the selected genes
sketched_mat_scaled[sketched_mat_scaled < 0] <- 0
pca <- prcomp(sketched_mat_scaled, center = FALSE, scale. = FALSE)

plot(pca$rotation[, 3], pca$rotation[, 4])
