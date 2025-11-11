library(Seurat)
library(BPCells)
library(dplyr)
setwd("/storage/scratch01/shared/projects/bc-meta")

# Load expression matrix
seu_lvl2 <- readRDS("./single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
seu_lvl2 <- JoinLayers(seu_lvl2)

# Load beyondcell matrix
bc <- open_matrix_dir("./beyondcell_immuno/full_mat_beyondcell")
metadata <- read.table(
    "./single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t", header = TRUE
)

# Subset breast cancer cells
metadata_brca <- metadata %>%
    filter(malignancy == "True") %>%
    mutate(cell = row_number()) %>%
    filter(refined_tumor_type == "BRCA" & patient != "ccl")

colnames(bc) <- paste0(c(1:ncol(bc)))
norm_mat_brca <- bc[, metadata_brca$cell]
norm_mat_brca <- as(norm_mat_brca, "sparseMatrix")
expr_mat_brca <- seu_lvl2[["RNA"]]$counts[, metadata_brca$old_barcode]
expr_mat_brca <- as(expr_mat_brca, "sparseMatrix")

colnames(expr_mat_brca) <- metadata_brca$cell
    
saveRDS(norm_mat_brca, "./beyondcell_immuno/brca_usecase/bc_sparsemat.rds")
saveRDS(expr_mat_brca, "./beyondcell_immuno/brca_usecase/expr_sparsemat.rds")
write.table(
    metadata_brca,
    file = "./beyondcell_immuno/brca_usecase/metadata_brca.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Create beyondcell object
library(beyondcell)
library(dplyr)

expr_mat_brca <- readRDS("./beyondcell_immuno/brca_usecase/expr_sparsemat.rds")
norm_mat_brca <- readRDS("./beyondcell_immuno/brca_usecase/bc_sparsemat.rds")
metadata_brca <- read.table(
    "./beyondcell_immuno/brca_usecase/metadata_brca.tsv",
    sep = "\t", header = TRUE
)
norm_mat_brca <- as.matrix(norm_mat_brca)
expr_mat_brca <- as.matrix(expr_mat_brca)
rownames(metadata_brca) <- metadata_brca$cell
bc_brca <- new("beyondcell",
    normalized = norm_mat_brca,
    data = norm_mat_brca,
    scaled = matrix(ncol = 0, nrow = 0),
    switch.point = numeric(0),
    ranks = list(),
    expr.matrix = expr_mat_brca,
    meta.data = metadata_brca,
    SeuratInfo = list(),
    background = matrix(ncol = 0, nrow = 0),
    reductions = list(),
    regression = list(),
    n.genes = nrow(norm_mat_brca),
    mode = c("up", "down"),
    thres = 0.1
)
scaled.matrix <- t(apply(bc_brca@normalized, 1, scales::rescale, to = c(0, 1)))
slot(bc_brca, "scaled") <- round(scaled.matrix, digits = 2)

# Then run CreatebcObject to compute SP:
bc_brca <- CreatebcObject(bc_brca)

# Subset primary BRCA samples
cells_primary <- bc_brca@meta.data %>%
    filter(sample_type == "p") %>%
    pull(cell)
bc_brca_primary <- bcSubset(bc_brca, cells = cells_primary)
bc_brca_primary <- bcRecompute(bc_brca_primary)

saveRDS(bc_brca_primary, "./beyondcell_immuno/brca_usecase/bc_brca_primary.rds")

# Clustering
png("./beyondcell_immuno/brca_usecase/pcs.png")
bc_brca_primary <- bcUMAP(bc_brca_primary, k.neighbors = 4, res = 0.2)
dev.off()

bc_brca_primary <- bcUMAP(bc_brca_primary, pc = 10, k.neighbors = 4, res = 0.2)

png("./beyondcell_immuno/brca_usecase/umap.png")
bcClusters(bc_brca_primary, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)
dev.off()

# Regress out nFeature_RNA
bc_brca_primary <- bcRegressOut(bc_brca_primary, vars.to.regress = "nFeature_RNA")

png("./beyondcell_immuno/brca_usecase/pcs2.png")
bc_brca_primary <- bcUMAP(bc_brca_primary, pc = 10, k.neighbors = 20, res = 0.1)
dev.off()

png("./beyondcell_immuno/brca_usecase/umap_after_regression.png")
bcClusters(bc_brca_primary, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)
dev.off()


png("./beyondcell_immuno/brca_usecase/umap_colored_tc.png")
bcClusters(bc_brca_primary, UMAP = "beyondcell", idents = "scTherapy_cluster", pt.size = 1.5, cells = bc_brca_primary@meta.data[!is.na(bc_brca_primary@meta.data$scTherapy_cluster), "cell"])
dev.off()

# Plot drug signatures
library(cowplot)
library(ggplot2)

drugs_oi <- c(
            "sig-20883", "sig-21356", "sig-20886", "sig-21369", "sig-21370",
            "sig-20922", "sig-20900", "sig-21212", "sig-20889", "sig-20921",
            "sig-20995", "sig-20879", "sig-21166", "sig-20911", "sig-21041",
            "sig-20898", "sig-20925", "sig-20926", "sig-20952", "sig-21405",
            "sig-20907", "sig-20932", "sig-21323", "sig-21324", "sig-20890",
            "sig-20941", "sig-20967", "sig-21141", "sig-20980", "sig-21250",
            "sig-21238", "sig-20958", "sig-21297", "sig-21120", "sig-20881",
            "sig-21144", "sig-20949", "sig-20924", "sig-21007", "sig-20904",
            "sig-21194", "sig-20896", "sig-21193", "sig-21039", "sig-20957",
            "sig-21038", "sig-20888", "sig-21182"
        )
brca_drugs <- bcSignatures(
    bc_brca_primary,
    UMAP = "beyondcell",
    signatures = list(
        values = drugs_oi
    ),
    pt.size = 1.5
)


combined <- plot_grid(plotlist = brca_drugs[25:48], ncol = 6, align = "hv")

ggsave(
    "./beyondcell_immuno/brca_usecase/drug_signatures_umap_part2.png",
    plot = combined,
    width = 26, height = 15,
    dpi = 300,
    bg = "white"
)

# Subset cells with scTherapy cluster annotation
bc <- CreatebcObject(bc_brca_primary)
cells_sctherapy <- bc_brca_primary@meta.data %>%
    filter(!is.na(scTherapy_cluster)) %>%
    pull(cell)

bc_brca_sctherapy <- bcSubset(bc, cells = cells_sctherapy)
bc_brca_sctherapy <- bcRecompute(bc_brca_sctherapy)

# Clustering
png("./beyondcell_immuno/brca_usecase/pcs_sctherapy.png")
bc_brca_sctherapy <- bcUMAP(bc_brca_sctherapy, k.neighbors = 4, res = 0.2)
dev.off()

bc_brca_sctherapy <- bcUMAP(bc_brca_sctherapy, pc = 10, k.neighbors = 4, res = 0.2)

png("./beyondcell_immuno/brca_usecase/umap2_sctherapy.png")
bcClusters(bc_brca_sctherapy, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)
dev.off()

# Regress out nFeature_RNA
bc_brca_sctherapy <- bcRegressOut(bc_brca_sctherapy, vars.to.regress = "nFeature_RNA")

png("./beyondcell_immuno/brca_usecase/pcs2_sctherapy.png")
bc_brca_sctherapy <- bcUMAP(bc_brca_sctherapy, pc = 10, k.neighbors = 20, res = 0.1)
dev.off()

png("./beyondcell_immuno/brca_usecase/umap_after_regression_sctherapy.png")
bcClusters(bc_brca_sctherapy, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)
dev.off()

bc_brca_sctherapy@meta.data$scTherapy_cluster <- factor(bc_brca_sctherapy@meta.data$scTherapy_cluster, levels = as.character(1:10))
png("./beyondcell_immuno/brca_usecase/umap_colored_tc_sctherapy.png", res = 300, unit = "in", width = 6, height = 5)
bcClusters(bc_brca_sctherapy, UMAP = "beyondcell", idents = "scTherapy_cluster", pt.size = 0.5, cols = sctherapy_colors)
dev.off()

saveRDS(bc_brca_sctherapy, "./beyondcell_immuno/brca_usecase/bc_brca_sctherapy.rds")


# Drug sensitivity plots for cells with scTherapy cluster annotation
brca_drugs <- bcSignatures(
    bc_brca_sctherapy,
    UMAP = "beyondcell",
    signatures = list(
        values = drugs_oi
    ),
    pt.size = 1.5
)


combined <- plot_grid(plotlist = brca_drugs[25:48], ncol = 6, align = "hv")

ggsave(
    "./beyondcell_immuno/brca_usecase/drug_signatures_umap_sctherapy_part2.png",
    plot = combined,
    width = 26, height = 15,
    dpi = 300,
    bg = "white"
)


# Boxplot of BCS by scTherapy cluster
cluster <- bc_brca_sctherapy@meta.data %>%
    dplyr::select(scTherapy_cluster) %>%
    mutate(scTherapy_cluster = factor(scTherapy_cluster, levels = as.character(1:10)))

bcs_mat <- t(bc_brca_sctherapy@normalized[drugs_oi, ])
df <- cbind(as.data.frame(bcs_mat), cluster)

df_long <- df %>%
    tidyr::pivot_longer(
        cols = starts_with("sig-"),
        names_to = "DrugSignature",
        values_to = "Score"
    )

source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
plot <- ggplot(df_long, aes(x = scTherapy_cluster, y = Score, fill = scTherapy_cluster)) +
    geom_boxplot(outlier.size = 0.5, width = 0.7) +
    facet_wrap(~DrugSignature, scales = "free_y", ncol = 6) +
    theme_bw(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) +
    scale_fill_manual(values = sctherapy_colors) +
    labs(x = "Cluster", y = "Score", title = "Distribución de sensibilidad a drogas por cluster")

ggsave(
    "./beyondcell_immuno/brca_usecase/boxplot_bcs_by_cluster.png",
    plot = plot,
    width = 26, height = 15,
    dpi = 300,
    bg = "white"
)