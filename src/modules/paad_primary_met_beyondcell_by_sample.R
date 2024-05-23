library(beyondcell)

mat <- readRDS("results/paad/normalized_seu_malignant_mat.rds")
meta.data <- read.table("results/paad/seu_metadata_malignants.tsv")

p2_primary <- meta.data[(meta.data$patient == "P02" & meta.data$sample_type == "m"), ]
mat <- mat[, rownames(p2_primary)]

#gs <- GetCollection(SSc, include.pathways = FALSE)
gs <- GenerateGenesets(
    x = "reference/drug_signatures_fold.gmt",
    perform.reversal = FALSE
)

bc <- bcScore(sc = mat, gs = gs, expr.thres = 0.1)
bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")

## load nFeatures from seuv5
bc@meta.data <- p2_primary

bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

## get elbow
bc <- bcUMAP(bc = bc, pc = NULL)
bc <- bcUMAP(bc = bc, pc = 30, k.neighbors = 20, npcs = 50, res = 0.2)

## get TC plot
tcs <- bcClusters(bc = bc, idents = "bc_clusters_res.0.2")

bcSignatures(bc, signatures = list(values = "Gemcitabine_GDSC_1190"))


## simple heatmap
library(ComplexHeatmap)

## get top RV
top_rv <- matrixStats::rowVars(x = bc@normalized)
top_rv <- sort(top_rv, decreasing = TRUE)

## get top RC
top_rcs <- matrixStats::colVars(x = bc@normalized)
top_rcs <- sort(top_rcs, decreasing = TRUE)

## get top 10
top_ten_drugs <- head(top_rv, n = 20)
top_100_cells <- head(top_rcs, n = 500)


## get modules
mods <- read.csv(
    file = "results/test.csv",
    sep = "\t"
    )

drugs_to_keep <- mods[order(mods$tm), "drug"]
drugs_to_keep <- intersect(drugs_to_keep, rownames(bc@normalized))


test <- fabias(X = bc@normalized, p = 5)
bis <- extractBic(fact = test) 

heattop <- ComplexHeatmap::Heatmap(
    matrix = t(scale(t(bc@normalized[names(top_ten_drugs), names(top_100_cells)]))),
    col = circlize::colorRamp2(colors = c("blue", "white", "red"), breaks = c(-3, 0, 3)),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    #column_split = bc@meta.data$bc_clusters_res.0.2,
    #cluster_row_slices = TRUE,
    row_split = 3,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8)
)

heat <- ComplexHeatmap::Heatmap(
    matrix = scale(bc@normalized[drugs_to_keep, ], center = TRUE),
    col = circlize::colorRamp2(colors = c("blue", "white", "red"), breaks = c(-3, 0, 3)),
    cluster_rows = FALSE,
    row_split = mods[mods$drug %in% drugs_to_keep, ]$tm,
    #cluster_columns = TRUE,
    #cluster_column_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    #column_split = bc@meta.data$bc_clusters_res.0.2,
    cluster_row_slices = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6)
)

heat2 <- ComplexHeatmap::Heatmap(
    matrix = cor(t(bc@normalized[drugs_to_keep, ])),
    col = circlize::colorRamp2(colors = c("blue", "white", "red"), breaks = c(-1, 0, 1)),
    cluster_rows = FALSE,
    row_split = mods[mods$drug %in% drugs_to_keep, ]$tm,
    column_split = mods[mods$drug %in% drugs_to_keep, ]$tm,
    #cluster_columns = TRUE,
    #cluster_column_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    #column_split = bc@meta.data$bc_clusters_res.0.2,
    cluster_row_slices = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6)
)
