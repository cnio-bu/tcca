library("beyondcell")
library("biclust")
library("pheatmap")
library("Seurat")
library("data.table")

all_data <- readRDS("../../results/biclustering/breast_sunny_wu.rds")

tumor <- all_data[["CID4515"]]

## get cell cycle scoring
tumor <- CellCycleScoring(
  object = tumor,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

gs <- GetCollection(x = SSc, n.genes = 250, include.pathways = FALSE)
bc <- bcScore(sc = tumor, gs = gs, expr.thres = 0.1)


bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")

bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA", "G2M.Score"))

bc <- bcUMAP(bc = bc, pc = 50, k.neighbors = 20)
bcClusters(bc = bc, idents = "bc_clusters_res.0.2")

full_norm_mat_rgrss <- bc@normalized

write.table(
  x = full_norm_mat_rgrss,
  file = "../../results/biclustering/full_mat_rgrss.tsv",
  sep = "\t"
  )


## same analysis only malignants
filter_malignant <- function(sc) {
  
  types_to_keep <- "Cancer Epithelial"
  if (sum(sc@meta.data$celltype_major == types_to_keep) > 1) {
    sc_filtered <- subset(x = sc, subset = celltype_major == types_to_keep)
    return(sc_filtered)
    
  } else {
    return(NULL)
  }
  
}

tumor_subset <- filter_malignant(tumor)

bc2 <- bcScore(sc = tumor_subset, gs = gs, expr.thres = 0.1)
bc2@normalized[is.na(bc2@normalized)] <- 0
bc2 <- bcRecompute(bc2, slot = "normalized")
bc2 <- bcRegressOut(bc = bc2, vars.to.regress = c("nFeature_RNA"))

#### Full characterization ####
#
# The sample is a breast cancer patient with ~1500 malignant cells. Let's
# thoroughly characterize it.
#
bc2 <- bcUMAP(bc = bc2, pc = 50, k.neighbors = 10, npcs = 50, res = 0.1)
bc_clusts <- bcClusters(bc = bc2, idents = "bc_clusters_res.0.1", spatial = FALSE)

## save the sample bc clusters at 0.1 resolution
ggplot2::ggsave(
  plot = bc_clusts,
  dpi = 100,
  filename = "results/biclustering/bc_clusters_0.1.png"
  )

## Now draw a heatmap of the norm. mat.
norm_mat <- bc2@normalized

heat <- pheatmap(
  mat = norm_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",
  main = "Heatmap of the norm. mat. w/out clustering",
  filename = "results/biclustering/heatmap_no_clustering_norm_mat.png"
)

heat2 <- pheatmap(
  mat = norm_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  cutree_rows = 4,
  cutree_cols = 2,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = bc2@meta.data[, c("celltype_minor", "bc_clusters_res.0.1")],
  scale = "row",
  main = "Heatmap of the norm. mat. with clustering",
  filename = "results/biclustering/heatmap_hc_clustering_norm_mat.png"
)

## Test with norm/rgrss mat
norm_mat <- bc2@normalized

results_biclust <- spectral(
  x = norm_mat,
  normalization = "irrc",
  minr = 5,
  minc = 50,
  withinVar = 3000,
  n_clusters = 3,
  n_best = 3
  )

norm_mat <- norm_mat[, rownames(bc2@meta.data[order(bc2@meta.data$bc_clusters_res.0.1), ])]

moas <- drugInfo$MoAs %>%
  select(IDs, main.MoAs) %>%
  mutate(
    moa_to_use = case_when(
      is.na(main.MoAs) ~ "Other",
      main.MoAs == "" ~ "Other",
      TRUE ~ main.MoAs
    )
  ) %>%
  select(-main.MoAs) %>%
  mutate(
    moa_to_use = as_factor(moa_to_use)
  ) %>%
  distinct(.keep_all = TRUE) %>%
  distinct(IDs, .keep_all = TRUE) %>%
  filter(IDs %in% rownames(norm_mat),) %>%
  mutate(
    is_egfr_inhibitor = case_when(
      moa_to_use == "ErbB inhibitor" ~ "EGFR inhibitor",
      moa_to_use == "PARP inhibitor" ~ "PPAR inhibitor",
      TRUE ~ "Other",
    ),
    is_egfr_inhibitor = as_factor(is_egfr_inhibitor),
  ) %>%
  as.data.frame()

rownames(moas) <- moas$IDs  

renacer_moas <- readr::read_csv(
  "src/spectral_clustering_sample_test/moa_annotation_manual.csv",
  ) %>%
  rename(
    "drug_signature" = ...1
  ) %>%
  mutate(
    drug_signature = toupper(drug_signature)
  )


drug_info_sig <- drugInfo$IDs %>%
  filter(collections == "SSc") %>%
  mutate(
    sig_id = paste(preferred.drug.names, studies, original.IDs, sep = "_")
  )

common_drugs <- intersect(renacer_moas$drug_signature, drug_info_sig$sig_id)

drugs_id_to_keep <- drug_info_sig %>%
  select(IDs, sig_id) %>%
  left_join(renacer_moas, by = c("sig_id" = "drug_signature")) %>%
  select(IDs, new_moa) %>%
  mutate(
    new_moa = case_when(
      is.na(new_moa) ~ "Other",
      TRUE ~ new_moa
    ),
    new_moa = as_factor(new_moa)
  ) %>%
  as.data.frame()

rownames(drugs_id_to_keep) <- drugs_id_to_keep$IDs

heat3 <- pheatmap(
  mat = norm_mat[
    order(results_biclust@info$row_labels), ],
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-5, 5, 0.2),
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = bc2@meta.data[, c("celltype_minor", "bc_clusters_res.0.1")],
#  annotation_row = drugs_id_to_keep,
  scale = "row",
  main = "Heatmap of biclustering",
  filename = "results/biclustering/heatmap_biclust_norm_test.png",
  height = 10,
  width = 10
)


## decompose the matrix to aggregate it to moa level
norm_mat_dt <- data.table(t(norm_mat))
norm_mat_dt[, "cell"] <- colnames(norm_mat)

mat_transposed_melted <- melt(
  norm_mat_dt,
  id.vars = c("cell"),
  variable.name = "signature",
  value.name = "enrichment"
  )

moa_aggregated_mat <- mat_transposed_melted %>%
  right_join(y = drugs_id_to_keep, by = c("signature" = "IDs")) %>%
  group_by(cell, new_moa) %>%
  summarise(enrichment = median(enrichment)) %>%
  ungroup() %>%
  pivot_wider(names_from = cell, values_from = enrichment) %>%
  as.data.frame()  

rownames(moa_aggregated_mat) <- moa_aggregated_mat$new_moa
moa_aggregated_mat$new_moa <- NULL
moa_aggregated_mat <- as.matrix(moa_aggregated_mat)

moa_aggregated_mat <- moa_aggregated_mat[, rownames(bc2@meta.data[order(bc2@meta.data$bc_clusters_res.0.1), ])]

heat4 <- pheatmap(
  mat = moa_aggregated_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-5, 5, 0.2),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = bc2@meta.data[, c("bc_clusters_res.0.1", "celltype_minor")],
  #  annotation_row = drugs_id_to_keep,
  scale = "row",
  main = "Heatmap of biclustering MoA",
  filename = "results/biclustering/heatmap_biclust_moa_norm_test.png",
  height = 10,
  width = 10
)


## second sample: CID4471
second_sample <- all_data[["CID4530N"]]

second_tumor <- filter_malignant(second_sample)

bc3 <- bcScore(sc = second_tumor, gs = gs, expr.thres = 0.1)
bc3@normalized[is.na(bc3@normalized)] <- 0
bc3 <- bcRecompute(bc3, slot = "normalized")
bc3 <- bcRegressOut(bc = bc3, vars.to.regress = c("nFeature_RNA"))

bc3 <- bcUMAP(bc = bc3, pc = 50, npcs = 50, k.neighbors = 10, res = 0.1)

bc_clusters_3 <- bcClusters(bc = bc3, idents = "bc_clusters_res.0.1")

## Now draw a heatmap of the norm. mat.
norm_mat_3 <- bc3@normalized

heat5 <- pheatmap(
  mat = norm_mat_3,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",
  main = "Heatmap of the norm. mat. w/out clustering",
  filename = "results/biclustering/heatmap_no_clustering_norm_mat_3.png"
)

norm_mat_3 <- norm_mat_3[, rownames(bc3@meta.data[order(bc3@meta.data$bc_clusters_res.0.1), ])]

heat6 <- pheatmap(
  mat = norm_mat_3,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-5, 5, 0.2),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  cutree_rows = 4,
  cutree_cols = 2,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = bc3@meta.data[, c("celltype_minor", "bc_clusters_res.0.1")],
  scale = "row",
  main = "Heatmap of the norm. mat. with clustering",
  filename = "results/biclustering/heatmap_hc_clustering_norm_mat_3.png"
)

results_biclust2 <- spectral(
  x = norm_mat_3,
  normalization = "irrc",
  minr = 5,
  minc = 50,
  withinVar = 3000,
  n_clusters = 3,
  n_best = 3
)

heat7 <- pheatmap(
  mat = norm_mat_3[
    order(results_biclust2@info$row_labels), ],
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = seq(-5, 5, 0.2),
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = bc3@meta.data[, c("celltype_minor", "bc_clusters_res.0.1")],
  #  annotation_row = drugs_id_to_keep,
  scale = "row",
  main = "Heatmap of biclustering",
  filename = "results/biclustering/heatmap_biclust_norm_test_2.png",
  height = 10,
  width = 10
)

test <- as.data.frame(results_biclust@info$row_labels)
test2 <- as.data.frame(results_biclust2@info$row_labels)
colnames(test2) <- c("assigned_clusters")

test$results2 <- test2[rownames(test), "assigned_clusters"]
colnames(test) <- c("tnbc", "luma")
test$drug <- rownames(test)

## mod 1 vs mod 1
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

mod1 <- test[test$tnbc == 1, "drug"]
mod2 <- test[test$luma == 1, "drug"]

jaccard(mod1, mod2)

module_test <- intersect(mod1, mod2)
mod1
mod2

res_test <- drugs_id_to_keep[drugs_id_to_keep$IDs %in% module_test, ]
