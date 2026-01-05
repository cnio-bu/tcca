library(Seurat)
library(BPCells)
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(cluster)
library(dynamicTreeCut)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy")
source(file = "~/bc-meta/figures/TCCA_palette.R")

data <- read.table("full_table_drug_prediction.tsv")
data$study_sample <- paste0(sub("\\..*", "", data$Subclone), "_", data$Sample)

drug_subclone_mat <- data %>%
  select(Subclone, Drug_Name, Response) %>%
  pivot_wider(names_from = Subclone,
              values_from = Response,
              values_fill = NA) %>%
  arrange(Drug_Name) %>%
  column_to_rownames(var = "Drug_Name")

# Keep drugs with predicted in more than 1% of subclones
keep_drugs <- apply(drug_subclone_mat, 1, function(drug) {
  keep <- sum(!is.na(drug)) > 0.05 * length(drug)
  return(keep)
})

drug_subclone_mat <- drug_subclone_mat[keep_drugs, ]
drug_subclone_mat[is.na(drug_subclone_mat)] <- "Not predicted"
response_levels <- c("High", "High-to-moderate", "Moderate", "Not predicted")
colors <- c("#cd322f", "#ebad55", "#4284b5", "#b9b9b9")
names(colors) <- response_levels

# Convert matrix into factors
drug_subclone_mat <- apply(drug_subclone_mat, 2, function(x)
  factor(x, levels = response_levels))

### Perform hierarchical clustering
drug_subclone_encoded <- drug_subclone_mat
drug_subclone_encoded[drug_subclone_encoded == "High"] <- 3
drug_subclone_encoded[drug_subclone_encoded == "High-to-moderate"] <- 2
drug_subclone_encoded[drug_subclone_encoded == "Moderate"] <- 1
drug_subclone_encoded[drug_subclone_encoded == "Not predicted"] <- 0

# Convert to numeric matrix
drug_subclone_encoded <- apply(drug_subclone_encoded, 2, as.numeric)
colnames(drug_subclone_encoded) <- colnames(drug_subclone_mat)
rownames(drug_subclone_encoded) <- rownames(drug_subclone_mat)

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/sctherapy")

# Distance matrix is computed betweet rows so we transpose the matrix
dist_matrix <- dist(t(drug_subclone_encoded), method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
png(
  file = "hc_subclones.png",
  res = 300,
  width = 20,
  height = 10,
  units = "in"
)
plot(hc, labels = FALSE, main = "Dendrogram")
abline(h = 25, col = "red")
dev.off()

# Compute silhoutte score for different cuts
sil_scores <- c()
for (h in seq(5, 150, by = 5)) {
  clusters <- cutree(hc, h = h)
  sil_scores <- c(sil_scores, mean(silhouette(clusters, dist_matrix)[, 3]))
}
names(sil_scores) <- seq(5, 150, by = 5)

# Gap statistic
# gap_stat <- clusGap(t(drug_subclone_encoded), FUN = hcut, K.max = 75, B = 100)
# fviz_gap_stat(gap_stat)

# Automatically detect clusters based on dendrogram shape and stabilit
clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist_matrix))
table(clusters)
png(
  file = "hc_subclones.png",
  res = 300,
  width = 20,
  height = 10,
  units = "in"
)
plot(hc, labels = FALSE, main = "Dynamic Tree Cut Clusters")
rect.hclust(hc, k = length(unique(clusters)), border = "red")
dev.off()

# Therefore, the best cut is for 16 clusters
clusters <- cutree(hc, k = 16)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy")
## Create top annotations for samples.
clinical <- data.table::fread("../seurat/v5/clinical_metadata_v4_clean.tsv")
clinical$study_sample <- paste0(clinical$study, "_", clinical$sample)

seu <- readRDS("../seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
new_sex <- seu@meta.data %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  select(study_sample, sex) %>%
  distinct()

rownames(data) <- NULL
subclones <- data %>%
  select(Sample, study_sample, Subclone) %>%
  distinct() %>%
  column_to_rownames(var = "Subclone")

# Add clincal annotation to subclones
clinical_subclones <- subclones %>%
  left_join(clinical, by = "study_sample") %>%
  select(-sex) %>%
  left_join(new_sex, by = "study_sample")
rownames(clinical_subclones) <- rownames(subclones)

translat_human_sites <- c(
  "bone_marrow" = "Bone marrow",
  "brain" = "Brain",
  "adrenal_gland" = "Adrenal gland",
  "breast" = "Breast",
  "skin" = "Skin",
  "esophagus" = "Esophagus",
  "oesophagus" = "Esophagus",
  "liver" = "Liver",
  "lung" = "Lung",
  "lymph_node" = "Lymph node",
  "other" = "Other",
  "ovary" = "Ovary",
  "pancreas" = "Pancreas",
  "prostate" = "Prostate",
  "soft_tissue" = "Soft tissue",
  "bladder" = "Bladder",
  "colon" = "Colon",
  "kidney" = "Kidney"
)

clinical_subclones <- clinical_subclones %>%
  mutate(
    summarised_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "Other"
    ),
    adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
    is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
    treated = ifelse(treated == "t", "Treated", "Untreated"),
    sex = ifelse(sex == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
  )

clinical_subclones$hc_clusters <- factor(clusters, levels = unique(clusters))
subclone_annot_df <- clinical_subclones %>%
  select(
    sex,
    adult_pediatric,
    is_blood,
    sample_type,
    summarised_tumor_site,
    treated,
    hc_clusters
  ) %>%
  as.data.frame()

subclone_annot_df$summarised_tumor_site <- translat_human_sites[subclone_annot_df$summarised_tumor_site]


colnames(subclone_annot_df) <- c(
  "Chromosomal sex",
  "Age group",
  "Solid/Liquid",
  "Sample type",
  "Sample site",
  "Treatment",
  "Hierarchical clusters"
)
names(mps_colors) <- unique(clusters)
pals <- list(
  "Chromosomal sex" = sex_colors,
  "Age group" = age_colors,
  "Solid/Liquid" = sl_colors,
  "Sample type" = pm_colors,
  "Sample site" = tumor_sites_colors,
  "Treatment" = treatment_colors,
  "Hierarchical clusters" = mps_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = subclone_annot_df,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_legend = c("Sample site" = FALSE, "Hierarchical clusters" = FALSE)
)


# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  ncol = 3,
  gap = unit(10, "mm"),
  title = "Sample site"
)
hc_clustering_legend <- Legend(
  at = names(pals$`Hierarchical clusters`),
  legend_gp = gpar(fill = pals$`Hierarchical clusters`),
  ncol = 3,
  gap = unit(10, "mm"),
  title = "Hierarchical clusters"
)


png(
  file = "/storage/scratch01/users/mgonzalezb/bc-meta/sctherapy/heatmap_drug_subclones_rm.png",
  res = 500,
  width = 19,
  height = 16,
  units = "in"
)

# Plot the heatmap
ht <- Heatmap(
  drug_subclone_mat,
  col = colors,
  na_col = "#b9b9b9",
  # Color for missing values
  top_annotation = top_annotation,
  # cluster_rows = TRUE,
  # cluster_row_slices = TRUE,
  # cluster_columns = TRUE,
  # cluster_column_slices = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_side = "left",
  column_dend_side = "bottom",
  column_dend_height = unit(3, "cm"),
  cluster_columns = hc,
  cluster_rows = hclust(dist(drug_subclone_encoded, 
                             method = "euclidean"), 
                        method = "ward.D2"),
  row_dend_width = unit(3, "cm"),
  column_split = 16,
  name = "Response Level",
  row_title = "Drugs",
  column_title = "Subclones",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  row_names_side = "right",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  heatmap_width = unit(16, "in"),
  heatmap_height = unit(14, "in")
)
#   cluster_rows = TRUE,
#   cluster_columns = TRUE,
#   cluster_row_slices = TRUE,
#   name = "Response Level",
#   row_title = "Drugs",
#   column_title = "Samples",
#   column_names_side = "top",
#   row_names_side = "left",
#   show_column_dend = TRUE
ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
  "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
  "legend_gap" = unit(1, "cm")
)
draw(
  ht,
  annotation_legend_side = "top",
  annotation_legend_list = list(tumor_site_legend)
)
dev.off()
