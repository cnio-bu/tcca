library(dplyr)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(dynamicTreeCut)
library(dendextend)
library(factoextra)
library(kernlab)
library(igraph)
set.seed(123)
setwd("/home/lmgonzalezb/Documents/bc-meta/sctherapy/")
source(file = "../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

### FUNCTIONS ###
jaccard_similarity <- function(vector1, vector2) {
  intersection <- length(intersect(vector1, vector2))
  union <- length(union(vector1, vector2))
  return(intersection / union)
}


data <- read.table("full_table_drug_prediction.tsv")
data$study_sample <- paste0(sub("\\..*", "", data$Subclone), "_", data$Sample)

drug_subclone <- data %>%
  select(Subclone, Drug_Name) %>%
  distinct()

subclones <- unique(data$Subclone)

# Create a list with drug predictions per subclone
subclones_drug_list <- lapply(subclones, function(subclone) {
  drugs <- drug_subclone$Drug_Name[drug_subclone$Subclone == subclone]
  return(drugs)
})

names(subclones_drug_list) <- subclones

# Compute similarity matrix
similarity_matrix <- matrix(0, nrow = length(subclones), ncol = length(subclones),
                            dimnames = list(subclones, subclones))

# Compute pairwise Jaccard similarity
for (i in seq_along(subclones)) {
  for (j in seq_along(subclones)) {
    similarity_matrix[i, j] <- jaccard_similarity(subclones_drug_list[[subclones[i]]], 
                                                  subclones_drug_list[[subclones[j]]])
  }
}

saveRDS(similarity_matrix, "jaccard_matrix.rds")

similarity_matrix <- readRDS("jaccard_matrix.rds")
# Hierarchical clustering
# Distance matrix is computed betweet rows so we transpose the matrix
hc <- hclust(as.dist(1 - similarity_matrix))
png(
    file = "figures/hc_jaccard.png",
    res = 300,
    width = 20, 
    height = 10,
    units = "in"
)
plot(hc, labels = FALSE, main = "Dendrogram")
abline(h = 0.8, col = "red")
dev.off()
clusters <- cutreeDynamic(hc, distM = as.matrix(as.dist(1 - similarity_matrix)), method = "tree")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, clusters = clusters)
png(
  file = "hc_jaccard_colored.png",
  res = 300,
  width = 20, 
  height = 10,
  units = "in"
)
labels(dend) <- rep("", length(labels(dend)))
plot(dend, main = "Dendrogram")
dev.off()


# Check spectral clustering
spectral_result <- specc(similarity_matrix, centers = 19)
cluster_assignment <- as.factor(spectral_result)
cluster_assignment <- saveRDS(cluster_assignment, "speclustering.rds")
names_subclones <- names(cluster_assignment)
cluster_assignment <- as.character(cluster_assignment)
cluster_assignment <- recode(cluster_assignment, "1" = "8", "2" = "5", "3" = "1", 
                             "4" = "2", "5" = "9", "6" = "6", "7" = "10", "8" = "4",
                             "9" = "3", "10" = "7")
cluster_assignment <- factor(cluster_assignment, levels = as.character(1:10))
names(cluster_assignment) <- names_subclones
saveRDS(cluster_assignment, "speclustering_reordered.rds")

# Plot heatmap of similarity matrix

## Create top annotations for samples.
clinical <- data.table::fread("../clinical_metadata_v4_clean.tsv")
clinical$study_sample <- paste0(clinical$study, "_", clinical$sample)

# Include the inferred sex
seu <- readRDS("../seu_lvl2_sex_inferred.rds")
new_sex <- seu@meta.data %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(study_sample, sex) %>%
    distinct()

# Include TME subtypes
tcca_annot <- read.table("/home/lmgonzalezb/Documents/bc-meta/cohort_statistics/tcca_annotation_raw.tsv",
                         header = TRUE)
tme <- tcca_annot %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  select(study_sample, tme_archetype) %>%
  distinct()

data$study_sample <- paste0(sub("\\..*", "", data$Subclone), "_", data$Sample)
rownames(data) <- NULL
subclones <- data %>%
    select(study_sample, Subclone) %>%
    distinct() %>%
    column_to_rownames(var = "Subclone")

# Add clinical and tme annotation to subclones
clinical_subclones <- subclones %>%
    left_join(clinical, by = "study_sample") %>% 
    select(-sex) %>%
    left_join(new_sex, by = "study_sample") %>%
    left_join(tme, by = "study_sample")
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
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        cluster = cluster_assignment)


subclone_annot_df <- clinical_subclones %>%
    select(
        sex,
        adult_pediatric,
        is_blood,
        sample_type,
        summarised_tumor_site,
        treated,
        tme_archetype,
        cluster) %>%
    as.data.frame()

subclone_annot_df$summarised_tumor_site <-  translat_human_sites[subclone_annot_df$summarised_tumor_site]

write.table(subclone_annot_df, "annotations_subclones.tsv")

colnames(subclone_annot_df) <- c(
    "Chromosomal sex",
    "Age group",
    "Solid/Liquid",
    "Sample type",
    "Sample site",
    "Treatment",
    "TME archetype",
    "Cluster"
)

pals <- list(
    "Chromosomal sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Sample type" = pm_colors,
    "Sample site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "TME archetype" = tme_colors,
    "Cluster" = mps_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df =  subclone_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = FALSE
    )



# Customize legends 
sex_legend <- Legend(
  at = names(pals$`Chromosomal sex`),
  legend_gp = gpar(fill = pals$`Chromosomal sex`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Chromosomal sex"
)
age_legend <- Legend(
    at = names(pals$`Age group`),
    legend_gp = gpar(fill = pals$`Age group`),
    ncol = 1,  # Split Group legend into 2 columns
    gap = unit(10, "mm"),
    title = "Age group"
  )
solid_liquid_legend <- Legend(
  at = names(pals$`Solid/Liquid`),
  legend_gp = gpar(fill = pals$`Solid/Liquid`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Solid/Liquid"
)
sample_type_legend <- Legend(
  at = names(pals$`Sample type`),
  legend_gp = gpar(fill = pals$`Sample type`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Sample type"
)

treatment_legend <- Legend(
  at = names(pals$`Treatment`),
  legend_gp = gpar(fill = pals$`Treatment`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Treatment"
)

tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  ncol = 2,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Sample site"
)

tme_legend <- Legend(
  at = names(pals$`TME archetype`),
  legend_gp = gpar(fill = pals$`TME archetype`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "TME archetype"
)

cluster_legend <- Legend(
  at = names(pals$`Cluster`),
  legend_gp = gpar(fill = pals$`Cluster`),
  ncol = 1,  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Cluster"
)
similarity_matrix <- round(similarity_matrix, 3)

# Plot the heatmap
jaccard_dist <- as.dist(1-similarity_matrix)
heat <- ComplexHeatmap::Heatmap(
    similarity_matrix,
    #col = colorRamp2(c(0, 1), hcl_palette = "Inferno", reverse = TRUE),
    top_annotation = top_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = rownames(subclone_annot_df[order(subclone_annot_df$Cluster),]),
    column_order = rownames(subclone_annot_df[order(subclone_annot_df$Cluster),]),
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    row_split = subclone_annot_df$Cluster,
    column_split = subclone_annot_df$Cluster,
    clustering_distance_rows = jaccard_dist,
    clustering_distance_columns = jaccard_dist,
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    row_dend_side = "left",
    column_dend_side = "bottom",
    name = "Jaccard Index",
    row_title = "Subclones",
    row_title_side = "left",
    column_title = "Subclones",
    column_title_side = "bottom",
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 8),
    row_names_side = "right",
    heatmap_legend_param = list(
        title = "Similarity\n(Jaccard index)",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        title_gap = unit(10, "mm")),
    heatmap_width = unit(10, "in"),
    heatmap_height = unit(10, "in")
)

png(
  file = "figures/heatmap_jaccard_spectralclust10_reordered.png",
  res = 500,
  width = 16,
  height = 14,
  units = "in"
)

ht_opt(
    "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
    "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
    "legend_gap" = unit(1, "cm")
)

pd <- packLegend(sex_legend, 
                 age_legend, 
                 solid_liquid_legend, 
                 sample_type_legend, 
                 treatment_legend, 
                 max_height = unit(8, "cm"), 
                 column_gap = unit(1, "cm"))
draw(heat, 
     annotation_legend_side = "top", 
     heatmap_legend_side = "right", 
     annotation_legend_list = list(pd, tumor_site_legend, tme_legend, cluster_legend))
dev.off()

# Extract orders and dendograms
row_order <- row_order(heat)
column_order <- column_order(heat)
column_dend <- column_dend(heat)

# Compute mean Jaccard Index per cluster
clusters <- lapply(column_order, function(x) colnames(similarity_matrix)[x])
mean_jaccard_index <- lapply(column_order, function(x){
  subclones <- colnames(similarity_matrix)[x]
  mean_jaccard <- mean(similarity_matrix[subclones, subclones])
  return(mean_jaccard)
})

# Get top 10 drugs per cluster of subclones
drugs_per_cluster <- lapply(clusters, function(subclones){
  drugs <- data[data$Subclone %in% subclones, "Drug_Name"]
  top_drugs <- data.frame(sort(table(drugs), decreasing = TRUE)[1:10])
})

  



### Parameters for customized clustering
# The minimal intersection cutoff for defining the first NMF program in a cluster
min_intersect_initial <- 0.3
# The minimal intersection cutoff for adding a new NMF to the forming cluster
min_intersect_cluster <- 0.3
# The minimal group size to consider for defining the first NMF program in a cluster
min_group_size <- 15

drug_subclone_encoded <- readRDS("drug_subclone_encoded.rds")
sorted_intersection <- sort(apply(similarity_matrix, 2, function(x) {
    (length(which(x >= min_intersect_initial)) - 1)
}), decreasing = TRUE)

# Select top 25 drugs for each cluster
# Every entry contains the NMFs of a chosen cluster
cluster_list <- list()
# Every entry contains the genes of the metaprogram
drug_list <- list()

k <- 1
curr_cluster <- c()
similarity_matrix_original <- similarity_matrix

while (sorted_intersection[1] > min_group_size) {
    curr_cluster <- c(curr_cluster, names(sorted_intersection[1]))

    # Intersection between all remaining subclones and drugs in cluster
    # Genes in the forming cluster are first chosen to be those in the first cluster
    # 'drugs_cluster' always has only 50 genes and evolves during cluster's formation
    drugs_cluster <- subclones_drug_list[[names(sorted_intersection[1])]]

    # Remove selected NMF
    subclones_drug_list[[names(sorted_intersection[1])]] <- NULL

    # Jaccard index between all other subclones and 'drugs_mp'
    subclones <- names(subclones_drug_list)
    intersection_with_drugs_cluster <- c()

    # Compute pairwise Jaccard similarity
    for (i in seq_along(subclones)) {
        intersection_with_drugs_cluster <- c(
            intersection_with_drugs_cluster,
            jaccard_similarity(
                subclones_drug_list[[subclones[i]]],
                drugs_cluster
            )
        )
    }

    names(intersection_with_drugs_cluster) <- subclones
    intersection_with_drugs_cluster <- sort(intersection_with_drugs_cluster, decreasing = TRUE)

    # 'nmf_history' contains drugs from all subclones in the current cluster,
    # for redefining 'drugs_cluster' after adding a new NMF
    subclone_history <- drugs_cluster

    # Create drug list composed of intersecting genes (in descending order by frequency)
    while (intersection_with_drugs_cluster[1] >= min_intersect_cluster) {
        curr_cluster <- c(curr_cluster, names(intersection_with_drugs_cluster)[1])

        # 'drugs_cluster' is newly defined each time according to all subclones in the current cluster
        drugs_cluster_temp <- sort(
            table(c(
                drugs_cluster,
                subclones_drug_list[[names(intersection_with_drugs_cluster)[1]]]
            )),
            decreasing = TRUE
        )

        # Drugs with overlap equal to the 50th gene
        drugs_at_border <- drugs_cluster_temp[which(drugs_cluster_temp == drugs_cluster_temp[20])]

        if (length(drugs_at_border) > 1) {
            # Sort last drugs in genes_at_border according to maximal response level
            # Run across all subclones in curr_cluster and extract response level for each drug
            drugs_curr_response_level <- c()

            for (i in curr_cluster) {
                q <- drug_subclone_encoded[names(drugs_at_border), i]
                drugs_curr_response_level <- c(drugs_curr_response_level, q)
            }

            drugs_curr_response_level_sort <- sort(drugs_curr_response_level, decreasing = TRUE)
            drugs_curr_response_level_sort <- drugs_curr_response_level_sort[unique(
                names(drugs_curr_response_level_sort)
            )]

            drugs_cluster_temp <- c(
                names(drugs_cluster_temp[which(drugs_cluster_temp > drugs_cluster_temp[20])]),
                names(drugs_curr_response_level_sort)
            )
        } else {
            drugs_cluster_temp <- names(drugs_cluster_temp)[1:20]
        }

        subclone_history <- c(
            subclone_history,
            subclones_drug_list[[names(intersection_with_drugs_cluster)[1]]]
        )
        drugs_cluster <- drugs_cluster_temp[1:25]

        # Remove selected NMF
        subclones_drug_list[[names(intersection_with_drugs_cluster)[1]]] <- NULL

        # Intersection between all other subclones and 'drugs_cluster'
        subclones <- names(subclones_drug_list)
        intersection_with_drugs_cluster <- c()

        # Compute pairwise Jaccard similarity
        for (i in subclones) {
            intersection_with_drugs_cluster <- c(
                intersection_with_drugs_cluster,
                jaccard_similarity(
                    subclones_drug_list[[i]],
                    drugs_cluster
                )
            )
        }

        names(intersection_with_drugs_cluster) <- subclones
        intersection_with_drugs_cluster <- sort(intersection_with_drugs_cluster, decreasing = TRUE)
        print(length(subclones_drug_list))
    }

    cluster_list[[paste0("Cluster", k)]] <- curr_cluster
    drug_list[[paste0("Cluster", k)]] <- drugs_cluster
    k <- k + 1

    # Remove current chosen cluster
    similarity_matrix <- similarity_matrix[
        -match(curr_cluster, rownames(similarity_matrix)),
        -match(curr_cluster, colnames(similarity_matrix))
    ]

    # Sort intersection of remaining subclones not included in any of the previous clusters
    sorted_intersection <- sort(apply(similarity_matrix, 2, function(x) {
        (length(which(x >= min_intersect_initial)) - 1)
    }), decreasing = TRUE)

    curr_cluster <- c()
    print(dim(similarity_matrix)[2])
}

saveRDS(cluster_list, "subclone_cluster.rds")
saveRDS(drug_list, "drug_list.rds")



# Plot heatmap with the new clustering
# Reorder MPs to group related ones together and remove MP11, linked only to
# brain metastases of melanoma (with skin pigmentation enrichment)
mp_names <- names(which(unlist(lapply(cluster_list, length)) >= 10))
cluster_list <- cluster_list[mp_names]


#  Sort Jaccard similarity plot according to new clusters:
inds_sorted <- c()
for (j in seq_along(cluster_list)) {
  inds_sorted <- c(inds_sorted, match(
    cluster_list[[j]],
    colnames(similarity_matrix_original)
  ))
}

similarity_matrix_sort <- similarity_matrix_original[inds_sorted, inds_sorted]
rownames(similarity_matrix_sort) <- 1: ncol(similarity_matrix_sort)
colnames(similarity_matrix_sort) <- 1: ncol(similarity_matrix_sort)
labels <- ifelse(1:ncol(similarity_matrix_sort) %in% seq(0, 3000, by = 150), 
                 rownames(similarity_matrix_sort), "")

# Transform cluster list into a vector
cluster_members <- unlist(lapply(seq_along(cluster_list), function(i) {
  group <- names(cluster_list)[i] # Get the name of the current group
  setNames(rep(group, length(cluster_list[[i]])), cluster_list[[i]])
}))

cluster_annotation <- data.frame(
  cluster = cluster_members,
  stringsAsFactors = FALSE
)

cluster_annotation <- cbind(cluster_annotation, subclone_annot_df[rownames(cluster_annotation),])

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = cluster_annotation[names(pals)],
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 14),
  show_legend = c(
    "Sample site" = FALSE
  )
)


# Plot the heatmap
heat <- ComplexHeatmap::Heatmap(
  similarity_matrix_sort,
  top_annotation = top_annotation,
  cluster_rows = FALSE,
  cluster_row_slices = TRUE, 
  row_split =  cluster_annotation$`cluster`,
  cluster_columns = FALSE,
  cluster_column_slices = TRUE,
  show_column_dend = FALSE,
  column_split = cluster_annotation$`cluster`,
  row_labels = labels,
  column_labels = labels,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  row_names_side = "left",
  column_names_side = "bottom",
  row_title = "Subclones",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Subclones",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "Similarity\n(Jaccard index)",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    title_gap = unit(10, "mm")),
  heatmap_width = unit(10, "in"),
  heatmap_height = unit(10, "in")
)

png(
  file = "figures/heatmap_jaccard_clustered.png",
  res = 500,
  width = 16,
  height = 14,
  units = "in"
)
draw(
  heat,
  annotation_legend_side = "top",
  annotation_legend_list = list(tumor_site_legend)
)
dev.off()
