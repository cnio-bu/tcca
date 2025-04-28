library(dplyr)
library(tidyverse)
library(uwot)
set.seed(123)
setwd("/home/lmgonzalezb/Documents/bc-meta/sctherapy/")
source(file = "../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

# Load the drug prediction per subclone
data <- read.table("full_table_drug_prediction.tsv")

# Load the jaccard similarity matrix between subclones
similarity_matrix <- readRDS("jaccard_matrix.rds")

# Load the clustering results of spectral clustering
cluster_assignment <- readRDS("speclustering_reordered.rds")


# Plot the clusters on the umap projection from the distance matrix
dist_matrix <- as.dist(1 - similarity_matrix)
umap <- umap(dist_matrix)
umap_df <- as.data.frame(umap)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$cluster <- cluster_assignment[rownames(umap_df)]

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(alpha = 0.8, size = 3) +
  theme_minimal() +
  labs(title = "UMAP Projection of Points Based on Distance Matrix",
       x = "UMAP 1",
       y = "UMAP 2")


# Plot the clusters on the umap projection from the binary matrix of drug predictions
subclone_drug_df <- data[c("Subclone", "Drug_Name")]
binary_df <- subclone_drug_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Drug_Name, values_from = value, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "Subclone")
umap <- umap(binary_mat)
umap_df <- as.data.frame(umap)
colnames(umap_df) <- c("UMAP1", "UMAP2")
rownames(umap_df) <- rownames(binary_df)
umap_df$cluster <- cluster_assignment[rownames(umap_df)]

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(values = sctherapy_colors) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line()) +
  labs(x = "UMAP 1",
       y = "UMAP 2", 
       colour = "Cluster")
