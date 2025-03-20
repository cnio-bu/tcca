library(Seurat)
library(beyondcell)
library(clustree)
library(tidyverse)
library(dplyr)

setwd("/home/lmgonzalezb/Documents/bc-meta/SCellBow/GBM/beyondcell_target/")

# Load normalized expression matrix for nourhan_abdelfattah
expr.matrix <- read.table("../seurat_target/normalized_expression.tsv")
expr.matrix <- as.matrix(expr.matrix)

# Generate geneset object with one of the ready to use signature collections
gs <- GetCollection(SSc, include.pathways = FALSE)

# Compute BCS for the SSc. This might take a few minutes depending on the size of your dataset.
bc <- bcScore(expr.matrix, gs, expr.thres = 0.1) 

saveRDS(bc, "bc_gbm_nourhan_abdelfattah.rds")

bc <- readRDS("bc_gbm_nourhan_abdelfattah.rds")

# Run the UMAP reduction. 
bc <- bcUMAP(bc, k.neighbors = 4, res = 0.2)
# Run the bcUMAP function again, specifying the number of principal components you want to use.
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.2)

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", factor.col = TRUE, pt.size = 1.5)

# Visualize whether the cells are clustered based on the number of genes detected per each cell.
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)

bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")
bc <- bcRegressOut(bc, vars.to.regress = c("nFeature_RNA"))

# Recompute the UMAP.
bc <- bcUMAP(bc, pc = 10, k.neighbors = 20, res = 0.1)
# Visualize the UMAP.
bcClusters(bc, UMAP = "beyondcell", idents = "nFeature_RNA", factor.col = FALSE, pt.size = 1.5)
# Visualize the therapeutic clusters.
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1.5)


# Load VAR scores
var_scores <- read.table("/home/lmgonzalezb/Documents/bc-meta/var_score/var_scores.tsv")
bc@meta.data <- merge(bc@meta.data, var_scores, by = "row.names", all = FALSE)
rownames(bc@meta.data) <- bc@meta.data$Row.names
bc@meta.data$Row.names <- NULL

# Group cells into clusters based on score quartiles
# Define the function to categorize scores
categorize_quartiles <- function(seu_metadata, score_name) {
  # Calculate the quartiles
  quartiles <- quantile(seu_metadata[[score_name]],
                        robs = c(0.25, 0.5, 0.75)
  )
  # Define a function to categorize scores based on quartiles
  categorize_score <- function(score) {
    if (score <= quartiles[2]) {
      return("L")
    } else if (score <= quartiles[3]) {
      return("IL")
    } else if (score <= quartiles[4]) {
      return("IH")
    } else {
      return("H")
    }
  }
  
  # Apply the categorization function to create a new column
  seu_metadata[[paste0(score_name, "_cat")]] <- sapply(
    seu_metadata[[score_name]],
    categorize_score
  )
  return(seu_metadata)
}


# Apply the categorization function to create new columns
bc@meta.data <- categorize_quartiles(bc@meta.data, "A_scaled")
bc@meta.data <- categorize_quartiles(bc@meta.data, "VR_scaled")
bc@meta.data <- categorize_quartiles(bc@meta.data, "VARscore_scaled_subs")
bc@meta.data <- categorize_quartiles(bc@meta.data, "VARscore_scaled_sum")

# Visualize A score.
bcClusters(bc, UMAP = "beyondcell", idents = "A_scaled", factor.col = FALSE, pt.size = 1.5)
bcClusters(bc, UMAP = "beyondcell", idents = "A_scaled_cat", factor.col = TRUE, pt.size = 1.5)
# bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.2", factor.col = TRUE, pt.size = 1.5)

# Obtain unextended therapeutic cluster-based statistics.
bc <- bcRanks(bc, idents = "A_scaled_cat", extended = FALSE)

bc4Squares(bc, idents = "A_scaled_cat", lvl = "H", top = 3)

vln.plot1 <- ggplot(bc@meta.data, aes(x = bc_clusters_res.0.1, y = A_scaled, color = bc_clusters_res.0.1)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "A score across Therapeutic Clusters", x = "Therapeutic clusters", y = "A score") +
  theme_bw()
vln.plot1

ggsave("plots/vlnplot_Ascore.png", plot = vln.plot1, width = 10, height = 5)

bc@meta.data$A_scaled_cat <- factor(bc@meta.data$A_scaled_cat, levels = c("H", "L", "IH", "IL"))
barplot <- ggplot(bc@meta.data, aes(x = bc_clusters_res.0.1, fill = A_scaled_cat)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("H" = "#f8766d",  # Blue  
                               "L" = "#7cae00",  # Green  
                               "IH" = "#00bfc4", # Red  
                               "IL" = "#c77cff")) +
  labs(x = "Therapeutic Clusters", y = "Cell fraction", fill = "A score groups") +
  ggtitle("A score groups across therapeutic clusters") +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
ggsave("plots/stacked_barplot_Ascore.png", plot = barplot, width = 7, height = 6)

FindDrugs(bc, "AA-COCF3")
drugs <- bcSignatures(bc, UMAP = "beyondcell", 
                      signatures = list(values = c("sig-21191", "sig-21040", "sig-21049")), 
                      pt.size = 1.5)
for (i in 1:length(drugs)){
  ggsave(paste0("plots/drug", i, ".png"), plot = drugs[[i]], width = 7, height = 6, dpi = 500)
}

