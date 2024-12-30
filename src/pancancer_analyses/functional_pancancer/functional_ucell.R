library(Seurat)
library(BPCells)
library(UCell)
library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional")
##---------------------------------Functions ---------------------------------##
read.gmt <- function(gmt_file) {
  sigs_list <- list()
  sigs <- scan(gmt_file, what = character(), sep = "\n")
  for (sig in sigs) {
    sig <- unlist(strsplit(sig, "\t"))
    sig <- unique(sig[nzchar(sig)])
    sigs_list[[sig[1]]] <- sig[3:length(sig)]
  }
  return(sigs_list)
}

# Read the gmt file with the 6 drug resistance mechanisms gene sets.
gmt_list <- read.gmt(
  "/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt"
)

# Join bidirectional signatures in the same vector adding '+' to the genes in the 
# upregulared gene set and '-' to the genes in the downregulated gene set.
sig_names <- sub("_(UP|DOWN|DN)$", "", names(gmt_list))
bisigs <- sig_names[duplicated(sig_names)]

for (sig in unique(sig_names)) {
  up <- paste0(sig, "_UP")
  down <-  if (paste0(sig, "_DOWN") %in% names(gmt_list))
    paste0(sig, "_DOWN")
  else
    paste0(sig, "_DN")
  
  if (sig %in% bisigs) {
    gmt_list[[up]] <- paste0(gmt_list[[up]], "+")
    gmt_list[[down]] <- paste0(gmt_list[[down]], "-")
    gmt_list[[sig]] <- c(gmt_list[[up]], gmt_list[[down]])
    gmt_list[[up]] <- NULL
    gmt_list[[down]] <- NULL
    
  } else{
    if (up %in% names(gmt_list)) {
      gmt_list[[up]] <- paste0(gmt_list[[up]], "+")
      gmt_list[[sig]] <- gmt_list[[up]]
      gmt_list[[up]] <- NULL
    }
    
    if (down %in% names(gmt_list)) {
      gmt_list[[down]] <- paste0(gmt_list[[down]], "-")
      gmt_list[[sig]] <- gmt_list[[down]]
      gmt_list[[down]] <- NULL
    }
  }
}


# Load level 2 Seurat object
seu <- readRDS(
  "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2.rds"
)
malignant <- subset(seu, subset = malignancy == TRUE)
colnames(malignant) <- paste0("c", c(1:ncol(malignant)))

# Load the beyondcell pancancer object to select only cells where BCS was computed
bc <- readRDS(
  "/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/results/beyondcell_pancancer_final_res.Rds"
)
cells <- colnames(bc)

# Subset the malignant cells with BCS score from the counts slot
malignant <- subset(malignant, cells = cells)

seu_ucell <- AddModuleScore_UCell(malignant, features = gmt_list)

# Save seurat object with UCell scores in metadata
saveRDS(
  seu_ucell,
  "/storage/scratch01/users/mgonzalezb/bc-meta/functional/seurat_ucell.rds"
)

# Save UCell scores separately as a matrix
sig_names <- paste0(names(gmt_list), "_UCell")
scores_ucell <- t(seu_ucell@meta.data[, sig_names])
scores_ucell <- as(scores_ucell, "sparseMatrix")
write_matrix_dir(mat = scores_ucell,
                 dir = "/storage/scratch01/users/mgonzalezb/bc-meta/functional/full_mat_ucell",
                 overwrite = TRUE)



## Plot UCell enrichment scores per Therapeutic Cluster
seu_ucell <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/functional/seurat_ucell.rds")
seu_ucell$Therapeutic_clusters <- bc@meta.data$therapeutic_clusters_k.300.res.0.5
scores_ucell <- as.matrix(seu_ucell@meta.data[, grepl("UCell", colnames(seu_ucell@meta.data))])

# keep signatures with high variance across cells
top_rv <- matrixStats::colVars(scores_ucell)
top_median <- median(top_rv)
top_rv <- top_rv[top_rv >= top_median]

scores_ucell <- scores_ucell[, names(top_rv)]
scores_ucell <- as.data.frame(scale(x = scores_ucell, center = TRUE, scale = TRUE))

metadata <- cbind(seu_ucell@meta.data[, "Therapeutic_clusters", drop = FALSE], scores_ucell)
long_metadata <- metadata %>%
  pivot_longer(cols = ends_with("UCell"),
               names_to = "GeneSet",
               values_to = "UCell_score") %>%
  as.data.frame()

# Group by Cluster and GeneSet, then calculate mean scores
cluster_means <- long_metadata %>%
  group_by(Therapeutic_clusters, GeneSet) %>%
  summarize(MeanScore = mean(UCell_score, na.rm = TRUE),
            .groups = "drop")

# Create bubble plot
bubble_plot <- ggplot(cluster_means, aes(x = GeneSet, y = Therapeutic_clusters)) +
  geom_point(aes(color = MeanScore), size = 5) + # Map size and color to scores
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) + # Gradient
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.margin = unit(c(1, 1, 1, 3), "cm"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(
      size = 12,
      face = "bold",
      margin = margin(
        t = 0,
        r = 20,
        b = 0,
        l = 0
      )
    ),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Mean Ucell scores per TC",
    x = "Functional gene sets",
    y = "TC",
    color = "UCell score"
  )

ggsave(
  bubble_plot,
  file = "results_ucell/bubble_mean_ucell.png",
  dpi = 500,
  width = 17,
  height = 10
)


# Compute UCell scores for meta-programs obtained by NMF
setwd("/home/lmgonzalezb/Documents/bc-meta/functional_mps/")
mp_list <- readRDS("mp_list_allsamples.rds")
print(mp_list)
mp_ucell <- AddModuleScore_UCell(malignant, features = mp_list)

# Save seurat object with UCell scores in metadata
saveRDS(mp_ucell, "seurat_mps_ucell.rds")

# Plot UCell scores 
mp_ucell <- readRDS("seurat_mps_ucell.rds")
mp_ucell$Therapeutic_clusters <- bc@meta.data$therapeutic_clusters_k.300.res.0.5
scores_ucell <- mp_ucell@meta.data[, grepl("UCell", colnames(mp_ucell@meta.data))]
colnames(scores_ucell) <- gsub("P_", "P", colnames(scores_ucell))
mps <- c("MP1", "MP8", "MP14", "MP3", "MP13", "MP2", "MP9", "MP4", "MP10", "MP12", 
         "MP5", "MP6",  "MP7")
select_mps <- match(paste0(mps, "_UCell"), colnames(scores_ucell))
scores_ucell <- as.matrix(scores_ucell[, select_mps])
scores_ucell <- as.data.frame(scale(x = scores_ucell, center = TRUE, scale = TRUE))
metadata <- cbind(mp_ucell@meta.data[, "Therapeutic_clusters", drop = FALSE], scores_ucell)
long_metadata <- metadata %>%
  pivot_longer(cols = ends_with("UCell"),
               names_to = "GeneSet",
               values_to = "UCell_score") %>%
  as.data.frame() %>%
  mutate(GeneSet = factor(GeneSet, levels = rev(paste0(mps, "_UCell"))))

# Group by Cluster and GeneSet, then calculate mean scores
cluster_means <- long_metadata %>%
  group_by(Therapeutic_clusters, GeneSet) %>%
  summarize(MeanScore = mean(UCell_score, na.rm = TRUE),
            .groups = "drop")

# Create bubble plot
bubble_plot <- ggplot(cluster_means, aes(x = Therapeutic_clusters, y = GeneSet)) +
  geom_point(aes(color = MeanScore), size = 5) + # Map size and color to scores
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    oob = scales::squish,
    limits = c(-0.5, 0.5)
  ) + # Gradient
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.margin = unit(c(1, 1, 1, 3), "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      margin = margin(
        t = 0,
        r = 20,
        b = 0,
        l = 0
      )
    ),
    plot.title = element_text(face = "bold", size = 15)
  ) +
  labs(
    title = "Mean Ucell scores per TC",
    x = "TC",
    y = "Meta-programs",
    color = "UCell score"
  )

ggsave(
  bubble_plot,
  file = "bubble_mean_ucell.pdf",
  width = 8,
  height = 8
)
