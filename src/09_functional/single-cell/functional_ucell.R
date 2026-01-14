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

# Read the gmt file
gmt_list <- read.gmt(
  "/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt"
)

# Join bidirectional signatures in the same vector adding '+' to the genes in the upregulared gene set
# and '-' to the genes in the downregulated gene set.
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
setwd("/storage/scratch01/shared/projects/bc-meta/")
seu <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2.rds")
malignant <- subset(seu, subset = malignancy == TRUE)

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
  pivot_longer(
    cols = ends_with("UCell"),
    # Adjust this to match your score column names
    names_to = "GeneSet",
    # New column for gene set names
    values_to = "UCell_score"             # New column for scores
  ) %>%
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
mp_list <- readRDS("functional_nmf/metaprograms/mp_list.rds")
print(mp_list)
mp_ucell <- AddModuleScore_UCell(malignant, features = mp_list)

# Save seurat object with UCell scores in metadata
saveRDS(mp_ucell, "functional_nmf/seurat_mps_ucell.rds")
write.table(
  seu@meta.data[, 
                grep("UCell", colnames(seu@meta.data), value = TRUE)], 
            "functional_nmf/mps_ucell_scores.tsv"
  )