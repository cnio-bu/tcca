library(Seurat)
library(BPCells)
library(ComplexHeatmap)
library(tidyverse)

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/")

## Source TCCA palette
source(file = "../TCCA_palette.R")
source(file = "Format.R")
top_rv <- readRDS("results/top_rv.Rds")
top_rv <- names(top_rv)
seu <- readRDS("results/beyondcell_pancancer_final_res.Rds")

DefaultAssay(seu) <- "sketch_50k"
full_mat <- as.matrix(seu[["sketch_50k"]]$data)
scaled.matrix <- t(apply(full_mat, 1, scales::rescale, to = c(0, 1)))
seu[["sketch_50k"]]$scale.data <- round(scaled.matrix, digits = 2)

# Plot umap with TCs
tcs_umap <- DimPlot(object = seu,
                    group.by = "therapeutic_clusters_k.300.res.0.5",
                    reduction = "umap")

tcs_umap_clean <- tcs_umap +
  ggtitle("") +
  scale_color_manual(name = "Therapeutic cluster", values = tcs_colors) +
  scale_shape_manual() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(title = "THERAPEUTIC CLUSTERS") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

tcs_umap_clean$layers[[1]]$aes_params$size <- 1
tcs_umap_clean$layers[[1]]$aes_params$alpha <- 0.7

ggsave(
  plot = tcs_umap_clean,
  filename = "results/figures/therapeutic_clusters_umap_all_k300_res0.5.png",
  dpi = 300,
  height = 7,
  width = 7
)

# Subset drugs of interest to compute the SP
seu <- subset(seu, features = top_rv)

# Compute beyondcell switch point for those drugs:
bcs_norm <- as.matrix(seu[["sketch_50k"]]$data)

sp <- lapply(top_rv, function(sig) {
  m <- bcs_norm[sig, ]
  if (any(m == 0)) {
    sp <- rep(which(m == 0)[1], times = 2)
  } else {
    lower.bound <- which(m == max(m[m <= 0]))[1]
    upper.bound <- which(m == min(m[m >= 0]))[1]
    sp <- c(lower.bound, upper.bound)
  }
  sp_scaled <- round(sum(seu[["sketch_50k"]]$scale.data[sig, sp]) / 2, digits = 2)
  return(sp_scaled)
})

sp <- unlist(sp)
names(sp) <- top_rv

## get drug names
drugs <- data.table::fread("../reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs) %>%
  mutate(
    collapsed.MoAs = case_when(
      preferred.drug.names == "VANDETANIB" ~ "VEGFR inhibitor",
      preferred.drug.names == "DASATINIB" ~ "Kinase inhibitor",
      preferred.drug.names == "RIGOSERTIB" ~ "Other",
      preferred.drug.names == "SORAFENIB" ~ "Kinase inhibitor",
      TRUE ~ collapsed.MoAs
    )
  ) %>%
  distinct() %>%
  as.data.frame()

drugs <- drugs[drugs$IDs %in% top_rv, ]
rownames(drugs) <- drugs$IDs

lapply(top_rv, function(drug_id) {
  drug_name <- drugs[drug_id, "preferred.drug.names"]
  bcs_umap <- FeaturePlot(
    object = seu,
    features = drug_id,
    alpha = 1,
    pt.size = 1,
    slot = "scale.data"
  ) +
    scale_colour_gradientn(
      colors = c("#1D61F2", "#83A8F7", "#F7F7F7", "#FF9CBB", "#DA0078"),
      values = c(0, sp[drug_id], 1),
      na.value = "grey50",
      guide = "colourbar"
    ) +
    labs(
      title = drug_name,
      color = "Scaled BCS",
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme(legend.position = "right",
          legend.title = element_text(margin = margin(b = 10)))
  
  ggsave(
    plot = bcs_umap,
    filename = paste0(
      "results/figures/scaled_bcs/bcs_",
      sub(" .*", "", drug_name),
      ".png"
    ),
    dpi = 300,
    height = 7,
    width = 7
  )
})
