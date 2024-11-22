library(BPCells)
library(Seurat)
library(Beyondcell)
library(clustree)
library(tidyverse)

setwd("/home/lmgonzalezb/Documents/bc-meta/nsclc_beyondcell_perception/")

# Load seurat object with malignant cells of nsclc_stefan_salcher study
seu <- readRDS("malignant_nsclc_stefan_salcher.rds")
seu <- merge(seu[[1]], seu[c(2:length(seu))])

# Keep only malignant cells remaining in the filtered level 2 object
tcca.metadata <- read.table("/home/lmgonzalezb/Documents/bc-meta/tcca_annotation_raw.tsv",
                            header = TRUE)
tcca.subset <- tcca.metadata[tcca.metadata$study == "nsclc_stefan_salcher" &
                               tcca.metadata$malignancy == TRUE, ]
seu <- RenameCells(seu, new.names = paste0(colnames(seu), "_", "22"))

seu <- subset(seu, cells = tcca.subset$cell)

# Load beyondcell object of the same study
bc <- open_matrix_dir("beyondcell_nsclc_stefan_salcher")

# Load perception object of the same study
pc <- open_matrix_dir("perception_nsclc_stefan_salcher")


# Perform Seurat analysis of malignant cell population
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 100)
ElbowPlot(seu, ndims = 100)


seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = seq(0.1, 1, 0.1))

seu <- RunUMAP(seu, dims = 1:50)
DimPlot(seu, group.by = "patient") + NoLegend()

# We need to perform integration before
seu@meta.data$patient_sample <- paste0(seu@meta.data$patient, "_", seu@meta.data$sample)

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient_sample)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

seu <- IntegrateLayers(
  object = seu,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca"
)

seu <- JoinLayers(seu)
seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:50)
seu <- FindClusters(seu, resolution = seq(0.1, 1, 0.1))

seu <- RunUMAP(
  seu,
  reduction = "integrated.cca",
  dims = 1:50,
  reduction.name = "umap.cca"
)
DimPlot(seu, reduction = "umap.cca", group.by = "patient") + NoLegend()

clustree(seu@meta.data[, grep("RNA_snn_res.", colnames(seu@meta.data))], prefix = "RNA_snn_res.")

clusters <- DimPlot(seu, reduction = "umap.cca", group.by = "RNA_snn_res.0.3")
ggsave(
  "nsclc_umap_clusters.png",
  clusters,
  dpi = 300,
  height = 7,
  width = 7
)
saveRDS(seu, "malignant_integrated_nsclc_stefan_salcher.rds")

# Common drugs between beyondcell aand perception
colnames(bc) <- paste0(colnames(bc), "_", "22")
colnames(pc) <- paste0(colnames(pc), "_", "22")
drug_models <- readRDS("FDA_approved_drugs_models.rds")
rownames(pc) <- names(drug_models)

bc <- bc[, colnames(seu)]
pc <- pc[, colnames(seu)]

bc_drug.names <- drugInfo$Synonyms[drugInfo$Synonyms$IDs %in% rownames(bc), "drugs"]
bc_drug.names <- tolower(unique(bc_drug.names))

common_drugs <- intersect(bc_drug.names, rownames(pc))

# Plot scores for drugs of interest
# Compute beyondcell switch point for those drugs:
bc <- as.matrix(bc)
scaled.matrix <- t(apply(bc, 1, scales::rescale, to = c(0, 1)))

drug_ids <- drugInfo$Synonyms %>%
  filter(drugs %in% toupper(common_drugs)) %>%
  pull(IDs)
drug_ids <- intersect(drug_ids, rownames(bc))

sp <- lapply(drug_ids, function(sig) {
  m <- bc[sig, ]
  if (any(m == 0)) {
    sp <- rep(which(m == 0)[1], times = 2)
  } else {
    lower.bound <- which(m == max(m[m <= 0]))[1]
    upper.bound <- which(m == min(m[m >= 0]))[1]
    sp <- c(lower.bound, upper.bound)
  }
  sp_scaled <- round(sum(scaled.matrix[sig, sp]) / 2, digits = 2)
  return(sp_scaled)
})

sp <- unlist(sp)
names(sp) <- drug_ids

## get drug names
drugs <- data.table::fread("final_moas - Collapsed.tsv") %>%
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

drugs <- drugs[drugs$IDs %in% drug_ids, ]
rownames(drugs) <- drugs$IDs

# Join BCS scores of each drug as columns in Seurat metadata
bc_subset <- as.data.frame(t(scaled.matrix[drug_ids, ]))
seu@meta.data <- cbind(seu@meta.data, bc_subset)
lapply(drug_ids, function(drug_id) {
  drug_name <- drugs[drug_id, "preferred.drug.names"]
  bcs_umap <- FeaturePlot(
    object = seu,
    features = drug_id,
    alpha = 1,
    pt.size = 1,
    reduction = "umap.cca"
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
      "figures/beyondcell/bcs_",
      drug_id,
      "_",
      sub(" .*", "", drug_name),
      ".png"
    ),
    dpi = 300,
    height = 7,
    width = 7
  )
})


# Get the plots of PERCEPTION predicted killing scores of drugs shared by Beyondcell and PERCEPTION
pc <- as.matrix(pc)
scaled.matrix_pc <- t(apply(pc, 1, function(drug) {
  scales::rescale(rank(-as.numeric(drug)), to = c(0, 1))
}))
colnames(scaled.matrix_pc) <- colnames(pc)
pc_subset <- as.data.frame(t(scaled.matrix_pc[common_drugs, ]))
colnames(pc_subset) <- paste0(colnames(pc_subset), "_", "perception")

seu@meta.data <- cbind(seu@meta.data, pc_subset)

lapply(colnames(pc_subset), function(drug_name) {
  perception_umap <- FeaturePlot(
    object = seu,
    features = drug_name,
    alpha = 1,
    pt.size = 1,
    reduction = "umap.cca"
  ) +
    scale_colour_gradientn(
      colors = c("#00BFC4", 'lightgrey' , "#F8766D"),
      #values = c(0, 0.5, 1),
      na.value = "grey50",
      guide = "colourbar"
    ) +
    labs(
      title = toupper(gsub("_.*", "", drug_name)),
      color = "Scaled Killing Score",
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme(legend.position = "right",
          legend.title = element_text(margin = margin(b = 10)))
  
  ggsave(
    plot = perception_umap,
    filename = paste0("figures/perception/bcs_", drug_name, ".png"),
    dpi = 300,
    height = 7,
    width = 7
  )
})
