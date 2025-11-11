library(BPCells)
library(Seurat)
library(clustree)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(reshape2)
library(ggpubr)


setwd("/storage/scratch01/shared/projects/bc-meta")

# Load seurat object with malignant cells of nsclc_stefan_salcher study
seu <- readRDS("./single_cell/seurat/malignant/nsclc_stefan_salcher.rds")
seu <- merge(seu[[1]], seu[c(2:length(seu))])

# Select only cells from study He_Fan_2021
seu <- subset(seu, subset = study == "He_Fan_2021")

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell_vs_perception/luad")
# Perform Seurat analysis of malignant cell population
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 100)
png(
  "plots/elbow.png",
  width = 10,
  height = 10,
  units = "in",
  res = 300
)
ElbowPlot(seu, ndims = 100)
dev.off()

seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = seq(0.1, 1, 0.1))

seu <- RunUMAP(seu, dims = 1:50)

# We need to perform integration before
seu@meta.data$patient_sample <- paste0(seu@meta.data$patient, "_", seu@meta.data$sample)

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient_sample)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

seu <- IntegrateLayers(
  object = seu,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony"
)

seu <- JoinLayers(seu)
seu <- FindNeighbors(seu, reduction = "integrated.harmony", dims = 1:50)
seu <- FindClusters(seu, resolution = seq(0.1, 1, 0.1))

seu <- RunUMAP(
  seu,
  reduction = "integrated.harmony",
  dims = 1:50,
  reduction.name = "umap.harmony"
)

# clustree <- clustree(seu@meta.data[, grep("integrated_snn_res.", colnames(seu@meta.data))],
#                      prefix = "integrated_snn_res.")
# ggsave("plots/clustree_integration.png", plot = clustree, dpi = 300, height = 7, width = 7)
# umap_patient <- DimPlot(seu, reduction = "umap.harmony", group.by = "patient") + NoLegend()
# ggsave(
#   "plots/umap_patient_integrated.png",
#   umap_patient,
#   dpi = 300,
#   height = 7,
#   width = 7
# )

# umap_treatment <- DimPlot(seu, reduction = "umap.harmony", group.by = "Treatment_Outcome")
# ggsave(
#   "plots/umap_treatment_outcome_integrated.png",
#   umap_treatment,
#   dpi = 300,
#   height = 7,
#   width = 7
# )

# seu$RNA_snn_res.0.4 <- factor(seu$RNA_snn_res.0.4, levels = 0:15)
# umap_clusters <- DimPlot(seu, reduction = "umap.harmony", group.by = "RNA_snn_res.0.4", label = TRUE)
# ggsave(
#   "plots/clusters_integrated.png",
#   umap_clusters,
#   dpi = 300,
#   height = 7,
#   width = 7
# )


saveRDS(seu, "integrated_luad_stefan_salcher_subset.rds")



setwd("/home/lmgonzalezb/Documents/bc-meta/luad_beyondcell_perception/")
#-------------------------------------------------------------------------------
### Check sensitivity to everolimus, sunitinib, ibrutinib and erlotinib
#-------------------------------------------------------------------------------
seu <- readRDS("integrated_luad_stefan_salcher_subset.rds")

# # Load beyondcell object of the same study
bc <- open_matrix_dir("beyondcell_nsclc_stefan_salcher/")
bc <- bc[, colnames(seu)]

# # Load perception object of the same study
pc <- open_matrix_dir("perception_nsclc_stefan_salcher/")
pc <- pc[, colnames(seu)]

# Rename cell names in seu to match with Beyondcell and Perception matrices
colnames(seu) <- paste0("c", seq_along(colnames(seu)))
colnames(bc) <- paste0("c", seq_along(colnames(bc)))
colnames(pc) <- paste0("c", seq_along(colnames(pc)))

# Common drugs between beyondcell aand perception
drug_models <- readRDS("../pdac_beyondcell_perception/new_drug_models.rds")
rownames(pc) <- names(drug_models)

# Subset cells from sample LUAD1
#seu <- subset(seu, subset = sample == "He_Fan_2021_LUAD1")
bc <- bc[, colnames(seu)]
pc <- pc[, colnames(seu)]

bc <- as.matrix(bc)
bc[is.na(bc)] <- 0
bc <- as.data.frame(t(as.matrix(bc)))
bc$patient <- seu$patient

bc <- bc %>%
  pivot_longer(cols = -patient,
               names_to = "Drug",
               values_to = "BCS")
bc_mean <- bc %>%
  group_by(patient, Drug) %>%
  summarize(
    Mean_BCS_patient = mean(BCS, na.rm = TRUE),
    Min_BCS_patient = min(BCS, na.rm = TRUE),
    Max_BCS_patient = max(BCS, na.rm = TRUE)
  )
load("../pdac_beyondcell_perception/drugInfo.RData")

# Select drugs used in ArmB and ArmC: ADE (cytarabine (Ara-C), daunorubicin hydrochloride,
# and etoposide phosphate) + asparaginase + mitoxantrone + bortezomib/sorafenib
drug_oi <- c("afatinib", "dacomitinib", "gefitinib")

# # Plot scores for drugs of interest
# Compute beyondcell switch point for those drugs:
bc <- as.matrix(bc)
bc[is.na(bc)] <- 0
scaled.matrix <- t(apply(bc, 1, scales::rescale, to = c(0, 1)))

drug_ids <- drugInfo$Synonyms %>%
  filter(drugs %in% toupper(drug_oi)) %>%
  pull(IDs)
drug_ids <- intersect(drug_ids, rownames(bc))
drug_ids <- c("sig-20987",
              "sig-21058",
              "sig-21185",
              "sig-21213",
              "sig-21408")
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
drugs <- data.table::fread("../bc-meta_repo/bc-meta/reference/final_moas - Collapsed.tsv") %>%
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

# Join normalized BCS scores and scaled scores of each drug as columns with Seurat metadata
scaled_bcs <- as.data.frame(t(scaled.matrix[drug_ids, ]))
colnames(scaled_bcs) <- paste0(colnames(scaled_bcs), "_scaled")
normalized_bcs <- t(bc[drug_ids, ])
colnames(normalized_bcs) <- paste0(colnames(normalized_bcs), "_normalized")
seu@meta.data <- cbind(seu@meta.data, cbind(scaled_bcs, normalized_bcs))


plot_bcs <- function(seurat_object,
                     features,
                     split.by,
                     filename_prefix,
                     png_height,
                     png_width) {
  lapply(features, function(feature) {
    drug_id <- strsplit(feature, split = "_")[[1]][1]
    drug_name <- drugs[drug_id, "preferred.drug.names"]
    bcs_umap <- FeaturePlot(
      object = seurat_object,
      features = feature,
      split.by = split.by,
      alpha = 1,
      pt.size = 0.8,
      reduction = "umap.harmony",
      raster = FALSE
    ) &
      scale_colour_gradientn(
        colors = c("#1D61F2", "#83A8F7", "#F7F7F7", "#FF9CBB", "#DA0078"),
        values = c(0, sp[drug_id], 1),
        limits = c(0, 1),
        na.value = "grey50",
        guide = guide_colorbar()
      ) &
      labs(
        title = drug_name,
        color = "Scaled BCS",
        x = "UMAP1",
        y = "UMAP2"
      ) &
      theme(
        legend.position = "right",
        legend.title = element_text(size = 8, margin = margin(b = 7)),
        legend.margin = margin(l = 5),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(size = 0.3)
      )
    ggsave(
      plot = bcs_umap,
      filename = paste0(
        filename_prefix,
        drug_id,
        "_",
        sub(" .*", "", drug_name),
        ".png"
      ),
      dpi = 300,
      height = png_height,
      width = png_width
    )
  })
}

plot_bcs(
  seurat_object = seu,
  features = paste0(drug_ids, "_scaled"),
  split.by = NULL,
  filename_prefix = "plots/beyondcell/bcs_luad_",
  png_height = 8,
  png_width = 10
)


# Get the plots of PERCEPTION predicted killing scores
pc <- as.matrix(pc)
scaled.matrix_pc <- t(apply(pc, 1, function(drug) {
  scales::rescale(rank(-as.numeric(drug)), to = c(0, 1))
}))
colnames(scaled.matrix_pc) <- colnames(pc)

# Join normalized BCS scores and scaled scores of each drug as columns with Seurat metadata
scaled_pcs <- as.data.frame(t(scaled.matrix_pc[drug_oi, ]))
colnames(scaled_pcs) <- paste0(colnames(scaled_pcs), "_scaled")
normalized_pcs <- apply(pc, 1, function(drug) {
  scale(-as.numeric(drug), center = TRUE, scale = TRUE)
})
normalized_pcs <- normalized_pcs[, drug_oi]
colnames(normalized_pcs) <- paste0(colnames(normalized_pcs), "_normalized")
seu@meta.data <- cbind(seu@meta.data, cbind(scaled_pcs, normalized_pcs))

plot_perception <- function(seurat_object,
                            features,
                            split.by,
                            filename_prefix,
                            png_height,
                            png_width) {
  lapply(features, function(feature) {
    drug_name <- strsplit(feature, split = "_")[[1]][1]
    perception_umap <- FeaturePlot(
      object = seurat_object,
      features = feature,
      split.by = split.by,
      alpha = 1,
      pt.size = 0.8,
      reduction = "umap.harmony",
      raster = FALSE
    ) &
      scale_colour_gradientn(
        colors = c("#00BFC4", 'lightgrey' , "#F8766D"),
        #values = c(0, 0.5, 1),
        limits = c(0, 1),
        na.value = "grey50",
        guide = "colourbar"
      ) &
      labs(
        title = toupper(gsub("_.*", "", drug_name)),
        color = "Scaled Killing Score",
        x = "UMAP1",
        y = "UMAP2"
      ) &
      theme(
        legend.position = "right",
        legend.title = element_text(size = 8, margin = margin(b = 7)),
        legend.margin = margin(l = 5),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(size = 0.3)
      )
    
    ggsave(
      plot = perception_umap,
      filename = paste0(
        filename_prefix,
        gsub("_perception", "", drug_name),
        ".png"
      ),
      dpi = 300,
      height = png_height,
      width = png_width
    )
  })
}

plot_perception(
  seurat_object = seu,
  features = paste0(drug_oi, "_scaled"),
  split.by = NULL,
  filename_prefix = "plots/perception/pcs_luad1",
  png_height = 8,
  png_width = 10
)


# Plot boxplots for quantificantion of Beyondcell and PERCEPTION scores differences
# across diagnostics, remission and relapse samples.
bc_drugs <- grep("^sig-.*_normalized", colnames(seu@meta.data), value = TRUE)
metadata <- seu@meta.data[, bc_drugs]
drug_names <- drugs[gsub("_normalized", "", bc_drugs), "preferred.drug.names"]
colnames(metadata) <- paste0(drug_names, "_", gsub("_normalized", "", bc_drugs))

comparisons <- combn(colnames(metadata), 2, simplify = FALSE)
long_df <- pivot_longer(metadata,
                        cols = everything(),
                        names_to = "Drug",
                        values_to = "Score")


metadata <- seu@meta.data[, c("patient", bc_drugs)]
colnames(metadata) <- gsub("-", "_", colnames(metadata))
metadata <- metadata %>%
  group_by(patient) %>%
  summarise(
    sig_21267 = mean(sig_21267_normalized, na.rm = TRUE),
    sig_21266 = mean(sig_21266_normalized, na.rm = TRUE),
    sig_20887 = mean(sig_20887_normalized, na.rm = TRUE),
    sig_21280 = mean(sig_21280_normalized, na.rm = TRUE),
    sig_21279 = mean(sig_21279_normalized, na.rm = TRUE),
    sig_21227 = mean(sig_21227_normalized, na.rm = TRUE)
  )
# %>%  column_to_rownames(var = "patient")
colnames(metadata) <- paste0(colnames(metadata), "_", drugs[gsub("_", "-", colnames(metadata)), "preferred.drug.names"])
long_data <- metadata %>%
  pivot_longer(
    cols = -patient,
    # Select drug columns
    names_to = "drug",
    # Name for the new column
    values_to = "value"
  )
long_data$drug <- paste0(long_data$drug, "_", drugs[gsub("_", "-", long_data$drug), "preferred.drug.names"])
(
  plot_line <- ggplot(long_data, aes(
    x = drug,
    y = value,
    group = patient,
    color = patient
  )) +
    geom_line() +              # Connect points with lines
    geom_point(size = 4) +     # Add points
    theme_minimal() +          # Use a minimal theme
    ylim(-0.025, 0.025) +
    labs(title = "Drug Response by Patient", x = "Drug", y = "BCS") +
    scale_color_discrete(name = "Patient") +
    theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12
    ))
)

ggsave(
  "plots/point_beyondcell_per_patient.png",
  plot = plot_line,
  dpi = 300,
  height = 4,
  width = 6
)

# Create boxplot
(
  boxplot_by_tool <- ggplot(long_df, aes(x = Drug, y = Score)) +
    geom_boxplot() +
    # stat_compare_means(comparisons = comparisons,
    #                    label = "p.signif",
    #                    method = "t.test",
    #                    label.y = seq(50, 60, length.out = 21),
    #                    label.x.npc = "center",
    #                    tip.length = 0.01,
    #                    vjust = 0.03,
    #                    size = 4) +
    labs(x = "Drug", y = "Score") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
)

ggsave(
  "plots/boxplot_beyondcell_all_sample.png",
  plot = boxplot_by_tool,
  dpi = 300,
  height = 7,
  width = 6
)

compare_means(Score ~ Drug, data = long_df, method = "t.test")

# Boxplot with PERCEPTION scores
metadata <- seu@meta.data[, c("patient", paste0(drug_oi, "_", "normalized"))]
metadata <- metadata %>%
  group_by(patient) %>%
  summarise(
    afatinib = mean(afatinib_normalized, na.rm = TRUE),
    dacomitinib = mean(dacomitinib_normalized, na.rm = TRUE),
    gefitinib = mean(gefitinib_normalized, na.rm = TRUE)
  )
long_data <- metadata %>%
  pivot_longer(
    cols = -patient,
    # Select drug columns
    names_to = "drug",
    # Name for the new column
    values_to = "value"
  )
(
  plot_line <- ggplot(long_data, aes(
    x = drug,
    y = value,
    group = patient,
    color = patient
  )) +
    geom_line() +              # Connect points with lines
    geom_point(size = 4) +     # Add points
    theme_minimal() +          # Use a minimal theme
    labs(title = "Drug Response by Patient", x = "Drug", y = "Response Value") +
    scale_color_discrete(name = "Patient")
)

ggsave(
  "plots/point_perception_per_patient.png",
  plot = plot_line,
  dpi = 300,
  height = 4,
  width = 6
)

comparisons <- combn(colnames(metadata), 2, simplify = FALSE)
long_df <- pivot_longer(metadata,
                        cols = everything(),
                        names_to = "Drug",
                        values_to = "Score")

# Create boxplot
(
  boxplot_by_tool <- ggplot(long_df, aes(x = Drug, y = Score)) +
    geom_boxplot() +
    # stat_compare_means(comparisons = comparisons,
    #                    label = "p.signif",
    #                    method = "t.test",
    #                    label.y = seq(1.2, 1.3, length.out = 3),
    #                    label.x.npc = "center",
    #                    tip.length = 0.01,
    #                    vjust = 0.03,
    #                    size = 4) +
    labs(x = "Drug", y = "Score") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
)

ggsave(
  "plots/boxplot_perception_per_sample.png",
  plot = boxplot_by_tool,
  dpi = 300,
  height = 7,
  width = 6
)


compare_means(Score ~ Drug, data = long_df, method = "t.test")

corr_matrix = cbind(
  bulk = unlist(lapply(drug_models, function(x)
    x$performance_in_bulk[2])),
  pseudo_bulk = unlist(lapply(drug_models, function(x)
    x$performance_in_pseudo_bulk[2])),
  sc = unlist(lapply(drug_models, function(x)
    x$performance_in_scRNA[2]))
)
rownames(corr_matrix) = names(drug_models)

pval_matrix = cbind(
  bulk = unlist(lapply(drug_models, function(x)
    x$performance_in_bulk[1])),
  pseudo_bulk = unlist(lapply(drug_models, function(x)
    x$performance_in_pseudo_bulk[1])),
  sc = unlist(lapply(drug_models, function(x)
    x$performance_in_scRNA[1]))
)
rownames(pval_matrix) = names(drug_models)
#Number of cell lines with AUC values used in validation
TotalcellLines_used_in_validation = unlist(lapply(drug_models, function(x)
  length(
    na.omit(x$predVSgroundTruth$pred_gt_mscRNA$Observed)
  )))

significant_drugs_names = names(which(corr_matrix[, 3] > 0.3 &
                                        pval_matrix[, 3] < 0.05))
