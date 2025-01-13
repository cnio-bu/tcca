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
seu <- readRDS("./single_cell/seurat/malignant/aml_sander_lambo.rds")


setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell_vs_perception/aml")
# Perform Seurat analysis of malignant cell population
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 100)
png("elbow.png", width = 10, height = 10, units = "in", res = 300)
ElbowPlot(seu, ndims = 100)
dev.off()

seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = seq(0.1, 1, 0.1))

seu <- RunUMAP(seu, dims = 1:50)
png("umap_patient.png", width = 10, height = 10, units = "in", res = 300)
DimPlot(seu, group.by = "patient") + NoLegend()
dev.off()


umap_treatment <- DimPlot(seu, reduction = "umap", group.by = "Treatment_Outcome")
ggsave(
  "plots/umap_treatment.png",
  umap_treatment,
  dpi = 300,
  height = 7,
  width = 7
)

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
DimPlot(seu, reduction = "umap.harmony", group.by = "patient") + NoLegend()

clustree <- clustree(seu@meta.data[, grep("RNA_snn_res.", colnames(seu@meta.data))], 
                     prefix = "RNA_snn_res.")
ggsave("plots/clustree_integration.png", plot = clustree, dpi = 300, height = 7, width = 7)
umap_patient <- DimPlot(seu, reduction = "umap.harmony", group.by = "patient") + NoLegend()
ggsave(
  "plots/umap_patient_integrated.png",
  umap_patient,
  dpi = 300,
  height = 7,
  width = 7
)

umap_treatment <- DimPlot(seu, reduction = "umap.harmony", group.by = "Treatment_Outcome")
ggsave(
  "plots/umap_treatment_outcome_integrated.png",
  umap_treatment,
  dpi = 300,
  height = 7,
  width = 7
)

seu$RNA_snn_res.0.4 <- factor(seu$RNA_snn_res.0.4, levels = 0:15)
umap_clusters <- DimPlot(seu, reduction = "umap.harmony", group.by = "RNA_snn_res.0.4", label = TRUE)
ggsave(
  "plots/clusters_integrated.png",
  umap_clusters,
  dpi = 300,
  height = 7,
  width = 7
)


saveRDS(seu, "integrated_aml_sander_lambo.rds")


# Rename cell names in seu to match with Beyondcell and Perception matrices
colnames(seu) <- paste0("c", seq_along(colnames(seu)))

# # Load beyondcell object of the same study
bc <- open_matrix_dir("/storage/scratch01/shared/projects/bc-meta/beyondcell/studywise_bpcells/aml_sander_lambo")
colnames(bc) <- paste0("c", seq_along(colnames(bc)))

# # Load perception object of the same study
pc <- open_matrix_dir("/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception/v5/aml_sander_lambo")
colnames(pc) <- paste0("c", seq_along(colnames(pc)))

# Common drugs between beyondcell aand perception
drug_models <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception/drug_models/FDA_approved_drugs_models.rds")
rownames(pc) <- names(drug_models)

bc <- bc[, colnames(seu)]
pc <- pc[, colnames(seu)]

load("drugInfo.RData")
bc_drug.names <- drugInfo$Synonyms[drugInfo$Synonyms$IDs %in% rownames(bc), "drugs"]
bc_drug.names <- tolower(unique(bc_drug.names))

# Select drugs used in ArmB and ArmC: ADE (cytarabine (Ara-C), daunorubicin hydrochloride, 
# and etoposide phosphate) + asparaginase + mitoxantrone + bortezomib/sorafenib 
drug_oi <- c("cytarabine", "daunorubicin", "etoposide", "mitoxantrone", "bortezomib", "sorafenib")
common_drugs <- intersect(drug_oi, rownames(pc))

# # Plot scores for drugs of interest
# Compute beyondcell switch point for those drugs:
bc <- as.matrix(bc)
scaled.matrix <- t(apply(bc, 1, scales::rescale, to = c(0, 1)))

drug_ids <- drugInfo$Synonyms %>%
  filter(drugs %in% toupper(drug_oi)) %>%
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
drugs <- data.table::fread("../../reference/final_moas - Collapsed.tsv") %>%
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
seu$Sample_type <- factor(sub("^[^_]*_([^_]*)_.*$", "\\1",seu$sample), levels = c("Dx", "Rem", "Rel"))


umap_sample_type <- DimPlot(seu, reduction = "umap.harmony", group.by = "Sample_type")
ggsave(
  "plots/umap_sample_type_integrated.png",
  umap_sample_type,
  dpi = 300,
  height = 7,
  width = 7
)

umap_sample_type <- DimPlot(seu, reduction = "umap", group.by = "Sample_type")
ggsave(
  "plots/umap_sample_type.png",
  umap_sample_type,
  dpi = 300,
  height = 7,
  width = 7
)

plot_bcs <-function(seurat_object, drug_ids, split.by, filename_prefix, png_height, png_width) {
  lapply(drug_ids, function(drug_id) {
    drug_name <- drugs[drug_id, "preferred.drug.names"]
    bcs_umap <- FeaturePlot(
      object = seurat_object,
      features = drug_id,
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
      theme(legend.position = "right",
            legend.title = element_text(size = 8, margin = margin(b = 7)),
            legend.margin = margin(l = 5),
            legend.text = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8), 
            axis.text.y = element_text(size = 8), 
            axis.line = element_line(size = 0.3))
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

plot_bcs(seurat_object = seu, 
         drug_ids = drug_ids,
         split.by = "Sample_type", 
         filename_prefix = "plots/beyondcell/split_sample_type/bcs_",
         png_height = 4, 
         png_width = 12)

# For bortezomid and sorafenib, we should separate cells from different arm treatments
# to study drug sensitivity to this drugs (these drugs are the difference between arms)
drug_armb <- drugs$IDs[drugs$preferred.drug.names == "BORTEZOMIB"] 
drug_armc <- drugs$IDs[drugs$preferred.drug.names == "SORAFENIB"] 

clinical <- read.xlsx("suplemmentary_aml_sander_lambo.xlsx", sheet = "Single cell Cohort") %>%
  mutate(patient = Patient.identifier) %>%
  select(patient, Treatment) %>%
  distinct()

seu@meta.data <- seu@meta.data %>%
  left_join(clinical, by = "patient")
rownames(seu@meta.data) <- paste0("c", seq_along(colnames(seu)))


seu_armb <- subset(seu, subset = Treatment == "Arm B")
seu_armc <- subset(seu, subset = Treatment == "Arm C")

# Plot Bortezomib and Sorafenib drug sensitivities for the corresponding subset of cells.
plot_bcs(seurat_object = seu_armb, 
         drug_ids = drug_armb,
         split.by = "Sample_type", 
         filename_prefix = "plots/beyondcell/split_sample_type/bcs_armB_",
         png_height = 4, 
         png_width = 12)
plot_bcs(seurat_object = seu_armc, 
         drug_ids = drug_armc,
         split.by = "Sample_type", 
         filename_prefix = "plots/beyondcell/split_sample_type/bcs_armC_",
         png_height = 4, 
         png_width = 12)  

## Analyze the sensitivity to all drugs splitting by treatment response
# First, subset Diagnosis samples
seu$Treatment_Outcome <- factor(seu$Treatment_Outcome, levels = c("Censored", "Relapsed"))
seu_dx <- subset(seu, subset = Sample_type == "Dx")
plot_bcs(seurat_object = seu_dx, 
         drug_ids = drug_ids,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Dx_treatment_response/bcs_",
         png_height = 4, 
         png_width = 8)


# Plot Bortezomib and Sorafenib drug sensitivities for the corresponding subset of cells.
seu_dx_armb <- subset(seu, subset = Sample_type == "Dx" & Treatment == "Arm B")
seu_dx_armc <- subset(seu, subset = Sample_type == "Dx" & Treatment == "Arm C")
plot_bcs(seurat_object = seu_dx_armb, 
         drug_ids = drug_armb,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Dx_treatment_response/bcs_armB_",
         png_height = 4, 
         png_width = 8)
plot_bcs(seurat_object = seu_dx_armc, 
         drug_ids = drug_armc,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Dx_treatment_response/bcs_armC_",
         png_height = 4, 
         png_width = 8)  

# Second, subset remission samples and anlyze drug response in patient with censored 
# and relapsed treatment outcomes
seu_rem <- subset(seu, subset = Sample_type == "Rem")
plot_bcs(seurat_object = seu_rem, 
         drug_ids = drug_ids,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Rem_treatment_response/bcs_",
         png_height = 4, 
         png_width = 8)

# Plot Bortezomib and Sorafenib drug sensitivities for the corresponding subset of cells.
seu_rem_armb <- subset(seu, subset = Sample_type == "Rem" & Treatment == "Arm B")
seu_rem_armc <- subset(seu, subset = Sample_type == "Rem" & Treatment == "Arm C")
plot_bcs(seurat_object = seu_rem_armb, 
         drug_ids = drug_armb,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Rem_treatment_response/bcs_armB_",
         png_height = 4, 
         png_width = 4)
plot_bcs(seurat_object = seu_rem_armc, 
         drug_ids = drug_armc,
         split.by = "Treatment_Outcome", 
         filename_prefix = "plots/beyondcell/split_Rem_treatment_response/bcs_armC_",
         png_height = 4, 
         png_width = 4)  




# Get the plots of PERCEPTION predicted killing scores of drugs shared by Beyondcell and PERCEPTION
pc <- as.matrix(pc)
scaled.matrix_pc <- t(apply(pc, 1, function(drug) {
  scales::rescale(rank(-as.numeric(drug)), to = c(0, 1))
}))
colnames(scaled.matrix_pc) <- colnames(pc)
pc_subset <- as.data.frame(t(scaled.matrix_pc[common_drugs, ]))
colnames(pc_subset) <- paste0(colnames(pc_subset), "_", "perception")

seu@meta.data <- cbind(seu@meta.data, pc_subset)

plot_perception <- function(seurat_object, split.by, filename_prefix, png_height, png_width){
  lapply(colnames(pc_subset), function(drug_name) {
    perception_umap <- FeaturePlot(
      object = seurat_object,
      features = drug_name,
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
      theme(legend.position = "right",
            legend.title = element_text(size = 8, margin = margin(b = 7)),
            legend.margin = margin(l = 5),
            legend.text = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            axis.text.x = element_text(size = 8), 
            axis.text.y = element_text(size = 8), 
            axis.line = element_line(size = 0.3))
    
    ggsave(
      plot = perception_umap,
      filename = paste0(filename_prefix, gsub("_perception", "", drug_name), ".png"),
      dpi = 300,
      height = png_height,
      width = png_width
    )
})
}

# Plot only for edaunorubicin and etoposide (the only drugs from ArmB and ArmC 
# available in PERCEPTION models)
plot_perception(seurat_object = seu, 
                split.by = "Sample_type",
                filename_prefix = "plots/perception/split_sample_type/killing_score_",
                png_height = 4, 
                png_width = 13)

## Analyze the sensitivity to the two drugs splitting by treatment response
# First, subset Diagnosis samples
seu_dx <- subset(seu, subset = Sample_type == "Dx")
plot_perception(seurat_object = seu_dx, 
                split.by = "Treatment_Outcome", 
                filename_prefix = "plots/perception/split_Dx_treatment_response/killing_score_",
                png_height = 4,
                png_width = 9)

# Second, subset remission samples and anlyze drug response in patient with censored 
# and relapsed treatment outcomes
seu_rem <- subset(seu, subset = Sample_type == "Rem")
plot_perception(seurat_object = seu_rem, 
                split.by = "Treatment_Outcome", 
                filename_prefix = "plots/perception/split_Rem_treatment_response/killing_score_",
                png_height = 4, 
                png_width = 9)

# Plot boxplots for quantificantion of Beyondcell and PERCEPTION scores differences
# across diagnostics, remission and relapse samples.
long_metadata <- melt(seu@meta.data, 
                      id.vars = "Sample_type", 
                      measure.vars = c("sig-21177", 
                                     "daunorubicin_perception",
                                     "sig-21164", 
                                     "etoposide_perception"),
                      variable.name = "Sensitivity",
                      value.name = "Score")

boxplot_bysample <- ggplot(long_metadata, aes(x = Sample_type, y = Score, fill = Sample_type)) +
                           geom_boxplot() +
                           facet_wrap(~ Sensitivity, scales = "free_y") +
                           stat_compare_means(label = "p.signif",
                                              comparisons = list(c("Dx", "Rem"),
                                                                  c("Rem", "Rel"),
                                                                  c("Dx", "Rel")),
                                              method = "t.test", 
                                              label.y = seq(0.9, 1.1, length.out = 5),
                                              label.x.npc = "center",
                                              tip.length = 0.01,
                                              vjust = 0.03,
                                              size = 4) +
                           labs(x = "Sample type", 
                                y = "Score") +
                           ylim(0, 1.5) +
                           theme_minimal() +
                           theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                 axis.text.y = element_text(size = 10),
                                 axis.title.x = element_text(size = 12),
                                 axis.title.y = element_text(size = 12))                      
ggsave("plots/quantification_common_drugs/boxplot_by_sample_type.png",
       plot = boxplot_bysample,
       dpi = 300,
       height = 12,
       width = 10)  



#Plot boxplots for quantificantion of Beyondcell and PERCEPTION scores differences
# between long-term remission patients and relapsed patient in diagnostic and relapsed samples.

# Reshape data: Wide to long format
long_metadata <- melt(seu_dx@meta.data, 
                      id.vars = "Treatment_Outcome", 
                      measure.vars = c("sig-21177", 
                                     "daunorubicin_perception",
                                     "sig-21164", 
                                     "etoposide_perception"),
                      variable.name = "Sensitivity",
                      value.name = "Score")

boxplot_dx <- ggplot(long_metadata, aes(x = Treatment_Outcome, y = Score, fill = Treatment_Outcome)) +
              geom_boxplot() +
              facet_wrap(~ Sensitivity, scales = "free_y") + # Separate plots for each signature
              stat_compare_means(label = "p.signif",
                                 comparisons = list(c("Censored", "Relapsed")),
                                 method = "t.test", 
                                 label.x.npc = "center",
                                 tip.length = 0.01,
                                 vjust = 0.03,
                                 size = 4) +
              labs(x = "Treatment Outcome", 
                  y = "Score") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12))                      
ggsave("plots/quantification_common_drugs/boxplot_diagnosis.png",
       plot = boxplot_dx,
       dpi = 300,
       height = 12,
       width = 10)  


long_metadata <- melt(seu_rem@meta.data, 
                      id.vars = "Treatment_Outcome", 
                      measure.vars = c("sig-21177", 
                                     "daunorubicin_perception",
                                     "sig-21164", 
                                     "etoposide_perception"),
                      variable.name = "Sensitivity",
                      value.name = "Score")

boxplot_dx <- ggplot(long_metadata, aes(x = Treatment_Outcome, y = Score, fill = Treatment_Outcome)) +
              geom_boxplot() +
              facet_wrap(~ Sensitivity, scales = "free_y") + # Separate plots for each signature
              stat_compare_means(label = "p.signif",
                                 comparisons = list(c("Censored", "Relapsed")),
                                 method = "t.test", 
                                 label.x.npc = "center",
                                 tip.length = 0.01,
                                 vjust = 0.03,
                                 size = 4) +
              labs(x = "Treatment Outcome", 
                  y = "Score") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12))                      
ggsave("plots/quantification_common_drugs/boxplot_remission.png",
       plot = boxplot_dx,
       dpi = 300,
       height = 12,
       width = 10)  