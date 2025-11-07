library(BPCells)
library(ComplexHeatmap)
library(tidyverse)

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception/")

## Source TCCA palette
source(file = "/storage/scratch01/users/mgonzalezb/bc-meta/TCCA_palette.R")

seu <- readRDS("results/perception_pancancer_FDA_approved.Rds")
sketched_mat <- open_matrix_dir(dir = "results/sketch_mat_perception_5k/")

sketched_mat <- as.matrix(sketched_mat)
scaled_mat <- t(apply(sketched_mat, 1, function(drug) {
  scales::rescale(rank(-as.numeric(drug)), to = c(0, 1))
}))

colnames(scaled_mat) <- colnames(sketched_mat)

# Clusters based on perception killing scores
perception_clusters <- seu@meta.data %>%
  filter(new_cell_id %in% colnames(sketched_mat))
colnames(perception_clusters) <- gsub("therapeutic", "perception", colnames(perception_clusters))

## Therapeutic clusters from beyondcell
tcs <- data.table::fread("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/results/tcs.tsv") %>%
  filter(new_cell_id %in% colnames(sketched_mat)) %>%
  select(c(1:32, 58)) %>%
  left_join(
    perception_clusters %>% select(new_cell_id, perception_clusters_k.300.res.0.2),
    by = "new_cell_id"
  )

## human readable origins
translat_human_sites <- c(
  adrenal_gland = "Adrenal gland",
  bladder = "Bladder",
  bone_marrow = "Bone marrow",
  brain = "Brain",
  breast = "Breast",
  colon = "Colon",
  skin = "Skin",
  esophagus = "Esophagus",
  oesophagus = "Esophagus",
  kidney = "Kidney",
  liver = "Liver",
  lung = "Lung",
  lymph_node = "Lymph node",
  other = "Other",
  ovary = "Ovary",
  pancreas = "Pancreas",
  prostate = "Prostate",
  soft_tissue = "Soft tissue"
)


clinical_features <- tcs %>%
  mutate(
    summarised_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "other"
    ),
    adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
    is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
    treated = ifelse(treated, "Treated", "Untreated"),
    sex = ifelse(sex == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
    therapeutic_clusters_k.300.res.0.5 = as_factor(therapeutic_clusters_k.300.res.0.5),
    perception_clusters_k.300.res.0.2 = factor(perception_clusters_k.300.res.0.2, 
                                               levels = as.character(seq(0, 5, 1)))
  )

cells_annot_df <- clinical_features %>%
  select(
    new_cell_id,
    sex,
    adult_pediatric,
    is_blood,
    sample_type,
    summarised_tumor_site,
    treated,
    therapeutic_clusters_k.300.res.0.5,
    perception_clusters_k.300.res.0.2
  ) %>%
  as.data.frame()

cells_annot_df$summarised_tumor_site <-  translat_human_sites[cells_annot_df$summarised_tumor_site]

rownames(cells_annot_df) <- cells_annot_df$new_cell_id
cells_annot_df$new_cell_id <- NULL

colnames(cells_annot_df) <- c(
  "Sex",
  "Age group",
  "Solid/Liquid",
  "Tumor type",
  "Tumor site",
  "Treatment",
  "Beyondcell Cluster",
  "Perception Cluster"
)

pcs_colors <- c(
  "5" = "#369CBB",
  "4" = "#406792",
  "3" = "#B46F9C",
  "2" = "#D05B61",
  "1" = "#FE7B47",
  "0" = "#FFA72C"
)
pals <- list(
  "Sex" = sex_colors,
  "Age group" = age_colors,
  "Solid/Liquid" = sl_colors,
  "Tumor type" = pm_colors,
  "Tumor site" = tumor_sites_colors,
  "Treatment" = treatment_colors,
  "Beyondcell Cluster" = tcs_colors,
  "Perception Cluster" = pcs_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df =  cells_annot_df,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_legend = c("Tumor site" = FALSE)
)


## get drug names
drugs <- data.table::fread(
  "/storage/scratch01/users/mgonzalezb/bc-meta/reference/final_moas - Collapsed.tsv"
) %>%
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

doi <- data.frame(drug_names = toupper(rownames(sketched_mat)))
MoAs <- doi %>%
  left_join(drugs, by = c("drug_names" = "preferred.drug.names")) %>%
  select(-IDs) %>%
  distinct() %>%
  mutate(
    MoAs = case_when(
      drug_names == "ABEMACICLIB" ~ "Cell cycle arrest",
      drug_names == "HOMOHARRINGTONINE" ~ "Protein synthesis inhibitor",
      drug_names == "NIRAPARIB" ~ "DNA related agent",
      drug_names == "PONATINIB" ~ "BCR-ABL inhibitor",
      drug_names == "THIOGUANINE" ~ "DNA related agent",
      drug_names == "VINDESINE" ~ "Microtubule agent",
      drug_names == "VINFLUNINE" ~ "Microtubule agent",
      TRUE ~ collapsed.MoAs
    )
  ) %>%
  select(MoAs)

colnames(MoAs) <- "Mechanism of action"

moa_pals <- list("Mechanism of action" = MoAs_colors)


right_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = MoAs,
  which = "row",
  col = moa_pals,
  show_annotation_name = FALSE
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Tumor site`),
  legend_gp = gpar(fill = pals$`Tumor site`),
  ncol = 2,
  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Tumor site"
)

png(
  file = "./results/figures/sketched_perception_with_tcs_killing_score.png",
  res = 300,
  width = 20,
  height = 14,
  units = "in"
)


heat <- ComplexHeatmap::Heatmap(
  mat = scaled_mat,
  #mat = t(sketched_mat),
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  row_split = 5,
  row_title = NULL,
  column_order = rownames(cells_annot_df[order(cells_annot_df$`Perception Cluster`), ]),
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  column_split =  cells_annot_df$`Perception Cluster`,
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  row_labels = rownames(sketched_mat),
  show_row_names = TRUE,
  column_names_rot = 0,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(title = "Killing score"),
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(8, "in")
)

ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
  "HEATMAP_LEGEND_PADDING" = unit(2, "cm"),
  legend_gap = unit(1, "cm")
)
draw(
  heat,
  annotation_legend_side = "top",
  annotation_legend_list = list(tumor_site_legend)
)
dev.off()


## Generate UMAP plot from raw UMAP
tcs_umap <- DimPlot(object = seu,
                    group.by = "therapeutic_clusters_k.300.res.0.2",
                    reduction = "umap")

tcs_umap_clean <- tcs_umap +
  ggtitle("") +
  scale_color_manual(name = "Perception cluster", values = pcs_colors) +
  scale_shape_manual() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

tcs_umap_clean$layers[[1]]$aes_params$size <- 0.1
tcs_umap_clean$layers[[1]]$aes_params$alpha <- 0.7

ggsave(
  plot = tcs_umap_clean,
  filename = "results/figures/therapeutic_clusters_umap_k300_res0.2.pdf",
  dpi = 300,
  height = 7,
  width = 7
)

# Plot distribution of beyondcell clusters into perception clusters and viceversa
colnames(seu@meta.data) <- gsub("therapeutic", "perception", colnames(seu@meta.data))
merged_all_df <- data.table::fread(
  "/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/results/tcs.tsv"
  ) %>%
  select(new_cell_id, therapeutic_clusters_k.300.res.0.5) %>%
  left_join(seu@meta.data %>% select(new_cell_id, perception_clusters_k.300.res.0.2),
            by = "new_cell_id") %>%
  mutate(
    beyondcell_clusters = factor(therapeutic_clusters_k.300.res.0.5, levels = seq(0, 4, 1)),
    perception_clusters = factor(perception_clusters_k.300.res.0.2, levels = seq(0, 5, 1))
  )


perception_barplot <- ggplot(merged_all_df,
                             aes(x = perception_clusters, fill = beyondcell_clusters)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tcs_colors) +
  theme_bw()

ggsave(
  "results/figures/perception_barplot_tcs.png",
  perception_barplot,
  width = 8,
  height = 8,
  dpi = 500
)


beyondcell_barplot <- ggplot(merged_all_df,
                             aes(x = beyondcell_clusters, fill = perception_clusters)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = pcs_colors) +
  theme_bw()

ggsave(
  "results/figures/beyondcell_barplot_tcs.png",
  beyondcell_barplot,
  width = 8,
  height = 8,
  dpi = 500
)
