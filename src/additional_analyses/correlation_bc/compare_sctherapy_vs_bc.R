library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta")

## Subclone-level annotation of TCCA cells
metadata <- read.table("./single_cell/seurat/tcca/tcca_metadata.tsv",
                       sep = "\t",
                       header = TRUE)
metadata <- metadata %>%
  rownames_to_column("original_barcode")

clonality <- read.table("./single_cell/cna_metadata/full_clonality_table_lvl2.tsv",
                        row.names = NULL)
duplicated <- clonality[clonality$original_barcode %in% clonality$original_barcode
                        [duplicated(clonality$original_barcode)], ]
clonality <- clonality %>%
  filter(!(original_barcode %in% duplicated$original_barcode)) %>% # Remove 16 cells with duplicated names
  dplyr::select(original_barcode, scevan_prediction)

full_metadata <- left_join(metadata, clonality, by = c("cell" = "original_barcode")) %>%
  filter(malignancy == "True")

## Add subclonal metadata info to beyondcell object
options(future.globals.maxSize = 20 * 1024^3)
options(Seurat.object.assay.version = 'v5')

# Load full beyondcell mat
mat <- open_matrix_dir(dir = "beyondcell/full_mat_beyondcell")
colnames(mat) <- paste0("c", c(1:ncol(mat)))

full_metadata$new_cell_id <-  paste0("c", c(1:nrow(full_metadata)))
rownames(full_metadata) <- full_metadata$new_cell_id

bc <- Seurat::CreateSeuratObject(
  counts = mat,
  assay = "RNA",
  project = "beyondcell_pancancer",
  meta.data = full_metadata
)


## Keep only cells in agreement (the ones used for scTherapy prediction)
bc <- subset(
  bc,
  subset = (malignancy == "False" & scevan_prediction == "normal") |
    (malignancy == "True" &
       scevan_prediction == "tumor")
)

normalized_bcs <- as(bc[["RNA"]]$counts, "sparseMatrix")

## Compute mean BCS per sample and subclone and residual's mean per subclone.
sigs <- rownames(normalized_bcs)
normalized.long <- normalized_bcs %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cells") %>%
  tidyr::pivot_longer(
    cols = all_of(sigs),
    names_to = "IDs",
    values_to = "enrichment",
    values_drop_na = FALSE
  )

# Add study, patient, sample, scevan_subclone
bc@meta.data <- bc@meta.data %>%
  mutate(study_sample = paste(study, sample, sep = "_"))

normalized.long <- normalized.long %>%
  left_join(select(
    bc@meta.data,
    c(
      "study",
      "patient",
      "sample",
      "study_sample",
      "scevan_subclone",
      "new_cell_id"
    )
  ), by = c("cells" = "new_cell_id"))
stats.long <- normalized.long %>%
  dplyr::group_by(IDs, study_sample) %>%
  dplyr::mutate(mean_sample = round(mean(enrichment, na.omit = TRUE), digits = 2)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(IDs, scevan_subclone) %>%
  dplyr::mutate(mean_subclone = round(mean(enrichment, na.omit = TRUE), digits = 2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(resid = enrichment - mean_sample) %>%
  dplyr::group_by(IDs, scevan_subclone) %>%
  dplyr::mutate(residuals.mean_subclone = round(mean(resid, na.rm = TRUE), digits = 2)) %>%
  dplyr::ungroup()

mean_bcs_subclone <- stats.long %>%
  dplyr::select(scevan_subclone, IDs, mean_subclone) %>%
  dplyr::distinct() %>% # para evitar duplicados
  tidyr::pivot_wider(names_from = scevan_subclone, values_from = mean_subclone) %>%
  column_to_rownames(var = "IDs")

# Guardar como tsv la matrix de drugs x subclones
write.table(
  mean_bcs_subclone,
  "single_cell/sctherapy_int_bc/mean_bcs_per_subclone.tsv",
  row.names = TRUE
)
mean_bcs_subclone <- read.table("/storage/scratch01/users/mgonzalezb/mean_bcs_per_subclone.tsv",
                                header = TRUE)

# Load scTherapy predictions per subclone and TC annotation
sctherapy <- read.table("single_cell/sctherapy/results/subclone_level_annotated.tsv",
                        header = TRUE)

drugs_sctherapy <- unique(sctherapy$Drug.Name)
drugs_bc <- data.table::fread("./reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs, studies) %>%
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
  as.data.frame() %>%
  filter(IDs %in% rownames(bc))


# Common drugs between scTherapy and Beyondcell
common_drugs <- intersect(drugs_sctherapy, drugs_bc$preferred.drug.names)
# For duplicated drugs (several IDs corresponding to the same drug) select based on study
drugs_bc_final <- drugs_bc %>%
  filter(IDs %in% rownames(bc)) %>%
  mutate(study_priority = factor(studies, levels = c("GDSC", "CTRP", "PRISM"))) %>%
  arrange(preferred.drug.names, study_priority) %>%
  # For duplicated drugs take the first one
  group_by(preferred.drug.names) %>%
  slice(1) %>%
  ungroup() %>%
  select(IDs, preferred.drug.names, collapsed.MoAs, studies)

# Remove duplicated drugs in BCS matrix of subclones
mean_bcs_subclone <- mean_bcs_subclone[rownames(mean_bcs_subclone) %in% drugs_bc_final$IDs, ]
rownames(mean_bcs_subclone) <- drugs_bc_final %>%
  filter(IDs %in% rownames(mean_bcs_subclone)) %>%
  arrange(match(IDs, rownames(mean_bcs_subclone))) %>%
  pull(preferred.drug.names)
#sctherapy$Subclone.Name <- gsub("-| ", ".", sctherapy$Subclone.Name)
mean_bcs_subclone_subset <- mean_bcs_subclone[common_drugs, unique(sctherapy$Subclone.Name)] %>%
  rownames_to_column("Drug.Name") %>%
  pivot_longer(cols = -Drug.Name,
               names_to = "Subclone.Name",
               values_to = "bcs")

bc_sctherapy <- mean_bcs_subclone_subset %>%
  left_join(select(sctherapy, c("Subclone.Name", "Response", "Drug.Name")), by = c("Drug.Name", "Subclone.Name"))

# Boxplot to visualize BCS score based on scTherapy predicted response
sctherapy_cluster <- sctherapy %>%
  select(Subclone.Name, scTherapy.Cluster) %>%
  distinct()
bc_sctherapy <- bc_sctherapy %>%
  mutate(
    Response = factor(Response, levels = c("High", "High-to-moderate", "Moderate")),
    Response = fct_na_value_to_level(Response, level = "Not predicted")
  ) %>%
  left_join(sctherapy_cluster, by = "Subclone.Name")

plot <- ggplot(bc_sctherapy, aes(x = Response, y = bcs, fill = Response)) +
  geom_boxplot() +
  #facet_wrap(~scTherapy.Cluster, scales = "free_y") + # un panel por cluster
  theme_bw() +
  labs(title = "Beyondcell scores por Subclone-Droga y Response", x = "scTherapy Response", y = "Beyondcell Score (bcs)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(
  "single_cell/sctherapy_int_bc/bcs_subclones.png",
  plot = plot,
  width = 15,
  height = 15
)

# Codify response of scTherapy as 0-3
bc_sctherapy <- bc_sctherapy %>%
  mutate(
    Response_num = case_when(
      Response == "High" ~ 3,
      Response == "High-to-moderate" ~ 2,
      Response == "Moderate" ~ 1,
      Response == "Not predicted" ~ 0
    )
  )

#
# Beyondcell: subclones en filas, drogas en columnas
bc_mat <- bc_sctherapy %>%
  select(Subclone.Name, Drug.Name, bcs) %>%
  pivot_wider(names_from = Drug.Name, values_from = bcs) %>%
  column_to_rownames("Subclone.Name") %>%
  as.matrix()

# scTherapy: subclones en filas, drogas en columnas
sct_mat <- bc_sctherapy %>%
  select(Subclone.Name, Drug.Name, Response_num) %>%
  pivot_wider(names_from = Drug.Name, values_from = Response_num) %>%
  column_to_rownames("Subclone.Name") %>%
  as.matrix()

drug_names <- colnames(bc_mat)
cor_mat <- matrix(
  NA,
  nrow = length(drug_names),
  ncol = length(drug_names),
  dimnames = list(drug_names, drug_names)
)

# Calcular correlación Spearman cruzada
for (i in seq_along(drug_names)) {
  for (j in seq_along(drug_names)) {
    cor_mat[i, j] <- cor(sct_mat[, i], bc_mat[, j], method = "spearman", use = "pairwise.complete.obs")
  }
}

library(pheatmap)
png(
  "single_cell/sctherapy_int_bc/correlation_plot.png",
  height = 15,
  width = 15,
  units = "in",
  res = 300
)
pheatmap(
  cor_mat,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Correlación cruzada: scTherapy vs Beyondcell",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  number_format = "%.2f"
)
dev.off()


## Comparison between Beyondcell and scTherapy at cell level.
source(file = "~/bc-meta/TCCA_palette.R")
library(ComplexHeatmap)
# First we need to create a Seurat object with BCS for cells included in
# scTherapy analysis
bc <- subset(bc, subset = scevan_subclone %in% unique(sctherapy$Subclone.Name))
sctherapy_subclone_lvl <- sctherapy %>%
  select(Subclone.Name, scTherapy.Cluster) %>%
  distinct()
bc@meta.data <- bc@meta.data %>%
  left_join(sctherapy_subclone_lvl,
            by = c("scevan_subclone" = "Subclone.Name"))
rownames(bc@meta.data) <- bc@meta.data$new_cell_id
bc$scTherapy.Cluster <- factor(bc$scTherapy.Cluster, levels = as.character(c(1:10)))
saveRDS(bc,
        "single_cell/sctherapy_int_bc/bc_obj_samecells_sctherapy.rds")

# Sketch of 5000 cells per scTherapy cluster
bc_c1 <- subset(bc, subset = scTherapy.Cluster == 1)
bc_c1[["RNA"]]$data <- bc_c1[["RNA"]]$counts

bc_c1 <- Seurat::SketchData(
  object = bc_c1,
  assay = "RNA",
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_5k"
)

sketched_mat <- bc_c1[["sketch_5k"]]$data
sketched_mat <- as.matrix(sketched_mat)
sketched_mat <- scale(x = sketched_mat, center = TRUE, scale = TRUE)

# Select drug names
sketched_mat_unique <- sketched_mat[drugs_bc_final$IDs, ]
rownames(sketched_mat_unique) <- drugs_bc_final$preferred.drug.names
sketched_mat_commondrugs <- sketched_mat_unique[common_drugs, ]

# Prepare data for heatmap
translat_human_sites <- c(
  "adrenal_gland" = "Adrenal gland",
  "bladder" = "Bladder",
  "bone_marrow" = "Bone marrow",
  "brain" = "Brain",
  "breast" = "Breast",
  "colon" = "Colon",
  "skin" = "Skin",
  "esophagus" = "Esophagus",
  "oesophagus" = "Esophagus",
  "kidney" = "Kidney",
  "liver" = "Liver",
  "lung" = "Lung",
  "lymph_node" = "Lymph node",
  "other" = "Other",
  "ovary" = "Ovary",
  "pancreas" = "Pancreas",
  "prostate" = "Prostate",
  "soft_tissue" = "Soft tissue"
)


clinical_features <- bc@meta.data %>%
  filter(new_cell_id %in% colnames(sketched_mat)) %>%
  mutate(
    summarised_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "Other"
    ),
    adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
    is_blood = ifelse(
      refined_tumor_type %in% c("ALL", "CLL", "LAML", "MM"),
      "Liquid",
      "Solid"
    ),
    treated = ifelse(
      treated == "t",
      "Treated",
      ifelse(treated == "f", "Untreated", treated)
    ),
    sex = ifelse(sex == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
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
    tme_archetype,
    refined_tumor_type
  ) %>%
  as.data.frame()

cells_annot_df$summarised_tumor_site <- translat_human_sites[cells_annot_df$summarised_tumor_site]

rownames(cells_annot_df) <- cells_annot_df$new_cell_id
cells_annot_df$new_cell_id <- NULL

colnames(cells_annot_df) <- c(
  "Sex",
  "Age group",
  "Solid/Liquid",
  "Sample type",
  "Sample site",
  "Treatment",
  "TME archetype",
  "Cancer type"
)

pals <- list(
  "Sex" = sex_colors,
  "Age group" = age_colors,
  "Solid/Liquid" = sl_colors,
  "Sample type" = pm_colors,
  "Sample site" = tumor_sites_colors,
  "Treatment" = treatment_colors,
  "TME archetype" = tme_colors,
  "Cancer type" = tumor_type_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = cells_annot_df,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_legend = c("Sample site" = FALSE)
)



## get drug names
MoAs <- drugs_bc_final %>%
  filter(preferred.drug.names %in% common_drugs) %>%
  select(preferred.drug.names, collapsed.MoAs) %>%
  as.data.frame() %>%
  column_to_rownames(var = "preferred.drug.names")
colnames(MoAs) <- "Mechanism of action"
MoAs <- MoAs[rownames(sketched_mat_commondrugs), , drop = FALSE]
moa_pals <- list("Mechanism of action" = MoAs_colors)

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = MoAs,
  which = "row",
  col = moa_pals,
  show_annotation_name = FALSE
)

png(
  file = "single_cell/sctherapy_int_bc/heatmap_beyondcell_with_tc1.png",
  res = 500,
  width = 20,
  height = 18,
  units = "in"
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  ncol = 2,
  # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Sample site"
)


heat <- ComplexHeatmap::Heatmap(
  mat = sketched_mat_unique[common_drugs, ],
  # mat = t(sketched_mat),
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  #row_split = 5,
  row_title = NULL,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  show_column_dend = TRUE,
  #column_split = cells_annot_df$`Cancer type`,
  #column_split = bc_c1@meta.data[colnames(sketched_mat), "scevan_subclone"],
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(title = "BCS score"),
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(12, "in")
)

ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
  "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
  "legend_gap" = unit(1, "cm")
)
draw(
  heat,
  annotation_legend_side = "top",
  annotation_legend_list = list(tumor_site_legend)
)
dev.off()
