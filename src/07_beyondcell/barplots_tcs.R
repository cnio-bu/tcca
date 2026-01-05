library(tidyverse)
library(BPCells)
library(Seurat)
library(clustree)
library(ComplexHeatmap)

## Set options
options(future.globals.maxSize = 20 * 1024^3)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/")

# Load color palette
source(file = "/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

# Load Seurat object with therapeutic cluster annotations
bc <- readRDS("results/beyondcell_pancancer_final_res.Rds")

seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
metadata <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_annotation_raw.tsv", 
sep = "\t", header = TRUE)
seu@meta.data <- cbind(seu@meta.data, metadata[, c(29:33)])

seu <- subset(seu, subset = malignancy == TRUE)
colnames(seu) <- paste0("c", c(1:ncol(seu)))

bc$sex <- seu@meta.data[bc$new_cell_id, "sex"]
bc$tme_archetype <- seu@meta.data[bc$new_cell_id, "tme_archetype"]

## human readable origins
translat_human_sites <- c(
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "adrenal_gland" = "Adrenal gland",
    "breast" = "Breast",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue",
    "bladder" = "Bladder",
    "colon" = "Colon",
    "kidney" = "Kidney"
)

clinical_features <- bc@meta.data %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated, "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        therapeutic_clusters_k300_res0.5 = factor(therapeutic_clusters_k.300.res.0.5, levels = sort(unique(therapeutic_clusters_k.300.res.0.5)))
    )
clinical_features$summarised_tumor_site <-  translat_human_sites[clinical_features$summarised_tumor_site]

# Sex barplot
sex_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = sex)) +
geom_bar(position = "fill") +
scale_fill_manual(values = sex_colors) + 
theme_bw()

ggsave("results/figures/sex_barplot_tc.png", sex_barplot, width = 8, height = 8, dpi = 500)

# Age group tumor barplot
age_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = adult_pediatric)) +
geom_bar(position = "fill") +
scale_fill_manual(values = age_colors) + 
theme_bw()

ggsave("results/figures/age_barplot_tc.png", age_barplot, width = 8, height = 8, dpi = 500)

# Liquid/solid tumor barplot
sl_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = is_blood)) +
geom_bar(position = "fill") +
scale_fill_manual(values = sl_colors) + 
theme_bw()

ggsave("results/figures/sl_barplot_tc.png", sl_barplot, width = 8, height = 8, dpi = 500)

# Sample type barplot
pm_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = sample_type)) +
geom_bar(position = "fill") +
scale_fill_manual(values = pm_colors) + 
theme_bw()

ggsave("results/figures/pm_barplot_tc.png", pm_barplot, width = 8, height = 8, dpi = 500)

# Tumor site barplot
site_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = summarised_tumor_site)) +
geom_bar(position = "fill") +
scale_fill_manual(values = tumor_sites_colors) + 
theme_bw()

ggsave("results/figures/site_barplot_tc.png", site_barplot, width = 8, height = 8, dpi = 500)

# Treated condition barplot
treated_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = treated)) +
geom_bar(position = "fill") +
scale_fill_manual(values = treatment_colors) + 
theme_bw()

ggsave("results/figures/treatment_barplot_tc.png", treated_barplot, width = 8, height = 8, dpi = 500)

# Treated condition barplot
clinical_features$ccl_patient <- ifelse(is.na(clinical_features$patient) | clinical_features$patient != "ccl", "Patient", "Cell line")
ccl_p_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = ccl_patient)) +
geom_bar(position = "fill") +
scale_fill_manual(values = ccl_p_colors) + 
theme_bw()

ggsave("results/figures/ccl_p_barplot_barplot_tc.png", ccl_p_barplot, width = 8, height = 8, dpi = 500)

# Cancer type
cancer_type <- list(
  "Brain Cancer" = c("GBM", "MB", "OGD"),
  "Neuroblastic Tumors" = c("GNB", "NB"),
  "Blood Cancer" = c("ALL", "LAML", "CLL", "MM"),
  "Skin Cancer" = c("BCC", "SKCM", "SKSC", "SKAM", "UVM"),
  "Sarcoma/Soft Tissue Cancer" = c("SARC", "GIST", "MESO"),
  "Breast Cancer" = c("BRCA"),
  "Lung Cancer" = c("SCLC", "NSCLC", "LUAD", "LUSC", "LCLC", "PLEU"),
  "Ovarian Cancer" = c("OV"),
  "Colon/Colorectal Cancer" = c("COAD", "READ"),
  "Endometrial/Uterine Cancer" = c("CESC", "UCEC", "UCS"),
  "Liver/Biliary Cancer" = c("LIHC", "CHOL"),
  "Bladder Cancer" = c("BLCA"),
  "Head and Neck Cancer" = c("HNSC"),
  "Prostate Cancer" = c("PRAD"),
  "Kidney Cancer" = c("KRCC", "KTCC", "KIRC", "KIRCH"),
  "Esophageal Cancer" = c("ESCA", "ESCC"),
  "Pancreatic Cancer" = c("PAAD"),
  "Thyroid Cancer" = c("THCA"),
  "Gastric Cancer" = c("STAD"),
  "Miscellaneous Cancer" = c("MISC"))

cancer_type <- enframe(cancer_type, name = "broad_cancer_type", value = "tumor_type") %>%
  unnest()

clinical_features <-  clinical_features %>% 
left_join(cancer_type, by = "tumor_type")


# Broad Cancer Type barplot
cancertype_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = broad_cancer_type)) +
geom_bar(position = "fill") +
scale_fill_manual(values = cancer_type_colors) + 
guides(fill = guide_legend(ncol = 1)) +
theme_bw()

ggsave("results/figures/cancertype_barplot_tc.png", cancertype_barplot, width = 8, height = 8, dpi = 500)

# Study barplot
study_barplot <- ggplot(clinical_features, aes(x = therapeutic_clusters_k300_res0.5, fill = study)) +
geom_bar(position = "fill") +
scale_fill_manual(values = study_colors) + 
theme_bw()

ggsave("results/figures/study_barplot_tc.png", study_barplot, width = 8, height = 8, dpi = 500)

# Plot TCs distribution among patients
for (study in unique(clinical_features$study)){
  df <- clinical_features[clinical_features$study == study, ]
  patient_barplot <- ggplot(df, aes(x = patient, fill = therapeutic_clusters_k300_res0.5)) +
  geom_bar(position = "fill") +
  facet_wrap(~ study, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = tcs_colors) + 
  guides(fill = guide_legend(title = "Therapeutic Cluster"))
  theme_bw()

  ggsave(paste0("results/figures/tcs_study/", study, "_barplot_tc.png"), patient_barplot, width = 8, height = 8, dpi = 500)
}

# Number of samples per specific cancer type
clinical <- read.table("clinical_metadata_v4_clean.tsv", header = TRUE, sep = "\t")

counts_cancer <- clinical %>%
  count(tumor_type)

barplot_counts_sample <- ggplot(counts_cancer, aes(x = reorder(tumor_type, n), y = n)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "Tumor Type", y = "Number of Samples") +
      theme_minimal()

ggsave("cohort_statistics/figures/samples_per_cancer.png", plot = barplot_counts_sample, width = 4, height = 6)


# Number of cells per spectific cancer type
counts_malignants <- metadata %>%
  filter(malignancy == TRUE) %>%
  count(tumor_type)

barplot_counts_cell <- ggplot(counts_malignants, aes(x = reorder(tumor_type, n), y = n)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "Tumor Type", y = "Number of malignant cells") +
      theme_minimal()

ggsave("cohort_statistics/figures/cells_per_cancer.png", plot = barplot_counts_cell, width = 4, height = 6)

# Plot distribution of broad cell types across cancer types
metadata <- metadata %>%
  left_join(cancer_type, by = "tumor_type")

cancer_by_sample <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  select(study_sample, broad_cancer_type) %>%
  distinct()

metadata <- metadata %>%
  mutate(
    cell_type_broad = ifelse(cell_type_broad == "", "Unknown", cell_type_broad),
    broad_cancer_type = factor(broad_cancer_type, 
                             levels = names(sort(table(cancer_by_sample$broad_cancer_type), 
                                                 decreasing = TRUE)))
  )

celltype_barplot <- ggplot(metadata, aes(x = broad_cancer_type, fill = cell_type_broad)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "Cancer types", y = "Cell fraction", fill = "Cell type") +
  ggtitle("Cell type fractions across cancer types") +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12)
    )

ggsave("/home/mgonzalezb/bc-meta/figures/celltype_cancertypes.png",
  celltype_barplot,
  width = 14, 
  height = 8, 
  dpi = 500
)

# Proportion of cells from each TC across samples
clinical_features <- clinical_features %>%
  mutate(study_sample = paste0(study, "_", sample))

sample_tc_summary <- clinical_features %>%
  group_by(study_sample, therapeutic_clusters_k300_res0.5) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(study_sample) %>%
  mutate(Proportion = Count / sum(Count))

sample_tc_mat <- sample_tc_summary %>%
   select(study_sample, therapeutic_clusters_k300_res0.5, Proportion) %>%
   tidyr::pivot_wider(names_from = study_sample, values_from = Proportion, values_fill = 0) %>%
     column_to_rownames("therapeutic_clusters_k300_res0.5")

samples_annot_df <- clinical_features %>%
  select(
    study_sample,
    sex,
    adult_pediatric,
    is_blood,
    sample_type,
    summarised_tumor_site,
    treated
  ) %>%
  distinct() %>%
  as.data.frame()

rownames(samples_annot_df) <- samples_annot_df$study_sample
samples_annot_df$study_sample <- NULL
samples_annot_df <- samples_annot_df[colnames(sample_tc_mat),]
colnames(samples_annot_df) <- c(
  "Sex",
  "Age group",
  "Solid/Liquid",
  "Sample type",
  "Sample site",
  "Treatment"
)

pals <- list(
  "Sex" = sex_colors,
  "Age group" = age_colors,
  "Solid/Liquid" = sl_colors,
  "Sample type" = pm_colors,
  "Sample site" = tumor_sites_colors,
  "Treatment" = treatment_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = samples_annot_df,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_legend = c("Sample site" = FALSE)
)


tcs_annot <- data.frame(TCs = factor(as.character(0:4), levels = as.character(0:4)))

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = tcs_annot,
  which = "row",
  col = list(TCs = tcs_colors),
  show_annotation_name = FALSE
)

png(
  file = "results/heatmap_tcs_per_sample.png",
  res = 500,
  width = 19,
  height = 14,
  units = "in"
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  ncol = 2, # Split Group legend into 2 columns
  gap = unit(10, "mm"),
  title = "Sample site"
)


heat <- ComplexHeatmap::Heatmap(
  mat = as.matrix(sample_tc_mat),
  # mat = t(sketched_mat),
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = tcs_annot$tcs,
  row_title = NULL,
  #column_order = rownames(samples_annot_df[order(samples_annot_df$`Sample site`), ]),
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  column_dend_side = "bottom", 
  clustering_distance_columns = "pearson",
  show_column_names = FALSE,
  row_labels = tcs_annot$TCs,
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(title = "Cell fractions\nacross TCs"),
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(8, "in")
)

ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
  "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top", annotation_legend_list = list(tumor_site_legend))
dev.off()
