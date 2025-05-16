library(tidyverse)
library(BPCells)
library(Seurat)
library(clustree)

## Set options
options(future.globals.maxSize = 20 * 1024 ^ 3)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load color palette
source(file = "/home/mgonzalezb/bc-meta/TCCA_palette.R")

# Load Seurat object with therapeutic cluster annotations
bc <- readRDS("beyondcell/beyondcell_pancancer.Rds")

seu <- readRDS(
  "./single_cell/seurat/v5/lvl2/seu_lvl2_sex_infered.rds"
)
metadata <- read.table(
  "./single_cell/seurat/tcca/tcca_annotation_raw.tsv",
  sep = "\t",
  header = TRUE
)
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
    is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
    treated = ifelse(treated, "Treated", "Untreated"),
    sex = ifelse(sex == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
    therapeutic_clusters_k300_res0.5 = factor(therapeutic_clusters_k.300.res.0.5, levels = sort(
      unique(therapeutic_clusters_k.300.res.0.5)
    ))
  )
clinical_features$summarised_tumor_site <-  translat_human_sites[clinical_features$summarised_tumor_site]

# Sex barplot
sex_barplot <- ggplot(clinical_features,
                      aes(x = therapeutic_clusters_k300_res0.5, fill = sex)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = sex_colors) +
  theme_bw()

ggsave(
  "results/figures/sex_barplot_tc.png",
  sex_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Age group tumor barplot
age_barplot <- ggplot(clinical_features,
                      aes(x = therapeutic_clusters_k300_res0.5, fill = adult_pediatric)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = age_colors) +
  theme_bw()

ggsave(
  "results/figures/age_barplot_tc.png",
  age_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Liquid/solid tumor barplot
sl_barplot <- ggplot(clinical_features,
                     aes(x = therapeutic_clusters_k300_res0.5, fill = is_blood)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = sl_colors) +
  theme_bw()

ggsave(
  "results/figures/sl_barplot_tc.png",
  sl_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Sample type barplot
pm_barplot <- ggplot(clinical_features,
                     aes(x = therapeutic_clusters_k300_res0.5, fill = sample_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = pm_colors) +
  theme_bw()

ggsave(
  "results/figures/pm_barplot_tc.png",
  pm_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Tumor site barplot
site_barplot <- ggplot(
  clinical_features,
  aes(x = therapeutic_clusters_k300_res0.5, fill = summarised_tumor_site)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tumor_sites_colors) +
  theme_bw()

ggsave(
  "results/figures/site_barplot_tc.png",
  site_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Treated condition barplot
treated_barplot <- ggplot(clinical_features,
                          aes(x = therapeutic_clusters_k300_res0.5, fill = treated)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = treatment_colors) +
  theme_bw()

ggsave(
  "results/figures/treatment_barplot_tc.png",
  treated_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Treated condition barplot
clinical_features$ccl_patient <- ifelse(
  is.na(clinical_features$patient) |
    clinical_features$patient != "ccl",
  "Patient",
  "Cell line"
)
ccl_p_barplot <- ggplot(clinical_features,
                        aes(x = therapeutic_clusters_k300_res0.5, fill = ccl_patient)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = ccl_p_colors) +
  theme_bw()

ggsave(
  "results/figures/ccl_p_barplot_barplot_tc.png",
  ccl_p_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

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
  "Miscellaneous Cancer" = c("MISC")
)

cancer_type <- enframe(cancer_type, name = "broad_cancer_type", value = "tumor_type") %>%
  unnest()

clinical_features <-  clinical_features %>%
  left_join(cancer_type, by = "tumor_type")


# Broad Cancer Type barplot
cancertype_barplot <- ggplot(
  clinical_features,
  aes(x = therapeutic_clusters_k300_res0.5, fill = broad_cancer_type)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cancer_type_colors) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw()

ggsave(
  "results/figures/cancertype_barplot_tc.png",
  cancertype_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Study barplot
study_barplot <- ggplot(clinical_features,
                        aes(x = therapeutic_clusters_k300_res0.5, fill = study)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = study_colors) +
  theme_bw()

ggsave(
  "results/figures/study_barplot_tc.png",
  study_barplot,
  width = 8,
  height = 8,
  dpi = 500
)

# Plot TCs distirbution among patients
for (study in unique(clinical_features$study)){
  df <- clinical_features[clinical_features$study == study, ]
  patient_barplot <- ggplot(df, aes(x = patient, fill = therapeutic_clusters_k300_res0.5)) +
  geom_bar(position = "fill") +
  facet_wrap(~ study, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = tcs_colors) + 
  guides(fill = guide_legend(title = "Therapeutic Cluster"))
  theme_bw()

  ggsave(paste0("results/figures/tcs_study/", study, "_barplot_tc.png"), 
  patient_barplot, width = 8, height = 8, dpi = 500)
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

