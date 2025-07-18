library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsankey)

setwd("/home/lmgonzalezb/Documents/bc-meta/")

# Load cohort metadata
tcca.metadata <- read.table("cohort_statistics/tcca_metadata.tsv", header = TRUE)
sex <- read.table("cohort_statistics/tcca_metadata_sex_inferred.tsv", header = TRUE)
tcca.metadata$sex <- sex$sex

cohort_features <- tcca.metadata %>%
  mutate(study_sample = paste(study, sample, sep = "_")) %>%
  mutate(sample_origin = ifelse(patient != "ccl", "Patient", "Cell line")) %>%
  select(
    study_sample,
    sample,
    study,
    sample_origin,
    refined_tumor_type,
    sample_type,
    refined_tumor_site,
    treated,
    age,
    sex,
    tme_archetype
  ) %>%
  distinct()

# Reset tumor site to human readable sample origins
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

# Add cancer type
cancer_type <- list(
  "Brain Cancer" = c("GBM", "LGG"),
  "Neuroblastic Tumors" = c("NB"),
  "Blood Cancer" = c("ALL", "LAML", "CLL", "MM"),
  "Skin Cancer" = c("BCC", "SKCM", "SKSC", "UVM"),
  "Sarcoma/Soft Tissue Cancer" = c("SARC", "MESO"),
  "Breast Cancer" = c("BRCA"),
  "Lung Cancer" = c("SCLC", "LUAD", "LUSC", "LCLC"),
  "Ovarian Cancer" = c("OV"),
  "Colon/Colorectal Cancer" = c("COAD", "READ"),
  "Endometrial/Uterine Cancer" = c("CESC", "UCEC", "UCS"),
  "Liver/Biliary Cancer" = c("LIHC", "CHOL"),
  "Bladder Cancer" = c("BLCA"),
  "Head and Neck Cancer" = c("HNSC"),
  "Prostate Cancer" = c("PRAD"),
  "Kidney Cancer" = c("KIRC"),
  "Esophageal Cancer" = c("ESCA"),
  "Pancreatic Cancer" = c("PAAD"),
  "Gastric Cancer" = c("STAD"),
  "Thyroid Cancer" = c("THCA")
)


cancer_type <- enframe(cancer_type, name = "broad_cancer_type", value = "refined_tumor_type") %>%
  unnest()

cohort_features <- cohort_features %>%
  left_join(cancer_type, by = "refined_tumor_type") %>%
  column_to_rownames(var = "study_sample") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(
    summarised_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "other"
    ),
    adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
    is_blood = ifelse(refined_tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
    treated = ifelse(treated == "t", "Treated", "Untreated"),
    sex = ifelse(sex == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
  ) %>%
  select(
    broad_cancer_type,
    sample_origin,
    sample_type,
    summarised_tumor_site,
    treated,
    adult_pediatric,
    sex
  ) %>%
  # Rename tumor sites and replace NA values by "Unknown"
  mutate(summarised_tumor_site = translat_human_sites[summarised_tumor_site]) %>%
  mutate(across(
    c(sample_type, treated, adult_pediatric, sex),
    ~ replace_na(.x, "Unknown")
  ))


# Set factor levels for each category in order based on sample number
cancer_type <- names(sort(table(cohort_features$broad_cancer_type)))
sample_origin_lvls <- names(sort(table(cohort_features$sample_origin)))
sample_type_lvls <- c("Unknown1", "Metastasis", "Primary")
summarised_tumor_site_lvls <- names(sort(table(cohort_features$summarised_tumor_site)))
treated_lvls <- c("Unknown2", "Treated", "Untreated")
adult_pediatric_lvls <- c("Unknown3", "Pediatric", "Adult")
sex_lvls <- c("Unknown4", "Female", "Male")

# Define the replacement for each column
unknown_mapping <- list(
  sample_type = "Unknown1",
  treated = "Unknown2",
  adult_pediatric = "Unknown3",
  sex = "Unknown4"
)

# Apply the mapping to update Unknown values in each relevant column
cohort_features <- cohort_features %>%
  mutate(
    broad_cancer_type = factor(broad_cancer_type, levels = cancer_type),
    sample_origin = factor(sample_origin, levels = sample_origin_lvls),
    sample_type = factor(
      ifelse(
        sample_type == "Unknown",
        unknown_mapping$sample_type,
        sample_type
      ),
      levels = sample_type_lvls
    ),
    summarised_tumor_site = factor(summarised_tumor_site, levels = summarised_tumor_site_lvls),
    treated = factor(
      ifelse(treated == "Unknown", unknown_mapping$treated, treated),
      levels = treated_lvls
    ),
    adult_pediatric = factor(
      ifelse(
        adult_pediatric == "Unknown",
        unknown_mapping$adult_pediatric,
        adult_pediatric
      ),
      levels = adult_pediatric_lvls
    ),
    sex = factor(ifelse(sex == "Unknown", unknown_mapping$sex, sex), levels = sex_lvls)
  )


# Reorder levels
ordered_levels <- c(
  c("Other", setdiff(
    levels(cohort_features$broad_cancer_type), "Other"
  )),
  c("Other", setdiff(
    levels(cohort_features$summarised_tumor_site), "Other"
  )),
  levels(cohort_features$sample_origin),
  levels(cohort_features$sample_type),
  levels(cohort_features$treated),
  levels(cohort_features$adult_pediatric),
  levels(cohort_features$sex)
)

# Move "Unknown" label to the below part of the plot
ordered_node <- unique(ordered_levels)
ordered_next.node <- c(setdiff(ordered_node, cohort_features$broad_cancer_type))

# Transform the sample features dataframe into long format
sankey_data <- cohort_features %>%
  make_long(
    broad_cancer_type,
    summarised_tumor_site,
    sample_origin,
    sample_type,
    treated,
    adult_pediatric,
    sex
  )


# Add colors for each category
source("bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

sankey_data <- sankey_data %>%
  mutate(
    node_color = case_when(
      x == "broad_cancer_type" ~ broad_cancer_type_colors[node],
      x == "sample_origin" ~ patient_ccl_colors[node],
      x == "sample_type" ~ pm_colors[node],
      x == "summarised_tumor_site" ~ tumor_sites_colors[node],
      x == "treated" ~ treatment_colors[node],
      x == "adult_pediatric" ~ age_colors[node],
      x == "sex" ~ sex_colors[node]
    )
  ) %>%
  mutate(
    node_color = ifelse(grepl("^Unknown", node), "#BBB9B7", node_color),
    node = factor(node, levels = ordered_node),
    next_node = factor(next_node, levels = ordered_next.node)
  )


# Create  labels column
sankey_data <- sankey_data %>%
  group_by(node) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(
    node_labels = factor(
      paste0(node, " (n = ", count, ")"),
      levels = paste0(levels(node), " (n = ", as.numeric(table(node)), ")")
    ),
    node_labels = factor(
      gsub("Unknown[0-9]+", "Unknown", node_labels),
      levels = gsub("Unknown[0-9]+", "Unknown", levels(node_labels))
    ),
    node = factor(gsub("Unknown[0-9]+", "Unknown", node), levels = unique(
      gsub("Unknown[0-9]+", "Unknown", ordered_levels)
    )),
    next_node = factor(gsub("Unknown[0-9]+", "Unknown", next_node), levels = unique(
      gsub("Unknown[0-9]+", "Unknown", levels(next_node))
    ))
  )

## Get sankey plot
colors <- sankey_data$node_color
names(colors) <- sankey_data$node_labels
colors <- colors[unique(names(colors))]

# Rename categories in the x axis of the plot
x_categories <- c(
  "broad_cancer_type" = "Primary cancer type",
  "summarised_tumor_site" = "Sample site",
  "sample_origin" = "Sample source",
  "sample_type" = "Sample type",
  "treated" = "Treatment",
  "adult_pediatric" = "Age",
  "sex" = "Sex"
)

sankey <- ggplot(
    sankey_data,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = node_labels
    )
  ) +
    geom_sankey(flow.alpha = .8, node.color = 'white') +
    scale_fill_manual(values = colors) +
    geom_sankey_label(
      data = sankey_data[sankey_data$x != "broad_cancer_type", ],
      aes(label = node_labels),
      size = 2.5,
      color = "black",
      fill = "white"
    ) +
    scale_x_discrete(breaks = names(x_categories), labels = x_categories) +
    theme_sankey(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color = "black", face = "bold")
    )


ggsave(
  sankey,
  file = "sankey_cohort_new.pdf",
  dpi = 500,
  width = 12,
  height = 8
)

