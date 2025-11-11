library(dplyr)
library(tidyverse)
library(ggplot2)


setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/cohort_statistics/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)
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

metadata <- metadata %>%
  mutate(
    tme_archetype = factor(
      tme_archetype,
      levels = c(
        "Immune_rich",
        "Immune_rich_Treg_cDC2_bias",
        "Tcell_centric",
        "Myeloid_centric",
        "Myeloid_centric_Mp_bias",
        "Myeloid_centric_Mo_bias(high)",
        "Myeloid_centric_Mo_bias(medium)",
        "Immune_stromal",
        "Immune_stromal_Endolike",
        "Immune_stromal_CAFlike_Mp_bias",
        "Immune_desert_CAFlike",
        "Immune_stromal_desert",
        "none"
      )
    ),
    tme_archetype_group = factor(
      tme_archetype_group,
      levels = c(
        "Immune_rich",
        "Tcell_centric",
        "Myeloid_centric",
        "Immune_stromal",
        "Immune_desert",
        "none"
      )
    ),
    sample_type = factor(
      case_when(
        sample_type == "p" ~ "Primary",
        sample_type == "m" ~ "Metastasis"
      ),
      levels = c("Primary", "Metastasis")
    ),
    sex = factor(
      case_when(
        sex == "m" ~ "Male",
        sex == "f" ~ "Female",
        sex == "" ~ "Unknown"
      ),
      levels = c("Male", "Female", "Unknown")
    ),
    treated = factor(
      case_when(
        treated == "t" ~ "Treated",
        treated == "f" ~ "Untreated",
        TRUE ~ "Unknown"
      ),
      levels = c("Untreated", "Treated", "Unknown")
    ),
    age_group = factor(
      case_when(
        age >= 0 & age <= 15 ~ "Pediatric",
        age >= 16 & age <= 39 ~ "Young adult",
        age >= 40 & age <= 64 ~ "Adult",
        age >= 65 ~ "Elderly",
        TRUE ~ "Unknown" # catch missing or invalid ages
      ),
      levels = c("Pediatric", "Young adult", "Adult", "Elderly", "Unknown")
    ),
    refined_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "Other"
    )
  )
metadata$refined_tumor_site <- translat_human_sites[metadata$refined_tumor_site]


theme_barplot <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(), 
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 14, margin = margin(t = 6), face = "bold"),
  axis.title.y = element_text(size = 14, margin = margin(r = 6), face = "bold"),
  axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
  axis.text.y = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12)
)

#### Focus on BRCA ###
subclones_metadata <- metadata %>%
  select(
    scevan_subclone, 
    study, 
    patient, 
    sample, 
    sex, 
    age, 
    age_group, 
    refined_tumor_type, 
    tumor_subtype,
    sample_type, 
    treated, 
    refined_tumor_site,
    tme_archetype,
    tme_archetype_group,
    scTherapy_cluster
    ) %>%
    distinct() %>% 
    filter(!is.na(scTherapy_cluster))

subclones_brca <- subclones_metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  filter(refined_tumor_type == "BRCA") %>%
  mutate(
    age_group_brca = factor(cut(age, breaks = c(0, 45, 55, Inf), right = FALSE,
                    labels = c("<45", "45-55", "+55")),
                    levels =  c("<45", "45-55", "+55"))
  ) %>%
  mutate(across(where(is.character), ~ if_else(. == "", "Unknown", .))) %>%
  distinct()


  
subclones_brca_long <- subclones_brca %>%
  select(sample_type, tme_archetype, tumor_subtype, age_group_brca, scTherapy_cluster) %>%
  pivot_longer(cols = c(sample_type, tme_archetype, tumor_subtype, age_group_brca),
               names_to = "Variable",
               values_to = "Category") %>%
  mutate(Category = factor(Category, levels = rev(c("<45", "45-55", "+55", 
                                                "Primary", "Metastasis",
                                                "Immune_rich",
                                                "Immune_rich_Treg_cDC2_bias",
                                                "Tcell_centric",
                                                "Myeloid_centric",
                                                "Myeloid_centric_Mp_bias",
                                                "Myeloid_centric_Mo_bias(high)",
                                                "Myeloid_centric_Mo_bias(medium)",
                                                "Immune_stromal",
                                                "Immune_stromal_Endolike",
                                                "Immune_stromal_CAFlike_Mp_bias",
                                                "Immune_desert_CAFlike",
                                                "Immune_stromal_desert",
                                                "none",
                                                "brca1_preneoplastic",
                                                "tnbc",
                                                "her2",
                                                "pr",
                                                "er",
                                                "luminal_a",
                                                "unknown"
                                              ))),
          scTherapy_cluster = factor(scTherapy_cluster))


barplot <- ggplot(subclones_brca_long, aes(x = Category,
                                           fill = scTherapy_cluster)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_grid(rows = vars(Variable), scales = "free_y", space = "free_y") +
  labs(
    x = "Number of subclones",   # al rotar, cambia de eje
    fill = "scTherapy cluster"
  ) +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(
    panel.spacing = unit(1.5, "lines"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  coord_flip()



subclones_pct <- subclones_brca_long %>%
  group_by(Variable, Category, scTherapy_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Variable, Category) %>%
  mutate(pct = count / sum(count)) %>%
  ungroup() %>%
  mutate(scTherapy_cluster = factor(scTherapy_cluster))

# Graficar barras agrupadas con proporciones
baplot <- ggplot(subclones_pct, aes(x = Category, y = pct, fill = scTherapy_cluster)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(~ Variable, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Percentage", fill = "scTherapy cluster") +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(
    panel.spacing = unit(1.5, "lines"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1)
  )
ggsave("./figures/treatede_cluster_tumortype.pdf", plot = barplot, width = 10, height = 10)


barplot <- ggplot(subclones_brca_long, aes(x = Category,
                                           fill = scTherapy_cluster)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_grid(~ Variable, scales = "free_x", space = "free_x") +
  labs(
    y = "Fraction of subclones",
    fill = "scTherapy cluster"
  ) +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1))


samples_brca <- subclones_brca %>%
  select(c(3:27)) %>%
  distinct()

barplot <- ggplot(samples_brca, aes(x = age_group,
                                           fill = tme_archetype)) +
  geom_bar(position = "fill", width = 0.7) +
  labs(
    y = "Fraction of subclones",
    fill = "TME archetype"
  ) +
  scale_fill_manual(values = tme_colors) +
  theme_barplot +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1))


## MoAs of drugs predicted for subclones within each age
subclones_brca_MoAs <- subclones %>%
  mutate(Study_Sample = paste0(Study, "_", Sample)) %>%
  filter(Refined.Tumor.Type == "BRCA") %>%
  mutate(
    age_group = factor(cut(Age, breaks = c(0, 45, 55, Inf), right = FALSE,
                           labels = c("<45", "45-55", "+55")),
                       levels =  c("<45", "45-55", "+55"))
  ) %>%
  select(1,2,4, age_group) %>%
  distinct()


barplot <- ggplot(subclones_brca_MoAs, aes(x = age_group,
                                       fill = Drug.Mechanism.of.Action)) +
  geom_bar(position = "fill", width = 0.7) +
  labs(
    y = "Fraction of subclones",
    fill = "Age group"
  ) +
  scale_fill_manual(values = MoAs_colors) +
  theme_barplot +
  guides(fill = guide_legend(ncol = 1))


