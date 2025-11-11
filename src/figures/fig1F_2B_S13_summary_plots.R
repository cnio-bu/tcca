library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggstatsplot)
library(ComplexHeatmap)

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

## Statistics using all cells in the cohort (TCCA)
# Pie chart with proportion of malignant cells versus TME cells
count <- as.data.frame(table(metadata$malignancy))
colnames(count) <- c("Malignancy", "Count")

pie <- ggplot(count, aes(x = "", y = Count, fill = Malignancy)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y")   # Transforma a gráfico circular
  scale_fill_manual(values = c("True" = "#ff9382", "False" = "#405169")) +
  theme_void() +
  theme(legend.position = "none")

ggsave("figures/pie_malignants.pdf", plot = pie, width = 8, height = 8)


# TME archetypes per tumor type
samples_with_tme <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  group_by(study_sample) %>%
  summarise(has_true = any(malignancy == "True"),
            has_false = any(malignancy == "False")) %>%
  filter(has_true & has_false) %>%
  pull(study_sample)

metadata_with_tme <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  filter(study_sample %in% samples_with_tme & tme_archetype != "none")

theme_barplot <- theme(
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 16, margin = margin(t = 6), face = "bold"),
  axis.title.y = element_text(size = 16, margin = margin(r = 6), face = "bold"),
  axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 14, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 14)
  )

sample_tme <- metadata_with_tme %>%
  select(sample, study, sex, age, treated, sample_type, refined_tumor_type, refined_tumor_site, tme_archetype_group, tme_archetype) %>%
  distinct()

barplot <- ggplot(sample_tme,
                    aes(x = refined_tumor_type, fill = tme_archetype)) +
     geom_bar(position = "fill") +
     facet_wrap(~ sample_type) +
     scale_fill_manual(values = tme_colors) +
     scale_y_continuous(labels = percent_format(accuracy = 1)) +
     labs(x = "Cancer type", y = "Proportion of samples", fill = "TME archeatype") +
     theme_bw() +
     theme_barplot +
     theme(legend.position = "None",
           axis.title = element_text(face = "bold"),
           panel.background = element_blank(),
           plot.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.line = element_line(color = "black", linewidth = 0.3)) +
    guides(fill = guide_legend(ncol = 2))

ggsave("figures/barplot_tme.pdf", plot = barplot, width = 10, height = 5)


barplot_tme <- ggbarstats(
  data = sample_tme,
  x = tme_archetype,
  y = refined_tumor_type,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = tme_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cancer type", y = "Percentage of samples", fill = "TME archetype") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 4, reverse = TRUE))

ggsave("figures/barplot_tme_archetype.pdf", plot = barplot_tme, width = 16, height = 8)

# Plot TME archetypes and other clinical variables
sample_tme <- sample_tme  %>%
  arrange(desc(tme_archetype))


ha <- HeatmapAnnotation(
  `TME` = sample_tme$tme_archetype,
  `Cancer Type` = sample_tme$refined_tumor_type,
  `Sample type` = sample_tme$sample_type,
  `Tumor site` = sample_tme$refined_tumor_site,
  col = list(
    `TME` = tme_colors,
    `Cancer Type` = tumor_type_colors,
    `Sample type` = pm_colors,
    `Tumor site` = tumor_sites_colors
  )
)

ht <- Heatmap(
  matrix(nrow = 0, ncol = nrow(sample_tme)), # matriz "dummy" sin filas
  top_annotation = ha,
  show_heatmap_legend = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = sample_tme$refined_tumor_type
)

pdf(file = "figures/top_annot_tme_per_cancer_vars.pdf", height = 5, width = 20)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Plot TME groups per cancer type
barplot_tme_groups <- ggbarstats(
  data = sample_tme,
  x = tme_archetype_group,
  y = refined_tumor_type,
  type = "proportion") +
  scale_fill_manual(values = tme_group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cancer type", y = "Percentage of samples", fill = "TME archetype") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 2, reverse = TRUE))

ggsave("figures/barplot_tme_groups.pdf", plot = barplot_tme_groups, width = 14, height = 9)

# Plot TME groups per sex, age and primary/met samples
barplot_tme_groups <- ggbarstats(
  data = sample_tme,
  x = tme_archetype_group,
  y = treated,
  type = "proportion") +
  scale_fill_manual(values = tme_group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Treatment condition", y = "Percentage of samples", fill = "TME group") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 4, reverse = TRUE))

ggsave("figures/barplot_tme_archetype_treated.pdf", plot = barplot_tme_groups, width = 4, height = 8)

# Plot distribution of TME archetypes and general groups across scTherapy clusters
subclone_metadata <- metadata %>%
  select(
    scTherapy_cluster,
    sample,
    study,
    sample_type,
    treated,
    age_group,
    refined_tumor_type,
    refined_tumor_site,
    tme_archetype_group,
    tme_archetype,
    scevan_subclone
  ) %>%
  filter(scTherapy_cluster != "" & tme_archetype != "none") %>%
  distinct()

barplot_tme_groups_tc <- ggbarstats(
  data = subclone_metadata,
  x = scTherapy_cluster,
  y = tme_archetype_group,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = sctherapy_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "TME group", y = "Percentage of samples", fill = "Therapeutic clusters") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 10, reverse = TRUE))

ggsave("figures/barplot_tme_groups_tcs.pdf", plot = barplot_tme_groups_tc, width = 6, height = 8)


barplot_tme_tc <- ggbarstats(
  data = subclone_metadata,
  x = scTherapy_cluster,
  y = tme_archetype,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = sctherapy_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "TME archetype", y = "Percentage of samples", fill = "Therapeutic clusters") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 3, reverse = TRUE))

ggsave("figures/barplot_tme_tcs.pdf", plot = barplot_tme_tc, width = 8, height = 10)

# Compare TME archetypes of samples distributed across different number of TCs
sample_tc_counts <- subclone_metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  group_by(study_sample) %>%
  summarise(n_TC = n_distinct(scTherapy_cluster),
            refined_tumor_type = first(refined_tumor_type),
            refined_tumor_site = first(refined_tumor_site),
            tme_archetype_group = first(tme_archetype_group),
            tme_archetype = first (tme_archetype)) %>%
  mutate(n_TC = factor(n_TC, levels = as.character(1:5)))


barplot_tme_ntc <- ggbarstats(
  data = sample_tc_counts,
  x = tme_archetype_group,
  y = n_TC,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = tme_group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Number of TCs", y = "Percentage of samples", fill = "TME archetype group") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave("figures/barplot_tme_groups_ntc.pdf", plot = barplot_tme_ntc, width = 8, height = 11)

# For samples with only one TC plot the distribution of TMEs across TCs
sample_one_tc <- sample_tc_counts %>%
  filter(n_TC == 1) %>%
  pull(study_sample)
sample_one_tc_metadata <- subclone_metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  filter(study_sample %in% sample_one_tc) %>%
  select(scTherapy_cluster,
         study_sample,
         sample_type,
         treated,
         age_group,
         refined_tumor_type,
         refined_tumor_site,
         tme_archetype_group,
         tme_archetype) %>%
  distinct()

barplot_tme_onetc <- ggbarstats(
  data = sample_one_tc_metadata,
  x = tme_archetype_group,
  y = scTherapy_cluster,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = tme_group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "scTherapy cluster", y = "Percentage of samples", fill = "TME archetype group") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave("figures/barplot_tme_groups_onetc.pdf", plot = barplot_tme_onetc, width = 15, height = 8)


ggbarstats(
  data = sample_one_tc_metadata,
  x = tme_archetype_group,
  y = scTherapy_cluster,
  type = "proportion"
) 
  scale_fill_manual(values = tme_group_colors) 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) 
  facet_wrap2(~ refined_tumor_type) 
  labs(x = "scTherapy cluster", y = "Percentage of samples",
       fill = "TME archetype group")

top_anno <- HeatmapAnnotation(
  TME = sample_one_tc_metadata$tme_archetype_group,
  Cluster = sample_one_tc_metadata$scTherapy_cluster,
  `Cancer type` = sample_one_tc_metadata$refined_tumor_type,
  `Sample type` = sample_one_tc_metadata$sample_type,
  `Treatment status` = sample_one_tc_metadata$treated,
  `Age group` = sample_one_tc_metadata$age_group,
  col = list(
    TME = tme_group_colors,
    Cluster = sctherapy_colors,
    `Cancer type` = tumor_type_colors,
    `Sample type` = pm_colors,
    `Treatment status` = treatment_colors,
    `Age group` = age_group_colors
  ),
  annotation_legend_param = list(
    TME = list(title = "TME archetype group"),
    Cluster = list(title = "Therapeutic cluster"),
    `Cancer type` = list(title = "Cancer type"),
    `Sample type` = list(title = "Sample type"),
    `Treatment status` = list(title = "Treatment status"),
    `Age group` = list(title = "Age group")
  )
)

pdf(
  file = "figures/nheatmap_test.pdf",
  width = 12,
  height = 10,
)

Heatmap(
  matrix(1, nrow = 1, ncol = nrow(sample_one_tc_metadata)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_anno,
  column_split = sample_one_tc_metadata$scTherapy_cluster,
  column_order = order(sample_one_tc_metadata$tme_archetype_group)
)

dev.off()

# Compare TME archetypes of samples with different number of subclones
sample_subclone_counts <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  select(study_sample, tme_archetype_group, tme_archetyp, scevan_subclone) %>%
  group_by(study_sample) %>%
  summarise(n_subclones = n_distinct(scevan_subclone),
            refined_tumor_type = first(refined_tumor_type),
            refined_tumor_site = first(refined_tumor_site),
            tme_archetype_group = first(tme_archetype_group),
            tme_archetype = first (tme_archetype)) %>%
  mutate(n_TC = factor(n_TC, levels = as.character(1:5)))

barplot_tme_ntc <- ggbarstats(
  data = sample_tc_counts,
  x = tme_archetype_group,
  y = n_TC,
  type = "proportion"          # show proportions instead of counts
) 
  scale_fill_manual(values = tme_group_colors) 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) 
  labs(x = "Number of TCs", y = "Percentage of samples", fill = "TME archetype group") 
  theme_bw() 
  theme_barplot 
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) 
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave("figures/barplot_tme_groups_ntc.png", plot = barplot_tme_ntc, width = 15, height = 12)


## Plot top drugs per cluster ##
subclones <- read.table("drug_response_subclone_sctherapy_pandrugs.tsv", header = TRUE)

# Get number of subclones per cluster
cluster_sizes <- subclones %>%
  select(Subclone.Name, scTherapy.Cluster) %>%
  distinct() %>%
  count(scTherapy.Cluster, name = "n_subclones")

# Count how often each drug appears in each cluster
drug_counts <- subclones %>%
  group_by(Drug.Name, scTherapy.Cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Frequency of each drug per subclone
drug_freq <- drug_counts %>%
  left_join(cluster_sizes, by = "scTherapy.Cluster") %>%
  mutate(freq = count / n_subclones)
         
# Create a Drug x Cluster matrix
drug_cluster_matrix <- drug_freq %>%
  select(Drug.Name, scTherapy.Cluster, freq) %>%
  pivot_wider(names_from = scTherapy.Cluster, values_from = freq, values_fill = 0)

drug_cluster_freq_mat <- drug_cluster_matrix %>%
  column_to_rownames("Drug.Name") %>%
  as.matrix()

# Order columns from 1 to 10
drug_cluster_freq_mat <- drug_cluster_freq_mat[, as.character(1:10)]


# # Compute top drugs for clusters
# top_drug_list <- list()
# for (cluster in c(1:10)){
#   specificity_score <- drug_cluster_freq_mat[, cluster] / 
#     (rowSums(drug_cluster_freq_mat)  eps)
#   
#   top_specific_drugs <- names(sort(specificity_score, decreasing = TRUE))[1:10]
#   top_drug_list[[paste0("Cluster_", cluster)]] <- top_specific_drugs
# }
# 
top_drug_list <- list()

# Build list of top drugs for each cluster
for (cluster in 1:10) {
  cluster_drugs <- drug_cluster_freq_mat[, cluster]
  top_drugs <- cluster_drugs[cluster_drugs > 0.6]
  top_drug_list[[paste0("Cluster_", cluster)]] <- top_drugs
}

# Get maximum length across all clusters
max_len <- max(sapply(top_drug_list, length))

# Pad each list element with NAs to make them equal length
top_drug_list_padded <- lapply(top_drug_list, function(x) {
  drugs <- names(x)
  c(drugs, rep(NA, max_len - length(x)))
})

# Convert to data frame
top_drugs_df <- as.data.frame(top_drug_list_padded, stringsAsFactors = FALSE)


write.table(top_drugs_df, "top_drugs_cluster.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)



# Plot MoAs of top 10 drugs
drug_moas <- subclones %>%
  select(Drug.Name, Drug.Mechanism.of.Action) %>%
  distinct()
top_drugs_moas <- top_drugs_df %>%
  pivot_longer(cols = everything(), names_to = "Cluster", values_to = "Drug") %>%
  left_join(select(drug_moas, Drug.Mechanism.of.Action, Drug.Name), 
            by = c("Drug" = "Drug.Name")) %>%
  mutate(Cluster = as.numeric(gsub("Cluster_", "", Cluster))) %>%
  arrange(Cluster, Drug.Mechanism.of.Action) %>% 
  mutate(Cluster = factor(Cluster, levels = 1:10),
         Drug = factor(Drug, levels = unique(.$Drug)),
         Drug.Mechanism.of.Action = factor(Drug.Mechanism.of.Action, 
                                           levels = c(setdiff(sort(Drug.Mechanism.of.Action), 
                                                              "Other"), 
                                                      "Other")))

# Remove NAs
top_drugs_moas <- top_drugs_moas %>%
  filter(!is.na(Drug) & !is.na(Drug.Mechanism.of.Action))

tile_plot <- ggplot(top_drugs_moas, aes(x = Cluster, y = Drug, fill = Drug.Mechanism.of.Action)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = MoAs_colors) +
  labs(fill = "MoA") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12, margin = margin(t = 6), face = "bold"),
    axis.title.y = element_text(size = 12, margin = margin(r = 6), face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2))

ggsave("figures/top_drug_tileplot.pdf", plot = tile_plot, width = 5, height = 12, dpi = 500)


## MoAs per TME archetype
drugs_intersection <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  group_by(study_sample, Subclone.Name) %>%
  summarise(drugs = list(Drug.Name), .groups = "drop") %>%
  group_by(study_sample) %>%
  summarise(intersected_drugs = Reduce(intersect, drugs), .groups = "drop")

moas <- subclones %>%
  select(Drug.Name, Drug.Mechanism.Of.Action
  ) %>%
  distinct()

tme_sample <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  select(study_sample, TME.Archetype.Group) %>%
  distinct()

# Join the moas, tme and common drugs per sample
drugs_intersection <- drugs_intersection %>%
  left_join(moas, by = c("intersected_drugs" = "Drug.Name")) %>%
  left_join(tme_sample, by = "study_sample")


# Plot the MoAs of drugs shared across subclones from the same sample, grouped 
# by TME archetype
drugs_intersection <- drugs_intersection %>%
  mutate(TME.Archetype = ifelse(TME.Archetype.Group == "none", "None", TME.Archetype.Group))

drugs_intersection$TME.Archetype.Group <- factor(drugs_intersection$TME.Archetype.Group, 
                                           levels = c(setdiff(unique(drugs_intersection$TME.Archetype.Group), 
                                                              "None"), "None"))
barplot <- ggplot(drugs_intersection,
                  aes(x = TME.Archetype.Group, fill = Drug.Mechanism.Of.Action)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = MoAs_colors) +
  labs(x = "Cancer type", y = "Proportion of samples", fill = "TME archeatype") +
  theme_bw() +
  theme_barplot +
  theme(legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        plot.margin = unit(c(1, 1, 3, 4), "cm")) +
  guides(fill = guide_legend(ncol = 3))

ggsave("figures/moas_tme.png", 
       plot = barplot, 
       width = 15, 
       height = 18, 
       dpi = 500)



# Redo the plot grouping MoAs into broader categories
drug_categories <- c(
  "BCR-ABL inhibitor" = "Targeted Kinase Inhibitors",
  "EGFR inhibitor" = "Targeted Kinase Inhibitors",
  "MAPK inhibitor" = "Targeted Kinase Inhibitors",
  "PI3K/AKT/mTOR signaling inhibitor" = "Targeted Kinase Inhibitors",
  "JAK-STAT signaling inhibitor" = "Targeted Kinase Inhibitors",
  "VEGFR inhibitor" = "Targeted Kinase Inhibitors",
  "BRAF inhibitor" = "Targeted Kinase Inhibitors",
  "SRC inhibitor" = "Targeted Kinase Inhibitors",
  "MET inhibitor" = "Targeted Kinase Inhibitors",
  "Kinase inhibitor" = "Targeted Kinase Inhibitors",
  "Cell cycle arrest;PI3K/AKT/mTOR signaling inhibitor" = "Targeted Kinase Inhibitors",
  
  "Chromatin agent" = "Epigenetic/Transcription Modulators",
  "Transcription inhibitor" = "Epigenetic/Transcription Modulators",
  
  "Microtubule agent" = "Cytotoxic Agents",
  "DNA related agent" = "Cytotoxic Agents",
  "Cell cycle arrest" = "Cytotoxic Agents",
  
  "Ubiquitin-proteasome system inhibitor" = "Proteostasis/Stress Modulators",
  "HSP inhibitor" = "Proteostasis/Stress Modulators",
  
  "Pro-apoptotic agent" = "Apoptosis/Redox Modulators",
  "ROS/RNS modulator" = "Apoptosis/Redox Modulators",
  "p53 activator/MDM2 inhibitor" = "Apoptosis/Redox Modulators",
  
  "NAMPT inhibitor" = "Metabolic Inhibitors",
  
  "Other" = "Other/Unclassified"
)

drugs_intersection <- drugs_intersection %>%
  mutate(MoAs_broad = recode(Drug.Mechanism.of.Action, !!!drug_categories))

barplot <- ggplot(drugs_intersection,
                  aes(x = TME.Archetype, fill = MoAs_broad)) 
  geom_bar(position = "fill") 
  #scale_fill_manual(values = MoAs_colors) 
  labs(x = "Cancer type", y = "Proportion of samples", fill = "TME archeatype") 
  theme_bw() 
  theme_barplot 
  theme(legend.position = "bottom",
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        plot.margin = unit(c(1, 1, 3, 4), "cm")) 
  guides(fill = guide_legend(ncol = 2))

ggsave("figures/moas_tme_broad.png", 
       plot = barplot, 
       width = 10, 
       height = 15, 
       dpi = 500)


tcs_sample <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  group_by(study_sample) %>%
  summarise(n_clusters = n_distinct(scTherapy.Cluster), .groups = "drop")

samples_in_three_or_more <- tcs_sample %>%
  filter(n_clusters >= 2) %>%
  pull(study_sample)

samples_with_tme <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  mutate(three_or_more = ifelse(study_sample %in% samples_in_three_or_more,
                                "YES", "NO")) %>%
  select(study_sample, TME.Archetype, three_or_more) %>%
  distinct()

barplot <- ggplot(samples_with_tme,
                  aes(x = three_or_more, fill = TME.Archetype)) 
  geom_bar(position = "fill") 
  scale_fill_manual(values = tme_colors) 
  labs(x = "Sample in three or more TCs?", 
       y = "Proportion of samples", 
       fill = "TME archeatype") 
  theme_bw() 
  theme_barplot 
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3)) 
  guides(fill = guide_legend(ncol = 1))


# For study eac_thomas_carroll
subclones_eac <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  filter(Study == "eac_thomas_carroll") %>%
  select(-c(2:6)) %>%
  distinct()

tabla_df <- subclones_eac %>%
  count(Sample, scTherapy.Cluster) %>%
  tidyr::pivot_wider(names_from = scTherapy.Cluster, values_from = n, values_fill = 0)
  
tabla_df <- tabla_df[, c("Sample", as.character(c(1, 2, 3, 7, 8, 9)))]
write.table(tabla_df, "eac_thomas_carrol_tab1.tsv", sep = "\t")


# Compute number of clusters each cancer type is distributed in 
tc_counts_cancertype <- metadata %>%
  select(refined_tumor_type, scTherapy_cluster, scevan_subclone) %>%
  filter(!is.na(scTherapy_cluster)) %>%
  distinct() %>%
  group_by(refined_tumor_type) %>%
  summarise(n_clusters = n_distinct(scTherapy_cluster)) %>%
  arrange(desc(n_clusters), refined_tumor_type) %>%
  mutate(refined_tumor_type = factor(refined_tumor_type, levels = rev(unique(refined_tumor_type))))


tc_count_plot <- ggplot(tc_counts_cancertype, aes(x = refined_tumor_type, y = n_clusters, fill = refined_tumor_type)) +
                        geom_col() +
                        scale_fill_manual(values = tumor_type_colors, name = "Cancer types") +
                        labs(x = "Cancer type", y = "Number of therapeutic clusters", title = "Number of TCs per cancer type") +
                        scale_y_continuous(breaks = seq(1, 10, by = 1)) +
                        theme_minimal() +
                        coord_flip() +
                        guides(fill = guide_legend(reverse = TRUE), ncol = 4) +
                        theme_barplot +
                        theme(axis.text.x = element_text(angle = 0),
                              legend.position = "top")

ggsave(
  "figures/tc_count_cancertype.pdf",
  plot = tc_count_plot,
  height = 8,
  width = 8,
  dpi = 300
)

# Plot TCs distribution per cancer type
subclones_clinical <- metadata %>%
  select(refined_tumor_type, scTherapy_cluster, scevan_subclone, sample_type, treated, age_group) %>%
  filter(!is.na(scTherapy_cluster)) %>%
  mutate(scTherapy_cluster = factor(scTherapy_cluster)) %>%
  distinct()

barplot_tc_cancertype <- ggbarstats(
  data = tc_cancertype ,
  x = scTherapy_cluster,
  y = refined_tumor_type,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = sctherapy_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cancer type", y = "Percentage of samples", fill = "Therapeutic cluster") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 10, reverse = TRUE))

ggsave("figures/barplot_tcs_cancertype.pdf", plot = barplot_tc_cancertype, width = 12, height = 7)


# Primary/Met subclones per cluster and cancer type
barplot_tc_clinical <- ggbarstats(
  data = subclones_clinical,
  x = scTherapy_cluster,
  y = age_group,
  type = "proportion"          # show proportions instead of counts
) +
  scale_fill_manual(values = sctherapy_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Age group", y = "Percentage of samples", fill = "Therapeutic cluster") +
  theme_bw() +
  theme_barplot +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

ggsave("./figures/barplot_tcs_age.pdf", plot = barplot_tc_clinical, width = 6, height = 8)

# Plot proportion of samples per cluster
seu@meta.data <- seu@meta.data %>%
  left_join(select(metadata_tcca, cell, refined_tumor_type), by = "cell")

barplot <- ggplot(seu@meta.data, aes(x = clusters, fill = refined_tumor_type)) 
  geom_bar(position = "fill") 
  scale_fill_manual(values = tumor_type_colors) 
  labs(x = "Cancer type", y = "Cell fraction", fill = "Tumor type") 
  ggtitle("Tumor types of cells across clusters") 
  theme_bw() 
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
ggsave("marker_genes/barplot_cell_cancertype.pdf", plot = barplot, height = 6, width = 10)


# Plot distribution of cancer types into TCs at subclone-level
barplot <- ggplot(subclones_clinical, aes(x = scTherapy_cluster, fill = refined_tumor_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tumor_type_colors) +
  labs(x = "Therapeutic clusters", y = "Subclone fraction", fill = "Cancer type") +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
ggsave("figures/barplot_subclone_cancertype.pdf", plot = barplot, height = 6, width = 10)

