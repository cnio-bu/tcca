library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/home/lmgonzalezb/Documents/bc-meta/cohort_statistics/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)

## Statistics using all cells in the cohort (TCCA)
# Pie chart with proportion of malignant cells versus TME cells
count <- as.data.frame(table(metadata$malignancy))
colnames(count) <- c("Malignancy", "Count")

pie <- ggplot(count, aes(x = "", y = Count, fill = Malignancy)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +  # Transforma a gráfico circular
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
  filter(study_sample %in% samples_with_tme & !(refined_tumor_type %in% c("ALL", "CLL", "LAML", "MM")))

theme_barplot <- theme(
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 14, margin = margin(t = 6), face = "bold"),
  axis.title.y = element_text(size = 14, margin = margin(r = 6), face = "bold"),
  axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12)
  )

barplot <- ggplot(metadata_with_tme,
                    aes(x = refined_tumor_type, fill = tme_archetype)) +
     geom_bar(position = "fill") +
     scale_fill_manual(values = tme_colors) +
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


## Plot top drugs per cluster ##
subclones <- read.table("subclone_level_annotated.tsv", header = TRUE)

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
#     (rowSums(drug_cluster_freq_mat) + eps)
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

tile_plot <- ggplot(top_drugs_moas, aes(x = Drug, y = Cluster, fill = Drug.Mechanism.of.Action)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = MoAs_colors) +
  theme_minimal() +
  labs(title = "Mechanism of Action of Top Drugs per Cluster", fill = "MoAs") +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12, margin = margin(t = 6), face = "bold"),
    axis.title.y = element_text(size = 12, margin = margin(r = 6), face = "bold"),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("figures/top_drug_tileplot.pdf", plot = tile_plot, width = 12, height = 5, dpi = 500)


## MoAs per TME archetype
drugs_intersection <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  group_by(study_sample, Subclone.Name) %>%
  summarise(drugs = list(Drug.Name), .groups = "drop") %>%
  group_by(study_sample) %>%
  summarise(intersected_drugs = Reduce(intersect, drugs), .groups = "drop")

moas <- subclones %>%
  select(Drug.Name, Drug.Mechanism.of.Action) %>%
  distinct()

tme_sample <- subclones %>%
  mutate(study_sample = paste0(Study, "_", Sample)) %>%
  select(study_sample, TME.Archetype) %>%
  distinct()

# Join the moas, tme and common drugs per sample
drugs_intersection <- drugs_intersection %>%
  left_join(moas, by = c("intersected_drugs" = "Drug.Name")) %>%
  left_join(tme_sample, by = "study_sample")


# Plot the MoAs of drugs shared across subclones from the same sample, grouped 
# by TME archetype
drugs_intersection <- drugs_intersection %>%
  mutate(TME.Archetype = ifelse(TME.Archetype == "none", "None", TME.Archetype))

drugs_intersection$TME.Archetype <- factor(drugs_intersection$TME.Archetype, 
                                           levels = c(setdiff(unique(drugs_intersection$TME.Archetype), 
                                                              "None"), "None"))
barplot <- ggplot(drugs_intersection,
                  aes(x = TME.Archetype, fill = Drug.Mechanism.of.Action)) +
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
                  aes(x = TME.Archetype, fill = MoAs_broad)) +
  geom_bar(position = "fill") +
  #scale_fill_manual(values = MoAs_colors) +
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
                  aes(x = three_or_more, fill = TME.Archetype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tme_colors) +
  labs(x = "Sample in three or more TCs?", 
       y = "Proportion of samples", 
       fill = "TME archeatype") +
  theme_bw() +
  theme_barplot +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3)) +
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
tc_counts_cancertype <- subclones %>%
  group_by(Refined.Tumor.Type) %>%
  summarise(n_clusters = n_distinct(scTherapy.Cluster)) %>%
  arrange(desc(n_clusters), Refined.Tumor.Type) %>%
  mutate(Refined.Tumor.Type = factor(Refined.Tumor.Type, levels = rev(unique(Refined.Tumor.Type))))


tc_count_plot <- ggplot(tc_counts_cancertype, aes(x = Refined.Tumor.Type, y = n_clusters, fill = Refined.Tumor.Type)) +
                        geom_col() +
                        scale_fill_manual(values = tumor_type_colors, name = "Cancer types") +
                        labs(x = "Cancer type", y = "Number of therapeutic clusters", title = "Number of TCs per cancer type") +
                        scale_y_continuous(breaks = seq(1, 10, by = 1)) +
                        theme_minimal() +
                        coord_flip() +
                        guides(fill = guide_legend(reverse = TRUE)) +
                        theme_barplot +
                        theme(axis.text.x = element_text(angle = 0))

ggsave(
  "figures/tc_count_cancertype.png",
  plot = tc_count_plot,
  height = 8,
  width = 8,
  dpi = 300
)



# Primary/Met subclones per cluster and cancer type
barplot <- ggplot(subclones, aes(x = factor(scTherapy_cluster), fill = treated)) +
  geom_bar(position = "stack") +
  facet_wrap(~ refined_tumor_type, ncol = 10) +
  labs(
    title = "Treated/Untreated per Cluster and Tumor Type",
    x = "Cluster",
    y = "Number of subclones",
    fill = "Treatment condition"
  ) +
  scale_fill_manual(values = c("t" = "#6ED1BC", "f" = "#D18B6E")) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("./figures/treatede_cluster_tumortype.png", plot = barplot, width = 16, height = 10)

# Plot proportion of samples per cluster
seu@meta.data <- seu@meta.data %>%
  left_join(select(metadata_tcca, cell, refined_tumor_type), by = "cell")

barplot <- ggplot(seu@meta.data, aes(x = clusters, fill = refined_tumor_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tumor_type_colors) +
  labs(x = "Cancer type", y = "Cell fraction", fill = "Tumor type") +
  ggtitle("Tumor types of cells across clusters") +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
ggsave("marker_genes/barplot_cell_cancertype.pdf", plot = barplot, height = 6, width = 10)
