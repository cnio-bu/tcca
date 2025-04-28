library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)

setwd("/home/lmgonzalezb/Documents/bc-meta/sctherapy/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

subclone_annot <- read.table("annotations_subclones.tsv")
subclone_annot$cluster <- as.factor(subclone_annot$cluster)
subclone_annot <- subclone_annot %>%
  mutate(across(everything(), ~replace_na(.x, "Unknown")))

run_chisq <- function(df, cluster_col, variable_col) {
  results <- data.frame(
    variable = character(),
    category = character(),
    cluster = character(),
    p_value = character(),
    FDR_Adjusted_P_value = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  all_p_values <- c()
  
  for (category in unique(df[[variable_col]])) {
    # Create a binary contingency table: Specific category vs. all others
    df$binary_var <- ifelse(df[[variable_col]] == category, 
                            category, 
                            paste0("Not_", category))
    
    # Perform Chi-square test for each cluster separately
    for (cl in unique(df[[cluster_col]])) {
      df$cluster_binary <- ifelse(df[[cluster_col]] == cl, 
                                  paste0("Cluster", cl), 
                                  "Other_Clusters")
      
      contingency_table <- table(df$cluster_binary, df$binary_var)
      
      # Only run if table is valid
      if (any(chisq.test(contingency_table)$expected < 5)) {
        test_result <- fisher.test(contingency_table)
      } else {
        test_result <- chisq.test(contingency_table)
      }
      p_value <- signif(test_result$p.value, 3)
      
      # Store results
      results <- rbind(
        results,
        data.frame(
          Variable = variable_col,
          Category = category,
          Cluster = cl,
          P_value = p_value,
          FDR_Adjusted_P_value = NA,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  results$FDR_Adjusted_P_value <- p.adjust(results$P_value, method = "fdr")
  results$significance <- ifelse(results$FDR_Adjusted_P_value <= 0.05, TRUE, FALSE)
  return(results)
}

plot_ggbarstats <- function(x, y, title, xlab, ylab, legend.title, colors, filename){
  plot <- ggbarstats(
    data = subclone_annot,
    x = !!sym(x),
    y = !!sym(y),
    title = title,
    xlab = xlab,
    ylab = ylab,
    proportion.test = FALSE,
    legend.title = legend.title,
    sample.size.label.args = list(size = 3),
  ) +
    scale_fill_manual(values = colors) +
    theme(axis.title = element_blank(), 
          axis.ticks.x = element_line(colour = "black"), 
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black", hjust = 0.5),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"))
  ggsave(paste0("figures/", filename, ".png"), plot = plot, width = 12, height = 10)
}

stats.sample_type <- run_chisq(subclone_annot, "cluster", "sample_type")
plot_ggbarstats(x = "sample_type",
                y = "cluster",
                title = "Primary/Metastasis subclones across clusters",
                xlab = "Cluster",
                ylab = "sample_type",
                legend.title = "Sample type",
                colors = pm_colors, 
                filename = "pm_plot"
                )

stats.sex <- run_chisq(subclone_annot, "cluster", "sex")
plot_ggbarstats(x = "sex",
                y = "cluster",
                title = "Sex of subclones across clusters",
                xlab = "Cluster",
                ylab = "sex",
                legend.title = "Sex",
                colors = sex_colors, 
                filename = "sex_plot"
)

stats.age <- run_chisq(subclone_annot, "cluster", "adult_pediatric")
plot_ggbarstats(x = "adult_pediatric",
                y = "cluster",
                title = "Adult/pediatric subclones across clusters",
                xlab = "Cluster",
                ylab = "adult_pediatric",
                legend.title = "Age",
                colors = age_colors, 
                filename = "age_plot"
)

stats.blood <- run_chisq(subclone_annot, "cluster", "is_blood")
plot_ggbarstats(x = "is_blood",
                y = "cluster",
                title = "Solid/liquid subclones across clusters",
                xlab = "Cluster",
                ylab = "is_blood",
                legend.title = "Solid/Liquid",
                colors = sl_colors, 
                filename = "solid-liquid_plot"
)

stats.treated <- run_chisq(subclone_annot, "cluster", "treated")
plot_ggbarstats(x = "treated",
                y = "cluster",
                title = "Treatment condition across clusters",
                xlab = "Cluster",
                ylab = "treated",
                legend.title = "Treatment condition",
                colors = treatment_colors, 
                filename = "treatment_plot"
)


stats.tme <- run_chisq(subclone_annot, "cluster", "tme_archetype")
plot_ggbarstats(x = "tme_archetype",
                y = "cluster",
                title = "TME archetypes across clusters",
                xlab = "Cluster",
                ylab = "TME archetype",
                legend.title = "TME archetype",
                colors = tme_colors, 
                filename = "tme_plot"
)


stats.tumor_site <- run_chisq(subclone_annot, "cluster", "summarised_tumor_site")
plot_ggbarstats(x = "summarised_tumor_site",
                y = "cluster",
                title = "Tumor site across clusters",
                xlab = "Cluster",
                ylab = "Tumor site",
                legend.title = "Tumor site",
                colors = tumor_sites_colors, 
                filename = "tumor_site_plot"
)


# Add study and cancer type to the table of subclones
subclone_annot$study_sample <- paste0(sub("\\.[0-9]+$","", rownames(subclone_annot)))
clinical <- data.table::fread("../clinical_metadata_v4_clean.tsv")
clinical$study_sample <- paste0(clinical$study, ".", clinical$sample)

clinical <- clinical %>%
  select(study_sample, study, tumor_type) %>%
  distinct()

subclone_annot <- subclone_annot %>%
  left_join(clinical, by = "study_sample")


# Add broad cancer types
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

subclone_annot <-  subclone_annot %>%
  left_join(cancer_type, by = "tumor_type")

# Cancer type barplot
sort_cancertypes <- subclone_annot %>%
  select(study_sample, broad_cancer_type) %>%
  distinct() %>%
  count(broad_cancer_type) %>%
  arrange(desc(n))

subclone_annot$broad_cancer_type <- factor(subclone_annot$broad_cancer_type,
                                           levels = sort_cancertypes$broad_cancer_type)

# Plot cancer type distribution across clusters
stats.cancer_type <- run_chisq(subclone_annot, "cluster", "broad_cancer_type")
plot_ggbarstats(x = "broad_cancer_type",
                y = "cluster",
                title = "Cancer type across clusters",
                xlab = "Cluster",
                ylab = "Cancer type",
                legend.title = "Cancer type",
                colors = broad_cancer_type_colors, 
                filename = "cancer_type_plot"
)

# Plot tumor site distribution across clusters
stats.tumor_site <- run_chisq(subclone_annot, "cluster", "study")
plot_ggbarstats(x = "study",
                y = "cluster",
                title = "Studies across clusters",
                xlab = "Cluster",
                ylab = "Study",
                legend.title = "Study",
                colors = study_colors, 
                filename = "study_plot"
)


# Compute the therapeutic heterogeneity of each sample (number of clusters per
# sample)
sample_annot <- subclone_annot %>%
  group_by(study_sample) %>%
  mutate(n_clusters = n_distinct(cluster),
         n_subclones = n()) %>%
  ungroup() %>%
  select(-cluster) %>%
  distinct() %>%
  as.data.frame()

# Plot the therapeutic heterogeneity regarding different clinical variables
plot <- ggplot(sample_annot, aes(x = factor(n_clusters), fill = broad_cancer_type)) +
  geom_bar(stat = "count") +  # "fill" makes it a proportional stacked bar+
  # geom_text(stat = "count", aes(label = ..count..), 
  #           position = position_stack(vjust = 0.5), size = 5)+
  labs(x = "Number of Clusters", y = "Proportion of Samples", 
       title = "Distribution of Clusters Across Cancer Types") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line()) +
  scale_fill_manual(values = broad_cancer_type_colors) 

ggsave("barplot_heterogeneity.png", plot = plot, width = 10, height = 10)


# Paired barplot with all variables

data <- data.frame(
  SampleID = paste0("Sample", 1:10),
  Metastasis = sample(c("Primary", "Metastasis"), 10, replace = TRUE),
  Sex = sample(c("Male", "Female"), 10, replace = TRUE),
  Treatment = sample(c("Treated", "Untreated"), 10, replace = TRUE),
  TumorSite = sample(c("Lung", "Liver", "Brain", "Bone"), 10, replace = TRUE)
)



# To keep the same order across all plots
data$SampleID <- factor(data$SampleID, levels = data$SampleID)

# Reshape to long format
long_data <- data %>%
  pivot_longer(cols = c(Metastasis, Sex, Treatment, TumorSite),
               names_to = "Category",
               values_to = "Annotation")

ggplot(long_data, aes(x = SampleID, fill = Annotation)) +
  geom_bar(position = "stack", stat = "count") +
  facet_wrap(~ Category, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Samples", y = "Annotation Presence", fill = "Annotation")

