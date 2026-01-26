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



# Fisher test table (specific enrichments)
tcca_metadata <- read.table("tcca_metadata.tsv", header = T, sep = "\t")
tcca_metadata[tcca_metadata$tme_archetype == "Immune_stromal_desert", "tme_archetype_group"] <- "Immune_desert"
subclone_level <- tcca_metadata %>%
  filter(!is.na(scevan_subclone) & !is.na(scTherapy_cluster)) %>%
  select(scevan_subclone, refined_tumor_type, scTherapy_cluster) %> %>% 
  distinct()

# GBM in TC5
x <- table(subclone_level$refined_tumor_type, subclone_level$scTherapy_cluster)
gbm_tc5 <- x["GBM", "5"]
gmb_not_tc5 <- sum(x["GBM", ]) - gbm_tc5

not_gbm_tc5 <- sum(x[, "5"]) - gbm_tc5
not_gbm_not_tc5 <- sum(x) - (gbm_tc5 + gmb_not_tc5 + not_gbm_tc5)

mat <- matrix(c(gbm_tc5, gmb_not_tc5, not_gbm_tc5, not_gbm_not_tc5), nrow = 2, ncol = 2)
fisher.test(mat)

# AML in TC6
laml_tc6 <- x["LAML", "6"]
laml_not_tc6 <- sum(x["LAML", ]) - laml_tc6

not_laml_tc6 <- sum(x[, "6"]) - laml_tc6
not_laml_not_tc6 <- sum(x) - (laml_tc6 + laml_not_tc6 + not_laml_tc6)
mat <- matrix(c(laml_tc6, laml_not_tc6, not_laml_tc6, not_laml_not_tc6), nrow = 2, ncol = 2)
fisher.test(mat)


# Fisher tests for TME archetypes in metastatic vs primary samples
mixed_samples <- tcca_metadata %>%
  group_by(study, sample) %>%
  summarise(
    has_malignant = any(malignancy == "True"),
    has_non_malignant = any(malignancy == "False"),
    .groups = "drop"
  ) %>%
  filter(has_malignant & has_non_malignant) %>%
  mutate(study_sample = paste0(study, "_", sample))

tme_samples <- tcca_metadata %>%
  filter(tme_archetype != "none") %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  filter(study_sample %in% mixed_samples$study_sample) %>%
  select(study_sample, sample, refined_tumor_type, tme_archetype_group,sample_type, treated, age, sex) %>%
  distinct()

# Cancer types
mat <- table(tme_samples$refined_tumor_type, tme_samples$tme_archetype)
fisher.test(mat, simulate.p.value = TRUE, B = 10000)

# Immune-desert
mat <- table(tme_samples$tme_archetype_group, tme_samples$sample_type)
met_desert <- mat["Immune_desert", "m"]
pri_desert <- mat["Immune_desert", "p"]
met_notdesert <- sum(mat[, "m"]) - met_desert
pri_notdesert <- sum(mat[, "p"]) - pri_desert

fisher.test(matrix(c(met_desert, met_notdesert,
                     pri_desert, pri_notdesert), nrow=2, byrow=T))

# Immune-stromal
met_stromal <- mat["Immune_stromal", "m"]
pri_stromal <- mat["Immune_stromal", "p"]
met_notstromal <- sum(mat[, "m"]) - met_stromal
pri_notstromal <- sum(mat[, "p"]) - pri_stromal

fisher.test(matrix(c(met_stromal, met_notstromal,
                     pri_stromal, pri_notstromal), nrow=2, byrow=T))

# T-cell-centric
met_tcell <- mat["Tcell_centric", "m"]
pri_tcell <- mat["Tcell_centric", "p"]
met_nottcell <- sum(mat[, "m"]) - met_tcell
pri_nottcell <- sum(mat[, "p"]) - pri_tcell
fisher.test(matrix(c(met_tcell, met_nottcell,
                     pri_tcell, pri_nottcell), nrow=2, byrow=T))

# Immune-rich
met_rich <- mat["Immune_rich", "m"]
pri_rich <- mat["Immune_rich", "p"]
met_notrich <- sum(mat[, "m"]) - met_rich
pri_notrich <- sum(mat[, "p"]) - pri_rich
fisher.test(matrix(c(met_rich, met_notrich,
                     pri_rich, pri_notrich), nrow=2, byrow=T))


# Fisher tests for TME archetypes in treated vs untreated samples
mat <- table(tme_samples$tme_archetype_group, tme_samples$treated)
# Immune-stromal
untreated_stromal <- mat["Immune_stromal", "f"]
treated_stromal <- mat["Immune_stromal", "t"]
untreated_notstromal <- sum(mat[, "f"]) - untreated_stromal
treated_notstromal <- sum(mat[, "t"]) - treated_stromal

fisher.test(matrix(c(untreated_stromal, untreated_notstromal,
                     treated_stromal, treated_notstromal), nrow=2, byrow=T))
# T-cell-centric
untreated_tcell <- mat["Tcell_centric", "f"]
treated_tcell <- mat["Tcell_centric", "t"]
untreated_nottcell <- sum(mat[, "f"]) - untreated_tcell
treated_nottcell <- sum(mat[, "t"]) - treated_tcell

fisher.test(matrix(c(untreated_tcell, untreated_nottcell,
                     treated_tcell, treated_nottcell), nrow=2, byrow=T))
# Myeloid-centric
treated_myeloid <- mat["Myeloid_centric", "t"]
untreated_myeloid <- mat["Myeloid_centric", "f"]
treated_notmyeloid <- sum(mat[, "t"]) - treated_myeloid
untreated_notmyeloid <- sum(mat[, "f"]) - untreated_myeloid
fisher.test(matrix(c(treated_myeloid, treated_notmyeloid,
                     untreated_myeloid, untreated_notmyeloid), nrow=2, byrow=T))
# Fisher tests for TME archetypes in treated vs untreated samples
mat <- table(tme_samples$sex,tme_samples$tme_archetype_group)

fisher.test(mat[c(2, 3), ], simulate.p.value = TRUE, B = 10000)
fisher.test(t(mat)[c(2, 3), ], simulate.p.value = TRUE, B = 10000)
# Fisher tests for TME archetypes in adult vs pediatric samples
tme_samples$age_group <- tme_samples %>%
  mutate(age_group = factor(
      case_when(
        age >= 0  & age <= 15 ~ "Pediatric",
        age >= 16 & age <= 39 ~ "Young adult",
        age >= 40 & age <= 64 ~ "Adult",
        age >= 65             ~ "Elderly",
        TRUE                  ~ "Unknown"
      ),
      levels = c("Pediatric", "Young adult", "Adult", "Elderly", "Unknown")
    )) %>%
  pull(age_group)

mat <- table(tme_samples$age_group, tme_samples$tme_archetype_group)
fisher.test(mat[1:4, ], simulate.p.value = TRUE, B = 10000)
