library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(ComplexUpset)

setwd("Documents/bc-meta/sctherapy/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

subclone_annot <- read.table("annotations_subclones.tsv")
subclone_annot$cluster <- as.factor(subclone_annot$cluster)


######################## CLINICAL CHARACTERIZATION #############################

stacked_barplot <- function(df, x, fill, title, colors, dir){
  barplot <- ggplot(df,
                    aes(x = get(x), fill = get(fill))) +
     geom_bar(position = "fill") +
     scale_fill_manual(values = colors) +
     labs(x = "Cancer type", y = "Sample fraction", fill = title) +
     ggtitle("TME archetypes across cancer types") +
     theme_bw() +
     theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
           axis.title.x = element_text(size = 14, margin = margin(t = 6)),
           axis.title.y = element_text(size = 14, margin = margin(r = 6)),
           axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
           axis.text.y = element_text(size = 12, color = "black"),
           legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12)) +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x))
    # guides(fill = guide_legend(ncol = 2))
  
  
  ggsave(
    paste0(dir, "/", fill,"_barplot.png"),
    barplot,
    width = 14,
    height = 8,
    dpi = 500
  )
}

# Sex barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "sex", 
                title = "Chromosomal sex",
                colors = sex_colors,
                dir = "figures")

# Age barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "adult_pediatric", 
                title = "Age",
                colors = age_colors,
                dir = "figures")

# Solid/Liquid barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "is_blood", 
                title = "Solid/Liquid",
                colors = sl_colors,
                dir = "figures")

# Tumor site barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "summarised_tumor_site", 
                title = "Sample type",
                colors = tumor_sites_colors,
                dir = "figures")

# Check for how many of the brain tumor sites are from primary and metastasis in each cluster
subclone_brain <- subclone_annot %>%
  filter(summarised_tumor_site == "Brain")
stacked_barplot(subclone_brain, 
                x = "cluster", 
                fill = "sample_type", 
                title = "Sample type",
                colors = pm_colors,
                dir = "figures")

# Sample type condition barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "sample_type", 
                title = "Sample type",
                colors = pm_colors,
                dir = "figures")

# Treated condition barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "treated", 
                title = "Treatment",
                colors = treatment_colors,
                dir = "figures")
# TME barplot
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "tme_archetype", 
                title = "TME archetype",
                colors = tme_colors,
                dir = "figures")


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
stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "broad_cancer_type", 
                title = "Cancer type",
                colors = broad_cancer_type_colors,
                dir = "figures")

stacked_barplot(subclone_annot, 
                x = "cluster", 
                fill = "study", 
                title = "Study",
                colors = study_colors,
                dir = "figures")

# Compute the distribution of TME archetypes across cancer types over the total 
# number of samples
tcca_annot <- read.table("../cohort_statistics/tcca_annotation_raw.tsv", 
                         header = TRUE)
sample_annot <- tcca_annot %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  select(study_sample, tumor_type, tme_archetype, sample_type, treated) %>%
  left_join(cancer_type, by = "tumor_type") %>%
  distinct() %>%
  mutate(broad_cancer_type = factor(broad_cancer_type, 
                                    levels = names(sort(table(broad_cancer_type), 
                                                        decreasing = TRUE))),
         tme_archetype = factor(tme_archetype, 
                                levels = c(setdiff(unique(tme_archetype), "none"), "none")),
         treated = ifelse(treated == "t", "Treated", ifelse(treated == "f", "Untreated", NA)),
         sample_type = ifelse(sample_type == "m", "Metastasis", ifelse(sample_type == "p", "Primary", NA)))

stacked_barplot(sample_annot, 
                x = "broad_cancer_type", 
                fill = "tme_archetype", 
                title = "TME archetype",
                colors = tme_colors,
                dir = "../cohort_statistics/figures/")

stacked_barplot(sample_annot, 
                x = "tme_archetype", 
                fill = "treated", 
                title = "Treatment",
                colors = treatment_colors,
                dir = "../cohort_statistics/figures/")


stacked_barplot(sample_annot, 
                x = "tme_archetype", 
                fill = "sample_type", 
                title = "Sample type",
                colors = pm_colors,
                dir = "../cohort_statistics/figures/")

######################## THERAPEUTIC CHARACTERIZATION ##########################
#drugInfo <- load("../drugInfo.RData")
data <- read.table("full_table_drug_prediction.tsv")
cluster <- readRDS("speclustering_reordered.rds")
cluster_counts <- table(cluster)
cluster <- data.frame(Subclone = names(cluster), cluster = cluster)

data <- data %>% 
  left_join(cluster, by = "Subclone")

# Contingency table: count the number of unique subclones for each drug prediction within a cluster
contingency_table <- data %>%
  group_by(cluster, Drug_Name) %>%
  summarise(subclone_counts = n_distinct(Subclone), .groups = "drop") %>%
  pivot_wider(names_from = Drug_Name, values_from = subclone_counts, values_fill = 0) %>%
  as.data.frame()

bubble_data <- data %>%
  group_by(cluster, Drug_Name) %>%
  summarise(subclone_counts = n_distinct(Subclone), .groups = "drop") 

# Get unique clusters and drug names
unique_clusters <- unique(bubble_data$cluster)
unique_drugs <- unique(bubble_data$Drug_Name)

# Create all possible combinations of clusters and drug names
all_combinations <- expand.grid(cluster = unique_clusters, Drug_Name = unique_drugs)

# Merge with bubble_data to include zeros for missing combinations
bubble_data <- all_combinations %>%
  left_join(bubble_data, by = c("cluster", "Drug_Name")) %>%
  replace_na(list(subclone_counts = 0))  # Replace NA with 0 for missing counts

#  Run the Chi-Squared/Fisher test for each drug
independence_test <- function(freq_table, clusters, cluster_counts, pairwise = FALSE) {
  results <- data.frame(matrix(ncol = 4, nrow = 0))  # Store results
  colnames(results) <- c("Cluster1", "Cluster2", "Drug", "p_value")
  freq_table <- freq_table[-1]
  all_p_values <- c()  # Store all p-values for FDR
  pairwise_keys <- c()  # Store identifiers for comparisons
  if (pairwise) {
    # Pairwise comparisons: Compare each cluster to every other cluster
    cluster_combinations <- combn(clusters, 2, simplify = FALSE)
  } else {
    # Compare each cluster to the rest of the clusters
    cluster_combinations <- lapply(clusters, function(k)
      list(k, setdiff(clusters, k)))
  }
  
  for (pair in cluster_combinations) {
    print(pair)
    cluster1 <- pair[[1]]
    cluster2 <- pair[[2]]
    
    print(paste(
      "Comparing Cluster",
      cluster1,
      "to",
      if (pairwise)
        paste("Cluster", cluster2)
      else
        "Rest"
    ))
    
    for (drug in colnames(freq_table)) {
      # Assuming first column is cluster labels
      # Get contingency table values
      subclones_in_with_drug <- freq_table[cluster1, drug]
      subclones_in_without_drug <- cluster_counts[cluster1] - subclones_in_with_drug
      
      if (pairwise) {
        subclones_out_with_drug <- freq_table[cluster2, drug]
        subclones_out_without_drug <- cluster_counts[cluster2] - subclones_out_with_drug
      } else {
        subclones_out_with_drug <- sum(freq_table[cluster2, drug])
        subclones_out_without_drug <- sum(cluster_counts[cluster2]) - subclones_out_with_drug
      }
      
      # Construct contingency table
      contab_per_cluster_drug <- matrix(
        c(
          subclones_in_with_drug,
          subclones_in_without_drug,
          subclones_out_with_drug,
          subclones_out_without_drug
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      # Choose test based on expected counts
      expected_counts <- chisq.test(contab_per_cluster_drug)$expected
      use_fisher <- any(expected_counts < 5)
      
      if (use_fisher) {
        test_result <- fisher.test(contab_per_cluster_drug)
      } else {
        test_result <- chisq.test(contab_per_cluster_drug)
      }
      
      p_value <- signif(test_result$p.value, 3)
      
      # Store p-values for FDR adjustment
      all_p_values <- c(all_p_values, p_value)
      pairwise_keys <- c(pairwise_keys, paste(cluster1, cluster2, drug, sep = "_"))
      
      # Store results
      cluster_ref <- ifelse(pairwise, cluster2, "Rest")
      
      # Create a data frame for new result and ensure correct structure
      new_result <- data.frame(Cluster1 = cluster1, 
                               Cluster2 = cluster_ref, 
                               Drug = drug, 
                               p_value = p_value, 
                               stringsAsFactors = FALSE)
      print(cluster2)
      results <- rbind(results, new_result)
    }
  }
  # Apply FDR correction
  adjusted_p_values <- p.adjust(all_p_values, method = "BH")
  
  # Store FDR-adjusted p-values
  results$FDR_Adjusted_P_value <- adjusted_p_values
  return(results)
}



## Bubble plot of the results
bubble_stats <- independence_test(contingency_table, 
                                  contingency_table$cluster, 
                                  cluster_counts, 
                                  pairwise = FALSE)
# Transform the FDR = 0
bubble_stats <- bubble_stats  %>%
  dplyr::mutate(FDR_Adjusted_P_value = dplyr::if_else(FDR_Adjusted_P_value == 0 , 
                                                      10e-20,
                                                      FDR_Adjusted_P_value))
# Significance variable
bubble_stats <- bubble_stats %>%
  dplyr::mutate(Significance = dplyr::if_else(FDR_Adjusted_P_value < 0.05, TRUE,
                                              FALSE))
# Join to the table with subclone counts
bubble_data <- bubble_data %>%
  left_join(bubble_stats, by = c("cluster" = "Cluster1", "Drug_Name" = "Drug"))


# Bubble plot theme
bubble.theme <- list(
  ggplot2::theme_classic(),
  ggplot2:: theme(axis.ticks.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 45,
                                                      hjust = 1),
                  axis.ticks.y = ggplot2::element_blank(),
                  axis.text.y = ggplot2::element_text(vjust = 0.5, hjust = 0),
                  axis.line = ggplot2::element_blank(),
                  legend.box.just = "bottom",
                  plot.margin = unit(c(1, 1, 0, 1), "inches")),
  ggplot2::scale_y_discrete(position = "right"),
  ggplot2::labs(x = NULL, y = NULL),
  ggplot2::scale_size_continuous(range = c(1, 10)))

# Bubble plot.
min.val <- min(bubble_data$subclone_counts, 0)
max.val <- max(bubble_data$subclone_counts, 0)

bubble <- ggplot2::ggplot(bubble_data, ggplot2::aes(x = Drug_Name, y = cluster)) +
  ggplot2::geom_point(aes(color = subclone_counts, size = pmin(-log10(FDR_Adjusted_P_value), 20)), shape = 16) +
  ggplot2::geom_point(data = subset(bubble_data, !Significance), stroke = 1, shape = 21,
                      ggplot2::aes(size = pmin(-log10(FDR_Adjusted_P_value), 20), color = subclone_counts,
                                   fill = Significance)) +
  # #ggplot2::scale_color_gradientn(colours = color.gradient, guide = "colorbar",
  #                                values = scales::rescale(unique(c(min.val, 0,
  #                                                                  max.val))),
  #                                limits = c(min.val, max.val)) +
  labs(size = "-log10(FDR) Score") +
  scale_color_viridis_c(option = "viridis") +
  ggplot2::scale_fill_manual(values = "white",
                             labels = paste("FDR >= ", 0.05)) +
  ggplot2::guides(colour = ggplot2::guide_colorbar(order = 1),
                  size = ggplot2::guide_legend(order = 2),
                  fill = ggplot2::guide_legend(
                    order = 3,
                    override.aes = list(size = 10))) +
  bubble.theme

ggsave("bubble_plot_scTherapy.png", bubble, height = 10, width = 33, dpi = 300)

# Add MoAs to the drugs
MoAs <- read.table("../bc-meta_repo/bc-meta/reference/final_moas - Collapsed.tsv", 
                   header = TRUE,
                   sep = "\t") %>%
  select(preferred.drug.names, collapsed.MoAs) %>%
  distinct()

MoAs_sctherapy <- data %>%
  select(Drug_Name, MoA) %>%
  distinct() %>%
  mutate(Drug_Name = toupper(Drug_Name)) %>%
  column_to_rownames(var = "Drug_Name")

bubble_data <- bubble_data %>%
  mutate(Drug_Name = toupper(Drug_Name)) %>%
  left_join(select(MoAs, preferred.drug.names, collapsed.MoAs),
            by = c("Drug_Name" = "preferred.drug.names"))
            
drug_na <- unique(bubble_data$Drug_Name[is.na(bubble_data$collapsed.MoAs)])

bubble_data <- bubble_data %>%
  mutate(collapsed.MoAs = ifelse(is.na(collapsed.MoAs), 
                                 MoAs_sctherapy[Drug_Name, "MoA"],
                                 collapsed.MoAs))

bubble_data <- bubble_data %>%
  mutate(
    collapsed.MoAs = case_when(
      collapsed.MoAs == "thymidylate_synthase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "tubulin_polymerization_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "MDM_inhibitor" ~ "p53 activator/MDM2 inhibitor",
      collapsed.MoAs == "ALK_tyrosine_kinase_receptor_inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "Brutons_tyrosine_kinase_BTK_inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "inosine_monophosphate_dehydrogenase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "RAF_inhibitor" ~ "BRAF inhibitor",
      collapsed.MoAs == "RAF_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "BRD-K24576554" ~ "Other",
      collapsed.MoAs == "BRD-K60237333" ~ "Other",
      collapsed.MoAs == "BRD-K95142244" ~ "Other",
      collapsed.MoAs == "MAP_kinase_inhibitor" ~ "MAPK inhibitor",
      collapsed.MoAs == "adenosine_receptor_antagonist" ~ "Other",
      collapsed.MoAs == "centromere_associated_protein_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "PLK_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "HDAC_inhibitor" ~ "Chromatin agent",
      collapsed.MoAs == "CDK_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "topoisomerase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "DNA_synthesis_inhibitor_microtubule_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "kinesin_inhibitor_kinesinlike_spindle_protein_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "VEGFR_inhibitor" ~ "VEGFR inhibitor",
      collapsed.MoAs == "potassium_channel_activator" ~ "Other",
      collapsed.MoAs == "KIT_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor_VEGFR_inhibitor" ~ "VEGFR inhibitor",
      collapsed.MoAs == "dehydrogenase_inhibitor_inositol_monophosphatase_inhibitor" ~ "Other",
      collapsed.MoAs == "FLT3_inhibitor_JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "Ttype_calcium_channel_blocker" ~ "Other",
      collapsed.MoAs == "BcrAbl_kinase_inhibitor_FLT3_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor" ~ "BCR-ABL inhibitor",
      collapsed.MoAs == "estrogen_receptor_antagonist_selective_estrogen_receptor_modulator_SERM" ~ "Other",
      collapsed.MoAs == "retinoid_receptor_agonist" ~ "Other",
      collapsed.MoAs == "HSP_inhibitor" ~ "HSP inhibitor",
      collapsed.MoAs == "RNA_polymerase_inhibitor" ~ "Transcription inhibitor",
      collapsed.MoAs == "PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
      collapsed.MoAs == "JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "mTOR_inhibitor_PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
      collapsed.MoAs == "Aurora_kinase_inhibitor_FLT3_inhibitor_VEGFR_inhibitor" ~ "Other",
      Drug_Name == "CYT387" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "-" ~ "Other",
      TRUE ~ collapsed.MoAs
    )
)

# Bubble plot (each cluster against the rest) splitting drugs by the MoAs
bubble <- ggplot2::ggplot(bubble_data, ggplot2::aes(x = Drug_Name, y = cluster)) +
ggplot2::geom_point(aes(color = subclone_counts, size = pmin(-log10(FDR_Adjusted_P_value), 20)), shape = 16) +
  ggh4x::facet_grid2(~ collapsed.MoAs, scales = "free_x", space = "free_x") +
  ggplot2::geom_point(data = subset(bubble_data, !Significance), stroke = 1, shape = 21,
                      ggplot2::aes(size = pmin(-log10(FDR_Adjusted_P_value), 20), color = subclone_counts,
                                   fill = Significance)) +
  # #ggplot2::scale_color_gradientn(colours = color.gradient, guide = "colorbar",
  #                                values = scales::rescale(unique(c(min.val, 0,
  #                                                                  max.val))),
  #                                limits = c(min.val, max.val)) +
  labs(size = "-log10(FDR) Score") +
  scale_color_viridis_c(option = "viridis") +
  ggplot2::scale_fill_manual(values = "white",
                             labels = paste("FDR >= ", 0.05)) +
  ggplot2::guides(colour = ggplot2::guide_colorbar(order = 1),
                  size = ggplot2::guide_legend(order = 2),
                  fill = ggplot2::guide_legend(
                    order = 3,
                    override.aes = list(size = 10))) +
  bubble.theme

ggsave("bubble_split_MoAscTherapy.png", bubble, height = 15, width = 33, dpi = 300)

library(ggplot2)

# Create separate plots for each `collapsed.MoAs` and save them independently
unique_moas <- unique(bubble_data$collapsed.MoAs)
subclone_df <- as.data.frame(table(cluster$cluster))
bubble_data <- bubble_data %>%
  left_join(subclone_df, by = c("cluster" = "Var1")) %>%
  mutate(subclones_proportion = round((subclone_counts/Freq) *100, 2))

for (moa in unique_moas) {
  # Filter data for the current MoA
  data_subset <- subset(bubble_data, collapsed.MoAs == moa)
  
  # Create the individual bubble plot
  bubble_plot <- ggplot(data_subset, aes(x = Drug_Name, y = cluster)) +
    geom_point(aes(color = subclones_proportion, size = pmin(-log10(FDR_Adjusted_P_value), 20)), shape = 16) +
    geom_point(data = subset(data_subset, !Significance), stroke = 1, shape = 21,
               aes(size = pmin(-log10(FDR_Adjusted_P_value), 20), color = subclones_proportion, fill = Significance)) +
    labs(title = moa, size = "-log10(FDR) Score") +
    scale_color_viridis_c(option = "viridis", limits = c(0, max(bubble_data$subclones_proportion))) +
    scale_fill_manual(values = "white", labels = paste("FDR >= ", 0.05)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
  
  # Define filename
  filename <- paste0("figures/bubble_", gsub("[ /]", "_", gsub("/", "|", moa)), ".png")
  
  # Save the plot
  ggsave(filename, bubble_plot, height = 5, width = max(5, length(unique(data_subset$Drug_Name)) * 0.6), dpi = 300)
  
  print(paste("Saved:", filename))  # Print message for confirmation
}



# Select drugs significantly associated with clusters (FDR <= 0.05 and >= 10% of subclones)
significant_drugs <- bubble_data %>%
  filter(FDR_Adjusted_P_value <= 0.05 & subclones_proportion >= 10)

# Prepare data for an UpSet plot
drugs_upset <- significant_drugs %>%
  select(cluster, Drug_Name, collapsed.MoAs) %>%
  mutate(Presence = 1) %>%
  distinct(Drug_Name, collapsed.MoAs, cluster, .keep_all = TRUE) %>%
  pivot_wider(names_from = cluster, values_from = Presence, values_fill = list(Presence = 0))

clusters <- colnames(drugs_upset)[3: ncol(drugs_upset)]
drugs_upset[clusters] <- drugs_upset[clusters] == 1
t(head(drugs_upset[clusters], 3))


set_size(8, 3)
source("/home/lmgonzalezb/Documents/bc-meta/bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")
upset(drugs_upset, 
      clusters, 
      mode = "inclusive_intersection", 
      name = "clusters", 
      base_annotations=list(
        'Intersection size' = intersection_size(
          counts=FALSE,
          mapping = aes(fill = collapsed.MoAs)
        )
      ),
      width_ratio = 0.1)

drugs_upset$collapsed.MoAs <- ifelse(
  drugs_upset$collapsed.MoAs == "BRAF inhibitor;VEGFR inhibitor", "Kinase inhibitor",
  ifelse(drugs_upset$collapsed.MoAs == "BCR-ABL inhibitor;SRC inhibitor", "BCR-ABL inhibitor", 
         ifelse(drugs_upset$collapsed.MoAs == "VEGFR inhibitor;MET inhibitor", "Multi-kinase inhibitor",
                drugs_upset$collapsed.MoAs)
))

unique(drugs_upset$collapsed.MoAs) %in% names(MoAs_colors)
