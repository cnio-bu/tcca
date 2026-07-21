library(BPCells)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(GSEABase)
library(UCell)
library(dplyr)
library(readxl)
library(stringr)

setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/therapeutic_analysis/gdsc")
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")


# Get cancer cell lines expression data
mat <- open_matrix_dir("cell_lines_gabriella_kinker_v5/")
metadata <- read.table("tcca_metadata_h5ad.tsv", header = T, sep = "\t") %>%
  column_to_rownames("cell")
seu <- CreateSeuratObject(mat, meta.data = metadata)


# Load cancer cell lines drug response data (GDSC2)
gdsc <- read_excel("GDSC2_fitted_dose_response_27Oct23.xlsx") %>%
  mutate(CELL_LINE_NAME = sub("-", "", CELL_LINE_NAME)) %>%
  dplyr::select(CELL_LINE_NAME, DRUG_NAME, PATHWAY_NAME, AUC) %>%
  group_by(CELL_LINE_NAME, DRUG_NAME, PATHWAY_NAME) %>%
  summarise(AUC = mean(AUC), .groups = "drop")## These steps averages duplicated AUC (same drug, different ID)


# Get common cell lines between expression and drug response data
seu@meta.data <- seu@meta.data %>%
  mutate(CELL_LINE_NAME = sub("_.*", "", sample))
common_cell_lines <- intersect(seu$CELL_LINE_NAME, gdsc$CELL_LINE_NAME)
seu <- subset(seu, subset = CELL_LINE_NAME %in% common_cell_lines)
gdsc <- subset(gdsc, subset = CELL_LINE_NAME %in% common_cell_lines)


# Read TC marker genes gmt file
markers <- getGmt("../sctherapy/marker_genes/survival_results/marker_sigs_filtered.gmt")
markers <- geneIds(markers)


# Adapt signatures to keep only genes expressed in the dataset
expr_matrix <- GetAssayData(seu, slot = "counts")
expr_genes <- rownames(expr_matrix)[Matrix::rowSums(expr_matrix > 0) > 0]
markers <- lapply(markers, function(genes_vector) {
  intersect(genes_vector, expr_genes)
})
seu <- AddModuleScore_UCell(seu, features=markers, name = NULL)

# Compute mean TC markers enrichment score for each cell line
ccl_markers <- seu@meta.data %>%
  dplyr::select(CELL_LINE_NAME, tumor_type, starts_with("Cluster")) %>%
  group_by(CELL_LINE_NAME, tumor_type) %>%
  summarise(across(starts_with("Cluster"), mean, na.rm = TRUE), .groups = "drop")

colnames(ccl_markers) <- c("CELL_LINE_NAME", "tumor_type", paste0(rep("TC", 10), 1:10))

ccl_markers <- pivot_longer(ccl_markers, cols = starts_with("TC"), names_to = "TC", values_to = "UCell_score")

# Join TC markers enrichment with drug data in a single table
global_data <- ccl_markers %>%
    left_join(gdsc, by = "CELL_LINE_NAME", relationship = "many-to-many") %>%
    dplyr::rename(TC_score = UCell_score) %>%
    relocate(CELL_LINE_NAME, tumor_type, DRUG_NAME, PATHWAY_NAME, AUC, TC, TC_score)


# Add a column indicating whether the drug is a top drug for the corresponding TC cluster
predicted_drugs <- read.table("../../cohort_statistics/drug_response_subclone_final.tsv", header = TRUE) %>%
    distinct(ScTherapy.Cluster, Drug.Name) %>%
    group_by(ScTherapy.Cluster) %>%
    mutate(row = row_number()) %>%
    pivot_wider(
        names_from = ScTherapy.Cluster,
        values_from = Drug.Name
    ) %>%
    dplyr::select(-row)
predicted_drugs <- predicted_drugs[, order(as.numeric(colnames(predicted_drugs)))]
colnames(predicted_drugs) <- paste0(rep("TC"), 1:10)
global_data <- global_data %>%
    mutate(DRUG_NAME = str_to_upper(DRUG_NAME)) %>%
    mutate(predicted_in_TC = case_when(
        TC == "TC1" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC1, TRUE, FALSE),
        TC == "TC2" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC2, TRUE, FALSE),
        TC == "TC3" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC3, TRUE, FALSE),
        TC == "TC4" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC4, TRUE, FALSE),
        TC == "TC5" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC5, TRUE, FALSE),
        TC == "TC6" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC6, TRUE, FALSE),
        TC == "TC7" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC7, TRUE, FALSE),
        TC == "TC8" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC8, TRUE, FALSE),
        TC == "TC9" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC9, TRUE, FALSE),
        TC == "TC10" ~ ifelse(DRUG_NAME %in% predicted_drugs$TC10, TRUE, FALSE),
        TRUE ~ FALSE
    ))

# Filter cancer types with less than 3 lines to compute the correlations
strata_ok <- global_data %>%
  group_by(TC, tumor_type) %>%
  summarise(n_lines = n_distinct(CELL_LINE_NAME), .groups = "drop") %>%
  filter(n_lines >= 8)

global_data_filtered <- global_data %>% semi_join(strata_ok, by = c("TC", "tumor_type"))

# # For a TC consider only cancer types represented in the TC according to scTherapy analysis
# tc_tumor_types <- read.table("../../cohort_statistics/drug_response_subclone_final.tsv", header = T, sep = "\t") %>%
#   dplyr::select(ScTherapy.Cluster, Refined.Tumor.Type) %>%
#   distinct()
# tumor_list <- tc_tumor_types %>%
#   group_by(ScTherapy.Cluster) %>%
#   summarise(clusters = list(Refined.Tumor.Type), .groups = "drop") %>%
#   tibble::deframe()
# names(tumor_list) <- paste0("TC", names(tumor_list))

# # Filter global_data_filtered to keep only tumor types represented in the TC according to scTherapy analysis
# global_data_filtered <- global_data_filtered %>%
#   rowwise() %>%
#   filter(tumor_type %in% tumor_list[[TC]]) %>%
#   ungroup()
  
# Correlation per TC x tumor_type x drug
# Spearman is safer than Pearson here: AUC-TC_score relationship is not guaranteed to be linear, and Spearman is robust to outlier cell lines
cor_strata <- global_data_filtered %>%
  group_by(TC, tumor_type, DRUG_NAME, PATHWAY_NAME, predicted_in_TC) %>%
  summarise(
    n  = n(),
    ct = list(suppressWarnings(cor.test(TC_score, AUC, method = "spearman"))),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(ct, ~ unname(.x$estimate)),
    p = map_dbl(ct, ~ .x$p.value)
  ) %>%
  dplyr::select(-ct) %>%
  filter(!is.na(r))

# # Meta-analysis across tumor-type strata (Fisher's z)
# fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
# fisher_z_inv <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)

# cor_drug_tc <- cor_strata %>%
#   mutate(
#     r_clamped = pmin(pmax(r, -0.999), 0.999),
#     z = fisher_z(r_clamped),
#     w = pmax(n - 3, 1)
#   ) %>%
#   group_by(TC, DRUG_NAME, PATHWAY_NAME) %>%
#   summarise(
#     z_meta   = weighted.mean(z, w),
#     n_strata = n(),
#     n_lines_total = sum(n),
#     .groups  = "drop"
#   ) %>%
#   mutate(
#     r_meta = fisher_z_inv(z_meta),
#     drug_norm = str_to_upper(DRUG_NAME)
#   )

# Percentage of predicted drugs per TC x drug
subclones_predictions <- read.table("../../cohort_statistics/drug_response_subclone_final.tsv", header = TRUE)
n_subclones_per_tc_cancertype <- subclones_predictions %>%
  distinct(ScTherapy.Cluster, Refined.Tumor.Type, Subclone.Name) %>% 
  count(ScTherapy.Cluster, Refined.Tumor.Type, name = "n_subclones")

pct_predicted_per_tc_cancertype_drug <- subclones_predictions %>%
  count(ScTherapy.Cluster, Refined.Tumor.Type, Drug.Name, name = "n_predicted") %>%
  left_join(n_subclones_per_tc_cancertype,
            by = c("ScTherapy.Cluster", "Refined.Tumor.Type")) %>%
  mutate(
    pct_subclones_predicted = n_predicted / n_subclones * 100,
    drug_norm = str_to_upper(Drug.Name),
    ScTherapy.Cluster = paste0("TC", ScTherapy.Cluster)
  ) %>%
  dplyr::rename(TC = ScTherapy.Cluster)

# Join final: concordance prediction (scTherapy) <-> validation (GDSC)
validation_continuous <- cor_strata %>%
    inner_join(
        pct_predicted_per_tc_cancertype_drug,
        by = c("TC", "DRUG_NAME" = "drug_norm", "tumor_type" = "Refined.Tumor.Type"),
        suffix = c("_gdsc", "_sctherapy")
    )

validation_continuous %>% count(TC) %>% print()

# Prepare data for plotting
plot_df <- validation_continuous %>%
  filter(predicted_in_TC) %>%
  mutate(
    sig_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )
  ) %>%
  mutate(
    TC = factor(TC, levels = paste0("TC", 1:10)),
    drug_ordered = tidytext::reorder_within(DRUG_NAME, -r, tumor_type)
  )

# Fixed order of drugs, grouping by pathway/MoA
drug_order <- plot_df %>%
  distinct(DRUG_NAME, PATHWAY_NAME) %>%
  arrange(PATHWAY_NAME, DRUG_NAME) %>%
  pull(DRUG_NAME)

plot_df <- plot_df %>%
  mutate(
    TC = factor(TC, levels = paste0("TC", rev(1:10))),
    DRUG_NAME = factor(DRUG_NAME, levels = drug_order)
  )

bubble_plot <- ggplot(plot_df, aes(x = DRUG_NAME, y = TC)) +
    geom_point(aes(size = pct_subclones_predicted, fill = r),
                shape = 21, color = "grey40", stroke = 0.4) +
    geom_text(aes(label = sig_label),
                color = "black", size = 3.5, fontface = "bold", vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                        name = "Spearman ρ") +
    scale_size_continuous(name = "% Subclones\npredicted to \nbe sensitive", range = c(2, 9)) +
    facet_wrap(~ tumor_type, ncol = 1) +
    labs(x = NULL, y = "Therapeutic Cluster (TC)") +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey70"),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold")
    )
ggsave("bubble_plot_correlations.png", bubble_plot, width = 14, height = 10, dpi = 300)


# Barplot
summary_by_tc <- validation_continuous %>%
  group_by(TC, tumor_type) %>%
  summarise(
    n_drugs = n(),
    n_negative = sum(r < 0),
    pct_negative = 100 * n_negative / n_drugs,
    n_sig_negative = sum(r < 0 & p < 0.05),
    .groups = "drop"
  ) %>%
  mutate(TC = factor(TC, levels = paste0("TC", 1:10)))

# Binomial test per TC-cancer_type: is % of negative correlations far way from the 50% expected by chance?
summary_by_tc <- summary_by_tc %>%
  rowwise() %>%
  mutate(
    p_binom = binom.test(n_negative, n_drugs, p = 0.5, alternative = "greater")$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_binom, method = "BH"),
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")
summary_plot <- ggplot(summary_by_tc, aes(x = TC, y = pct_negative, fill = tumor_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30", linewidth = 0.2) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
    geom_text(aes(label = sig_label, group = tumor_type),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 4) +
    scale_fill_manual(values = tumor_type_colors) +
    labs(
        x = "Therapeutic Cluster (TC)",
        y = "% of drugs with Spearman ρ < 0\n(TC signature vs. AUC)",
        fill = "Tumor type",
        title = "Prediction–validation concordance across therapeutic clusters",
        caption = "Dashed line marks the 50% expectation under the null. Asterisks indicate BH-adjusted one-sided binomial test significance."
        ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))

ggsave("summary_plot_validation.pdf", summary_plot, width = 7, height = 6, dpi = 300)


# Fisher z-transformation and inverse transformation for correlation coefficients.
# Correlations are transformed to an approximately normal scale before meta-analysis.
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
fisher_z_inv <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)

# Pan-cancer meta-analysis of prediction–validation correlations.
# For each Therapeutic Cluster (TC):
#   1. Transform correlation coefficients (r) to Fisher z values.
#   2. Compute a fixed-effect meta-analysis using inverse-variance weights
#      (w = n - 3, where Var(z) ≈ 1/(n - 3)).
#   3. Estimate the combined effect size (z_meta) and its standard error.
#   4. Test whether the pooled correlation differs from zero.
#   5. Adjust p-values across TCs using the Benjamini–Hochberg procedure.

summary_pancancer_sig <- validation_continuous %>%
  mutate(
    # Avoid infinite Fisher z values when r approaches ±1
    r_clamped = pmin(pmax(r, -0.999), 0.999),
    # Fisher z-transformed correlation
    z = fisher_z(r_clamped),
    # Fixed-effect meta-analysis weight
    w = pmax(n - 3, 1)
  ) %>%
  group_by(TC) %>%
  summarise(
    # Weighted mean Fisher z across drugs
    z_meta = weighted.mean(z, w),
    # Standard error of the pooled Fisher z
    se_meta = 1 / sqrt(sum(w)),
    # Number of drug-level observations contributing to the estimate
    n_drugs = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Convert pooled Fisher z back to the correlation scale
    r_meta = fisher_z_inv(z_meta),
    # Wald test statistic under H0: pooled correlation = 0
    z_stat = z_meta / se_meta,
    # Two-sided p-value from the standard normal distribution
    p_value = 2 * pnorm(-abs(z_stat)),
    # Benjamini–Hochberg FDR correction across TCs
    p_adj = p.adjust(p_value, method = "BH"),
    TC  = factor(TC, levels = paste0("TC", 1:10)),
    # Significance annotation for plotting
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Display TCs ranked by statistical significance
print(summary_pancancer_sig %>% arrange(p_adj), n = 10)

# Visualization of pooled prediction–validation correlations across TCs
ggplot(summary_pancancer_sig, aes(x = TC, y = r_meta, fill = r_meta)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sig_label, vjust = ifelse(r_meta >= 0, -0.3, 1.2)), size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                        name = "r meta") +
  labs(x = "Therapeutic Cluster", y = "Meta-analyzed correlation (SKCM, BRCA, and LUAD)",
       title = "Pan-cancer prediction–validation concordance across therapeutic clusters") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")