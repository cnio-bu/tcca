library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggstatsplot)
library(ComplexHeatmap)
library(ggpubr)

setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/cohort_statistics/")
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")

metadata <- read.table("tcca_metadata_h5ad.tsv", sep = "\t", header = TRUE)

# Exclude TME cells and samples (~5%) with unsuccessful CNV inference by SCEVAN
samples_with_subclones <- metadata %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    filter(malignancy == "True" & !(scevan_subclone %in% c("non_tumor", "nan")))

# Compute number of subclones and clonal diversity for each sample
samples_with_subclones <- samples_with_subclones %>%
    group_by(
        study_sample, patient, tumor_type, sex, age, treated,
        sample_type, tme_archetype_group, scevan_subclone
    ) %>%
    summarise(n_cells_subclone = n(), .groups = "drop") %>%
    group_by(
        study_sample, patient, tumor_type, sex, age, treated,
        sample_type, tme_archetype_group
    ) %>%
    mutate(freq = n_cells_subclone / sum(n_cells_subclone)) %>%
    summarise(
        clonal_dominance = max(freq), # proportion of the most abundant subclone
        shannon_entropy = -sum(freq * log2(freq)), # clonal diversity index
        simpson_index = sum(freq^2), # Simpson's index
        inverse_simpson = 1 / sum(freq^2), # clonal diversity based on inverse simpson
        n_subclones = n_distinct(scevan_subclone),
        n_cells = sum(n_cells_subclone),
        .groups = "drop"
    ) %>%
    as.data.frame()


# Plot clonal diversity for each sample across cancer types, excluding cell line samples
clonality_patients <- samples_with_subclones %>%
    filter(patient != "ccl") %>%
    mutate(
        # Reorder tumor types based on clonal diversity index
        tumor_type = fct_reorder(
            tumor_type,
            inverse_simpson,
            .fun = median,
            .desc = TRUE
        ),

        # Convert sample type to factor with meaningful labels
        sample_type = factor(
            case_when(
                sample_type == "p" ~ "Primary",
                sample_type == "m" ~ "Metastasis"
            ),
            levels = c("Primary", "Metastasis")
        ),

        # Convert sex to factor with standardized labels
        sex = factor(
            case_when(
                sex == "m" ~ "Male",
                sex == "f" ~ "Female",
                sex == "nan" ~ "Unknown"
            ),
            levels = c("Male", "Female", "Unknown")
        ),

        # Convert treatment status to factor
        treated = factor(
            case_when(
                treated == "t" ~ "Treated",
                treated == "f" ~ "Untreated",
                TRUE ~ "Unknown"
            ),
            levels = c("Untreated", "Treated", "Unknown")
        ),

        # Group ages into categories
        age_group = factor(
            case_when(
                age >= 0 & age <= 15 ~ "Pediatric",
                age >= 16 & age <= 39 ~ "Young adult",
                age >= 40 & age <= 64 ~ "Adult",
                age >= 65 ~ "Elderly",
                TRUE ~ "Unknown"
            ),
            levels = c("Pediatric", "Young adult", "Adult", "Elderly", "Unknown")
        ),

        # Convert TME archetype group to factor with defined levels
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
        )
    )

# Boxplot
theme_size <- theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 14, margin = margin(t = 6), face = "bold"),
    axis.title.y = element_text(size = 14, margin = margin(r = 6), face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
)

# Check the correlation between the number of cells and the number of subclones
# per sample
cor_test <- cor.test(clonality_patients$n_cells,
    clonality_patients$n_subclones,
    method = "spearman"
)

# Create label with statistics for the plot
cor_label <- paste0(
    "Spearman ρ = ", round(cor_test$estimate, 2), "\n",
    "R² = ", round(cor_test$estimate^2, 2), "\n",
    "p < 0.001"
)

cor_plot <- ggplot(clonality_patients, aes(x = n_cells, y = n_subclones)) +
    geom_jitter(aes(color = tumor_type),
        alpha = 0.7,
        size = 2,
        height = 0.2,
        width = 0
    ) +
    annotate(
        "text",
        x = 10,
        y = 8.5,
        label = cor_label,
        size = 4,
        hjust = 0,
        color = "black"
    ) +
    scale_x_log10(
        breaks = c(100, 500, 1000, 5000, 10000, 20000),
        labels = scales::comma
    ) +
    annotation_logticks(sides = "b") +
    scale_color_manual(values = tumor_type_colors) +
    scale_y_continuous(
        breaks = 1:9,
        limits = c(0.5, 9)
    ) +
    labs(
        x = "Number of cells per sample",
        y = "Number of subclones per sample",
        color = "Tumor type"
    ) +
    theme_bw() +
    theme_size
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/correlation_cells_subclones.pdf",
    plot = cor_plot,
    width = 8,
    height = 6
)

boxplot <-
    clonality_patients %>%
    mutate(
        tumor_type = fct_reorder(tumor_type, n_subclones, .fun = median, .desc = TRUE)
    ) %>%
    ggplot(
        aes(
            x = tumor_type,
            y = n_subclones,
            fill = tumor_type
        )
    ) +
    geom_boxplot() +
    scale_fill_manual(values = tumor_type_colors) +
    labs(
        x = "Cancer type",
        y = "Number of subclones"
    ) +
    theme_bw() +
    theme_size +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/num_subclones_per_cancertype.pdf",
    plot = boxplot,
    width = 10,
    height = 6
)

# Boxplot of clonal diversity for each sample across cancer types
boxplot <- ggplot(
    clonality_patients,
    aes(
        x = tumor_type,
        y = inverse_simpson,
        fill = tumor_type
    )
) +
    geom_boxplot() +
    scale_fill_manual(values = tumor_type_colors) +
    labs(
        x = "Cancer type",
        y = "Clonal diversity index (Inverse Simpson)"
    ) +
    theme_bw() +
    theme_size +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonal_diversity_per_cancertype_simpson.pdf",
    plot = boxplot,
    width = 10,
    height = 6
)


# Plot number of subclones for clinical variables
boxplot_function <- function(
  df = clonality_patients,
  x,
  y = "inverse_simpson",
  fill,
  colors,
  x_lab,
  y_lab = "Clonal diversity index (Inverse Simpson)"
) {
    comp_list <- combn(as.character(unique(df[[x]])), 2, simplify = FALSE)
    # Filter out comparisons with the "Unknown" or "None" groups
    comp_list <- comp_list[!sapply(comp_list, function(pair) {
        any(pair %in% c("Unknown", "none", "NA"))
    })]
    boxplot <- ggplot(
        df,
        aes(
            x = !!sym(x),
            y = !!sym(y),
            fill = !!sym(fill)
        )
    ) +
        geom_boxplot() +
        scale_fill_manual(values = colors) +
        labs(
            x = x_lab,
            y = y_lab
        ) +
        theme_bw() +
        theme_size +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) + # add all pairwise comparisons
        stat_compare_means(
            method = "wilcox.test", # or "t.test"
            comparisons = comp_list,
            label = "p.signif" # "p.format" for exact values
        )
    return(boxplot)
}

# Sex
boxplot <- boxplot_function(x = "sex", fill = "sex", colors = sex_colors, x_lab = "Sex")
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonality_per_sex_simpson.pdf",
    plot = boxplot,
    width = 4,
    height = 6
)

# Age
age_group_colors <- c(
    "Pediatric" = "#a9dce3",
    "Young adult" = "#6CC2BD",
    "Adult" = "#7689DE",
    "Elderly" = "#4B4A73",
    "Unknown" = "#BEBEBE"
)
boxplot <- boxplot_function(
    x = "age_group",
    fill = "age_group",
    colors = age_group_colors,
    x_lab = "Age group"
)
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonality_age_group_simpson.pdf",
    plot = boxplot,
    width = 4,
    height = 6
)

# Primary/Metastasis
boxplot <- boxplot_function(
    x = "sample_type",
    fill = "sample_type",
    colors = pm_colors,
    x_lab = "Sample type"
)
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonality_sample_type_simpson.pdf",
    plot = boxplot,
    width = 2,
    height = 6
)

# Treated/untreated
boxplot <- boxplot_function(
    x = "treated",
    fill = "treated",
    colors = treatment_colors,
    x_lab = "Treatment condition"
)
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonality_treated_simpson.pdf",
    plot = boxplot,
    width = 3,
    height = 6
)

# TME group archetype
boxplot <- boxplot_function(
    x = "tme_archetype_group",
    fill = "tme_archetype_group",
    colors = tme_group_colors,
    x_lab = "TME group"
)
ggsave(
    "/Users/mariagb/Documents/new_figures_tcca/clonality_tme_group_simpson.pdf",
    plot = boxplot,
    width = 4,
    height = 6
)

# Print real p-values
# Sex comparison
sex_stats <- compare_means(
    inverse_simpson ~ sex,
    data = clonality_patients,
    method = "wilcox.test"
)
print(sex_stats)

# Age group comparison
age_stats <- compare_means(
    inverse_simpson ~ age_group,
    data = clonality_patients,
    method = "wilcox.test",
    ref.group = "Pediatric"
)
print(age_stats)

# Sample type comparison (Primary vs Metastasis)
sample_type_stats <- compare_means(
    inverse_simpson ~ sample_type,
    data = clonality_patients,
    method = "wilcox.test"
)
print(sample_type_stats)

# Treatment comparison
treatment_stats <- compare_means(
    inverse_simpson ~ treated,
    data = clonality_patients,
    method = "wilcox.test"
)
print(treatment_stats)

# TME archetype comparison
tme_stats <- compare_means(
    inverse_simpson ~ tme_archetype_group,
    data = clonality_patients,
    method = "wilcox.test"
)
print(tme_stats)
