library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggstatsplot)
library(ComplexHeatmap)
library(ggpubr)
  
setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/cohort_statistics/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)

# Exclude TME cells and samples (~5%) with unsuccessful CNV inference by SCEVAN
samples_with_subclones <- metadata %>% 
    filter(malignancy == "True" & !(scevan_subclone %in% c("", "non_tumor")))

# Compute number of subclones per 1000 cells for each sample
samples_with_subclones <- samples_with_subclones %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    group_by(study_sample) %>%
    mutate(subclones_per_1000_cells = n_distinct(scevan_subclone)/ n() * 1000,
        n_cells = n()) %>%
    ungroup() %>%
    as.data.frame() %>%
    select(
        sample, 
        patient, 
        study, 
        refined_tumor_type, 
        sex, 
        age, 
        treated, 
        sample_type, 
        tme_archetype_group,
        subclones_per_1000_cells, 
        n_cells) %>%
    distinct()
    

# Plot subclones per 1000 cells for each sample across cancer types, excluding cell line samples
clonality_patients <- samples_with_subclones %>%
  filter(patient != "ccl") %>%
  mutate(
    study_sample = paste0(study, "_", sample),
    # Reorder tumor types based on median number of subclones per 1000 cells
    refined_tumor_type = fct_reorder(
      refined_tumor_type,
      subclones_per_1000_cells,
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
        sex == ""  ~ "Unknown"
      ),
      levels = c("Male", "Female", "Unknown")
    ),
    
    # Convert treatment status to factor
    treated = factor(
      case_when(
        treated == "t" ~ "Treated",
        treated == "f" ~ "Untreated",
        TRUE           ~ "Unknown"
      ),
      levels = c("Untreated", "Treated", "Unknown")
    ),
    
    # Group ages into categories
    age_group = factor(
      case_when(
        age >= 0  & age <= 15 ~ "Pediatric",
        age >= 16 & age <= 39 ~ "Young adult",
        age >= 40 & age <= 64 ~ "Adult",
        age >= 65             ~ "Elderly",
        TRUE                  ~ "Unknown"
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

boxplot <- ggplot(
    clonality_patients, 
    aes(
        x = refined_tumor_type, 
        y = subclones_per_1000_cells, 
        fill = refined_tumor_type
        )
    ) +
    geom_boxplot()+
    scale_fill_manual(values = tumor_type_colors) +
    labs(
        x = "Cancer type",
        y = "Number of subclones per 1000 cells"
    ) +
    theme_bw() +
    theme_size +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
ggsave("../figures/clonality/clonality_per_cancertype.pdf", plot = boxplot, width = 10, height = 6)


# Plot number of subclones for clinical variables
boxplot_function <- function(
    df = clonality_patients,
    x, 
    y = "subclones_per_1000_cells", 
    fill, 
    colors, 
    x_lab, 
    y_lab = "Number of subclones per 1000 cells"
    ) {
    comp_list <- combn(as.character(unique(df[[x]]), 2, simplify = FALSE))
    boxplot <- ggplot(
        df, 
        aes(
            x = !!sym(x), 
            y = !!sym(y), 
            fill = !!sym(fill)
            )
        ) +
        geom_boxplot()+
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
        )  + # add all pairwise comparisons 
        stat_compare_means( 
          method = "wilcox.test", # or "t.test" 
          comparisons = comp_list, 
          label = "p.signif" # "p.format" for exact values 
          )
    return(boxplot)
}

# Sex
boxplot <- boxplot_function(x = "sex", fill = "sex", colors = sex_colors, x_lab = "Sex")
ggsave("../figures/clonality/clonality_per_sex.pdf", plot = boxplot, width = 4, height = 6)

# Age
age_group_colors <- c(
  "Pediatric"   = "#a9dce3",
  "Young adult" = "#6CC2BD",
  "Adult"       = "#7689DE",
  "Elderly"     = "#4B4A73",
  "Unknown"     = "#BEBEBE"
)
boxplot <- boxplot_function(x = "age_group", fill = "age_group", colors = age_group_colors, x_lab = "Age group")
ggsave("../figures/clonality/clonality_age_group.pdf", plot = boxplot, width = 4, height = 6)

# Primary/Metastasis
boxplot <- boxplot_function(x = "sample_type", fill = "sample_type", colors = pm_colors, x_lab = "Sample type")
ggsave("../figures/clonality/clonality_sample_type.pdf", plot = boxplot, width = 4, height = 6)

# Treated/untreated
boxplot <- boxplot_function(x = "treated", fill = "treated", colors = treatment_colors, x_lab = "Treatment condition")
ggsave("../figures/clonality/clonality_treated.pdf", plot = boxplot, width = 4, height = 6)

# TME group archetype
boxplot <- boxplot_function(x = "tme_archetype_group", fill = "tme_archetype_group", colors = tme_group_colors, x_lab = "TME group")
ggsave("../figures/clonality/clonality_tme_group.pdf", plot = boxplot, width = 4, height = 6)

# Print real p-values
library(ggpubr)
library(dplyr)

# Comparación por sexo
sex_stats <- compare_means(
  subclones_per_1000_cells ~ sex,
  data = clonality_patients,
  method = "wilcox.test"
)
print(sex_stats)

# Comparación por edad
age_stats <- compare_means(
  subclones_per_1000_cells ~ age_group,
  data = clonality_patients,
  method = "wilcox.test",
  ref.group = "Pediatric" # opcional si quieres comparar contra Pediatric
)
print(age_stats)

# Comparación por sample type (Primary vs Metastasis)
sample_type_stats <- compare_means(
  subclones_per_1000_cells ~ sample_type,
  data = clonality_patients,
  method = "wilcox.test"
)
print(sample_type_stats)

# Comparación por tratamiento
treatment_stats <- compare_means(
  subclones_per_1000_cells ~ treated,
  data = clonality_patients,
  method = "wilcox.test"
)
print(treatment_stats)

# Comparación por TME archetype
tme_stats <- compare_means(
  subclones_per_1000_cells ~ tme_archetype_group,
  data = clonality_patients,
  method = "wilcox.test"
)
print(tme_stats)
