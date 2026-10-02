library(dplyr)
library(tidyverse)

setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta")
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")

# Load cell metadata
metadata <- read.table("cohort_statistics/tcca_metadata_h5ad.tsv", header = TRUE, sep = "\t")

# Load clonality information to select the same high-confident malignant and non-malignant cells than in scTherapy.
clonality <- read.table("cohort_statistics/full_clonality_table_lvl2.tsv", row.names = NULL)

metadata <- metadata %>%
    mutate(sample_id = paste0(study, "_", sample)) %>%
    left_join(select(clonality, c("original_barcode", "scevan_prediction")), by = c("cell" = "original_barcode"))

## !! Remember there are some cells from some subclones (a total of 60782 cells) that are predicted as malignant by SCEVAN but based on author annotation there are not malignant). We do not take into account the TC assignation for those cells, they were not taken into account during the drug prediction.
## !! Remember, on the other hand, there are 15 subclones (total of 698 cells), that are filtered out during scTherapy analysis because they raised errors due to low cell numbers, so they do not have predicted TC.

# Subset malignant cells with TC assignation (637,062 cells) and high-confident non-malignant cells from the same samples.
samples_with_TC <- metadata %>%
    filter(malignancy == "True" & !is.na(therapeutic_cluster)) %>%
    pull(sample_id) %>%
    unique()

metadata_hq <- metadata %>%
    filter(
        (malignancy == "True" & !is.na(therapeutic_cluster)) |
            (malignancy == "False" & sample_id %in% samples_with_TC & scevan_prediction == "normal" & scevan_subclone %in% c("non_tumor", ""))
    )

# Check number of non-malignant epithelial and TME cells per sample
epithelial_types <- c("Epithelial")

tme_types <- c(
    "Erythrocyte", "CD8+ T-cell", "Monocyte/Macrophage",
    "B-cell", "NK cell", "CD4+ T-cell", "Unconventional T-cells",
    "Stem", "Plasma cell", "Stromal cell", "Dendritic cell",
    "Plasmacytoid dendritic cell", "Mast", "Innate lymphoid cells",
    "Granulocyte", "Endothelial", "Regulatory T-cell",
    "Glial cell", "Neuron", "Unknown"
)

# Number of non-malignant epithelial and TME cells per sample
nonmal_per_sample <- metadata_hq %>%
    filter(malignancy == "False") %>%
    mutate(nonmal_group = case_when(
        cell_type_broad %in% epithelial_types ~ "NonMal_Epithelial",
        cell_type_broad %in% tme_types ~ "NonMal_TME",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(nonmal_group)) %>%
    group_by(sample_id, nonmal_group) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
        names_from  = nonmal_group,
        values_from = n,
        values_fill = 0
    )

# Number of malignant cells per subclone per sample
malignant_per_subclone <- metadata_hq %>%
    filter(malignancy == "True") %>%
    group_by(sample_id, scevan_subclone) %>%
    summarise(n_malignant = n(), .groups = "drop")

# Combine both tables
final_table <- malignant_per_subclone %>%
    left_join(nonmal_per_sample, by = "sample_id") %>%
    mutate(
        NonMal_Epithelial = replace_na(NonMal_Epithelial, 0),
        NonMal_TME = replace_na(NonMal_TME, 0)
    ) %>%
    select(
        sample_id,
        subclone = scevan_subclone,
        n_malignant,
        NonMal_Epithelial,
        NonMal_TME
    ) %>%
    arrange(sample_id, subclone)


write.table(
    final_table,
    "cohort_statistics/tcca_metadata_h5ad_nonmal_counts.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Barplot proportion ofNonMal_Epithelial and NonMal_TME per sample
liquid_tumors <- c("LAML", "MM", "CLL", "ALL")

data_plot <- nonmal_per_sample %>%
    left_join(
        metadata_hq %>%
            select(sample_id, tumor_type) %>%
            distinct(),
        by = "sample_id"
    ) %>%
    mutate(
        tumor_class = ifelse(tumor_type %in% liquid_tumors,
            "Liquid tumors",
            "Solid tumors"
        ),
        epi_category = case_when(
            NonMal_Epithelial == 0 ~ "0 cells",
            NonMal_Epithelial < 20 ~ "1-19 cells",
            NonMal_Epithelial < 100 ~ "20-99 cells",
            TRUE ~ "≥100 cells"
        ),
        tme_category = case_when(
            NonMal_TME == 0 ~ "0 cells",
            NonMal_TME < 20 ~ "1-19 cells",
            NonMal_TME < 100 ~ "20-99 cells",
            TRUE ~ "≥100 cells"
        )
    ) %>%
    pivot_longer(
        cols = c(epi_category, tme_category),
        names_to = "cell_group",
        values_to = "category"
    ) %>%
    mutate(
        category = factor(category,
            levels = c(
                "0 cells", "1-19 cells",
                "20-99 cells", "≥100 cells"
            )
        ),
        cell_group = recode(cell_group,
            "epi_category" = "Non-malignant Epithelial",
            "tme_category" = "Non-malignant TME"
        ),
        tumor_class = factor(tumor_class,
            levels = c("Solid tumors", "Liquid tumors")
        )
    ) %>%
    group_by(tumor_class, cell_group, category) %>%
    summarise(n_samples = n(), .groups = "drop")

barplot <- ggplot(
    data_plot,
    aes(x = category, y = n_samples, fill = cell_group)
) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    geom_text(aes(label = n_samples),
        position = position_dodge(width = 0.9),
        vjust = -0.5, size = 3
    ) +
    scale_fill_manual(values = c(
        "Non-malignant Epithelial" = "#F4A9A8", # rosa pastel
        "Non-malignant TME" = "#A5D6A7" # verde pastel
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~tumor_class, ncol = 2, scales = "free_y") +
    labs(
        x = "Cell count range",
        y = "Number of samples",
        fill = NULL,
        title = "Non-malignant reference cells used in scTherapy drug predictions"
    ) +
    theme_classic(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(size = 11, face = "bold")
    )

ggsave("cohort_statistics/nonmal_reference_cells.png", plot = barplot, width = 10, height = 5)

write.table(
    metadata_hq,
    "cohort_statistics/tcca_metadata_h5ad_sctherapy_cells.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
