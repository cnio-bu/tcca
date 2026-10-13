library(dplyr)
library(tidyverse)

setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/therapeutic_analysis/")
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")

# Load drugs predicted for subclones in samples with enough non-malignant
# epithelial cells
drug_pred_epi <- read.table("therapeutic_epi_non_mal/full_table_drug_prediction_epithelial.tsv", sep = "\t", header = TRUE)

# Load drug prediction for all subclones
drug_pred <- read.table(
    "../cohort_statistics/drug_response_subclone_final.tsv",
    sep = "\t", header = TRUE
)

# Compute Jaccard Index and overlap coefficient
results <- lapply(subclones, function(s) {
    set_epi <- drug_pred_epi %>%
        filter(subclone == s) %>%
        pull(Drug_Name) %>%
        toupper()
    set_all <- drug_pred %>%
        filter(Subclone.Name == s) %>%
        pull(Drug.Name) %>%
        toupper()
    tumor_type <- drug_pred %>%
        filter(Subclone.Name == s) %>%
        pull(Refined.Tumor.Type) %>%
        unique() %>%
        sub(" .*", "", .)

    n_shared <- length(intersect(set_epi, set_all))
    n_union <- length(union(set_epi, set_all))
    n_epi <- length(set_epi)
    n_all <- length(set_all)
    n_min <- min(n_epi, n_all)

    data.frame(
        subclone = s,
        n_epi = n_epi,
        n_all = n_all,
        n_shared = n_shared,
        tumor_type = tumor_type,
        jaccard = ifelse(n_union > 0, n_shared / n_union, NA),
        overlap_coeff = ifelse(n_min > 0, n_shared / n_min, NA),
        recall = ifelse(n_epi > 0, n_shared / n_epi, NA)
    )
}) %>% bind_rows()

cat("\nJaccard index:\n")
cat("Median: ", round(median(results$jaccard, na.rm = TRUE), 3), "\n")
cat("Mean: ", round(mean(results$jaccard, na.rm = TRUE), 3), "\n")

cat("\nOverlap coefficient:\n")
cat("Median: ", round(median(results$overlap_coeff, na.rm = TRUE), 3), "\n")
cat("Mean: ", round(mean(results$overlap_coeff, na.rm = TRUE), 3), "\n")

cat("\nNumber of predicted drugs:\n")
cat("Epithelial median: ", round(median(results$n_epi, na.rm = TRUE), 1), "\n")
cat("All-non-malignant median: ", round(median(results$n_all, na.rm = TRUE), 1), "\n")
cat("Shared median: ", round(median(results$n_shared, na.rm = TRUE), 1), "\n")

cat("Recall (epithelial as reference):\n")
cat("Median: ", round(median(results$recall, na.rm = TRUE), 3), "\n")
cat("Mean: ", round(mean(results$recall, na.rm = TRUE), 3), "\n")

# Mean and median per cancer type
results %>%
    group_by(tumor_type) %>%
    summarise(
        n_subclones = n(),
        median_jaccard = round(median(jaccard, na.rm = TRUE), 3),
        median_recall = round(median(recall, na.rm = TRUE), 3)
    ) %>%
    arrange(desc(median_jaccard)) %>%
    print(n = 30)

# Plot the results to visualized the overlap metrics
results_long <- results %>%
    select(subclone, jaccard, recall) %>%
    pivot_longer(
        cols = c(jaccard, recall),
        names_to = "metric",
        values_to = "value"
    ) %>%
    mutate(metric = recode(metric,
        "jaccard" = "Jaccard index",
        "recall" = "Recall (non-malignant epithelial as reference)"
    ))

# Medians for vertical lines
medians <- results_long %>%
    group_by(metric) %>%
    summarise(med = median(value, na.rm = TRUE))

hist <- ggplot(results_long, aes(x = value, fill = metric)) +
    geom_histogram(bins = 40, alpha = 0.8, color = "white", linewidth = 0.2) +
    geom_vline(
        data = medians,
        aes(xintercept = med, color = metric),
        linetype = "dashed",
        linewidth = 0.8
    ) +
    geom_text(
        data = medians,
        aes(x = med, y = Inf, label = paste0("Median = ", round(med, 3)), color = metric),
        vjust = 1.5, hjust = -0.1, size = 3.5
    ) +
    facet_wrap(~metric, ncol = 1) +
    scale_fill_manual(values = c("#4393c3", "#d6604d")) +
    scale_color_manual(values = c("#2166ac", "#b2182b")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
        x = "Score",
        y = "Number of subclones",
        title = "Concordance between epithelial and all non-malignant references"
    ) +
    theme_classic(base_size = 11) +
    theme(
        legend.position  = "none",
        strip.text       = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = NA)
    )

ggsave(
    "therapeutic_epi_non_mal/concordance_epi_vs_all_nonmal.png",
    width = 5,
    height = 6,
    dpi = 300
)

# Concordance by tumor type
metadata_hq <- read.table(
    "../cohort_statistics/tcca_metadata_h5ad_sctherapy_cells.tsv",
    header = TRUE,
    sep = "\t"
) %>%
    mutate(study_sample = paste0(study, "_", sample))

samples_with_epithelial <- metadata_hq %>%
    filter(malignancy == "False" & cell_type_broad == "Epithelial") %>%
    group_by(study_sample) %>%
    summarise(n_epithelial = n()) %>%
    filter(n_epithelial >= 20) %>%
    pull(study_sample)

cell_counts_per_sample <- metadata_hq %>%
    filter(
        study_sample %in% samples_with_epithelial &
            (malignancy == "True" | (malignancy == "False" & cell_type_broad == "Epithelial"))
    ) %>%
    group_by(study_sample, tumor_type) %>%
    summarise(
        n_malignant = sum(malignancy == "True"),
        n_epithelial = sum(malignancy == "False"),
        .groups = "drop"
    )

results_with_sample <- results %>%
    mutate(study_sample = sub("\\.", "_", sub("\\.[^.]*$", "", subclone))) %>% # extraer muestra del subclón
    group_by(study_sample, tumor_type) %>%
    summarise(
        median_jaccard = median(jaccard, na.rm = TRUE),
        n_subclones = n(),
        .groups = "drop"
    ) %>%
    left_join(cell_counts_per_sample, by = "study_sample") %>%
    mutate(
        tumor_type = coalesce(tumor_type.x, tumor_type.y)
    ) %>%
    select(-tumor_type.x, -tumor_type.y) %>%
    arrange(desc(median_jaccard))

barplot <- ggplot(
    results_with_sample,
    aes(
        x = reorder(study_sample, -median_jaccard),
        y = median_jaccard,
        fill = tumor_type
    )
) +
    geom_col(width = 0.7) +
    geom_point(
        aes(y = n_epithelial / max(n_epithelial, na.rm = TRUE)),
        color = "black", size = 0.8, alpha = 0.6
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        sec.axis = sec_axis(
            ~ . * max(results_with_sample$n_epithelial, na.rm = TRUE),
            name = "N non-malignant epithelial cells"
        )
    ) +
    geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.5
    ) +
    labs(
        x = "Sample",
        y = "Median Jaccard index",
        fill = "Tumor type",
        title = "Concordance between reference strategies per sample"
    ) +
    scale_fill_manual(values = tumor_type_colors) +
    theme_classic(base_size = 10) +
    theme(
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "right",
        plot.title      = element_text(hjust = 0.5, size = 11)
    )

ggsave(
    "therapeutic_epi_non_mal/concordance_by_tumor_type.png",
    plot = barplot,
    width = 7,
    height = 5,
    dpi = 300
)
