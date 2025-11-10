library(dplyr)
library(tidyverse)
library(entropy)
library(ggrepel)
library(htmlwidgets)
setwd("/storage/scratch01/shared/projects/bc-meta/")
source(file = "~/bc-meta/figures/TCCA_palette.R")

# Load the therapeutic similarity matrix
similarity_matrix <- readRDS("./single_cell/sctherapy/results/jaccard_matrix.rds")
colnames(similarity_matrix) <- gsub("[^a-zA-Z0-9]", "", colnames(similarity_matrix))
rownames(similarity_matrix) <- gsub("[^a-zA-Z0-9]", "", rownames(similarity_matrix))

## Compute genomic heterogeneity per sample
cnv_mat <- read.table("single_cell/cna_metadata/cnv_segments_clones_lvl2_cytobands.tsv")
cnv_mat <- as.matrix(cnv_mat)
colnames(cnv_mat) <- gsub("subclone", "", gsub("[^a-zA-Z0-9]", "", colnames(cnv_mat)))
cnv_mat <- apply(cnv_mat, c(1, 2), as.numeric)
cnv_mat[is.na(cnv_mat)] <- 2

# Compute frequency of subclones per sample
tcca_metadata <- read.table(
    "single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t", header = TRUE
)

subclone_freq <- tcca_metadata %>%
    mutate(scevan_subclone = gsub("[^a-zA-Z0-9]", "", scevan_subclone)) %>%
    filter(malignancy == "True" & scevan_subclone %in% colnames(similarity_matrix)) %>%
    group_by(study_sample, scevan_subclone) %>%
    summarise(n_cells = n()) %>%
    group_by(study_sample) %>%
    mutate(p = n_cells / sum(n_cells)) %>%
    ungroup()

samples <- unique(subclone_freq$study_sample)
entropy_sample <- setNames(numeric(length(samples)), samples)


prop_by_subclon <- lapply(colnames(cnv_mat), function(subclone) {
    cnv_subclone <- cnv_mat[, subclone]
    bins <- discretize(cnv_subclone, numBins = 5)
    prop <- bins / sum(bins)
    return(prop)
})

names(prop_by_subclon) <- colnames(cnv_mat)

for (s in samples) {
    subclones <- subclone_freq$scevan_subclone[subclone_freq$study_sample == s]
    fracs <- subclone_freq$p[subclone_freq$study_sample == s]
    mat <- do.call(rbind, prop_by_subclon[subclones])
    weighted_prop <- colSums(mat * fracs)
    entropy_sample[s] <- entropy(weighted_prop, unit = "log2")
}

# Genomic ITH
genomic_ith <- data.frame(
    study_sample = names(entropy_sample),
    genomic_ith = entropy_sample,
    row.names = NULL,
    stringsAsFactors = FALSE
)

## Compute therapeutic heterogeneity per sample (Rao's Quadratic Entropy, Rao Q)
dist_mat <- 1 - similarity_matrix

compute_rao <- function(p_vec, dist_mat) {
    subclones <- names(p_vec)
    dist_sub <- dist_mat[subclones, subclones, drop = FALSE]
    p <- matrix(p_vec, ncol = 1)
    as.numeric(t(p) %*% dist_sub %*% p)
}

therapeutic_ith <- subclone_freq %>%
    group_by(study_sample) %>%
    summarise(
        rao_q = {
            p_vec <- setNames(p, scevan_subclone)
            compute_rao(p_vec, dist_mat)
        }
    )


# Join genomic and therapeutic heterogeneity in a single table and add clinical variables
corr_genomic_therapeutic <- genomic_ith %>%
    select(study_sample, genomic_ith) %>%
    inner_join(therapeutic_ith, by = "study_sample")

# Add number of cells per sample
num_cells_sample <- subclone_freq %>%
    group_by(study_sample) %>%
    summarise(n_cells = sum(n_cells))
corr_genomic_therapeutic <- corr_genomic_therapeutic %>%
    inner_join(num_cells_sample, by = "study_sample")

# Add number of drugs per sample
sctherapy <- read.table(
    "single_cell/sctherapy/results/drug_response_subclone_sctherapy.tsv",
    sep = "\t", header = TRUE
)
num_drugs_sample <- sctherapy %>%
    mutate(study_sample = paste0(Study, "_", Sample)) %>%
    group_by(study_sample) %>%
    summarise(n_drugs = n_distinct(Drug.Name))
corr_genomic_therapeutic <- corr_genomic_therapeutic %>%
    inner_join(num_drugs_sample, by = "study_sample")


# Join clinical data to genomic and therapeutic ITH table
clinical_data <- tcca_metadata %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(
        study_sample,
        refined_tumor_type,
        sample_type,
        treated
    ) %>%
    distinct()
    
corr_genomic_therapeutic <- corr_genomic_therapeutic %>%
    left_join(clinical_data, by = "study_sample")

# Plot distribution of number of drugs per sample
plot <- ggplot(corr_genomic_therapeutic, aes(x = n_drugs)) +
    geom_histogram(binwidth = 5, fill = "#3182bd", color = "black", alpha = 0.7) +
    geom_vline(xintercept = c(20, 40), linetype = "dashed", color = "red") +
    labs(
        title = "Distribución del número de fármacos efectivos por muestra",
        x = "Número de fármacos (n_drugs)",
        y = "Frecuencia"
    ) +
    theme_bw(base_size = 14)
ggsave(
    "single_cell/genomic_ith/figures/num_drugs_distribution.png",
    plot = plot,
    width = 15,
    height = 6,
    dpi = 300
)

# Plot correlation between genomic and therapeutic heterogeneity
theme <- theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", face = "bold"),
        legend.title = element_text(color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

corr_plot_faceted <- corr_genomic_therapeutic %>%
    filter(sample_type %in% c("p", "m")) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("p", "m"),
            labels = c("Primary", "Metastasis")
        ),
        n_drugs_cat = cut(
            n_drugs,
            breaks = c(0, 30, 40, 50, 60, Inf),
            labels = c("<30", "30-40", "40-50", "50-60", ">60"),
            right = TRUE
        )
    ) %>%
    ggplot(aes(
        x = rao_q,
        y = genomic_ith,
        color = refined_tumor_type,
        size = n_drugs_cat
    )) +
    geom_point(alpha = 0.8) +
    facet_wrap(~sample_type, ncol = 2) +
    labs(
        x = "Therapeutic ITH",
        y = "Genomic ITH",
        color = "Cancer type",
        size = "Number of effective drugs"
    ) +
    scale_color_manual(values = tumor_type_colors) +
    scale_size_manual(
        values = c(2, 3, 4.5, 6, 8), # tamaños crecientes por categoría
        drop = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 9)) +
    theme_bw(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.background = element_blank(),
        legend.box.background = element_blank()
    )

ggsave(
    "single_cell/genomic_ith/figures/genomic_vs_therapeutic_ith_primary_vs_met.pdf",
    plot = corr_plot_faceted,
    width = 15,
    height = 6,
    dpi = 300
)


# Interactive version of the plot
library(ggplot2)
library(plotly)
library(dplyr)

corr_plot_faceted_interactive <- corr_genomic_therapeutic %>%
    filter(sample_type %in% c("p", "m")) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("p", "m"),
            labels = c("Primary", "Metastasis")
        ),
        n_drugs_cat = cut(
            n_drugs,
            breaks = c(0, 30, 40, 50, 60, Inf),
            labels = c("<30", "30-40", "40-50", "50-60", ">60"),
            right = TRUE
        )
    ) %>%
    ggplot(aes(
        x = rao_q,
        y = genomic_ith,
        color = refined_tumor_type,
        size = n_drugs_cat,
        text = paste0(
            "Sample: ", study_sample, "<br>",
            "Cancer type: ", refined_tumor_type, "<br>",
            "Genomic ITH: ", round(genomic_ith, 2), "<br>",
            "Therapeutic ITH: ", round(rao_q, 2), "<br>",
            "Effective drugs: ", n_drugs
        )
    )) +
    geom_point(alpha = 0.8) +
    facet_wrap(~sample_type, ncol = 2) +
    labs(
        x = "Therapeutic ITH",
        y = "Genomic ITH"
    ) +
    scale_color_manual(values = tumor_type_colors) +
    scale_size_manual(
        values = c(2, 3, 4.5, 6, 8),
        drop = FALSE
    ) +
    guides(color = "none", size = "none") + # Esto quita las leyendas
    theme_bw(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_blank(),
        legend.position = "none" # También puedes usar esto
    )

# Convertir a plotly
corr_plot_faceted_interactive <- ggplotly(corr_plot_faceted_interactive, tooltip = "text")

htmlwidgets::saveWidget(
    corr_plot_faceted_interactive,
    "single_cell/genomic_ith/figures/corr_plot_faceted_interactive.html",
    selfcontained = TRUE
)
