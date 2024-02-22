library(RColorBrewer)
library(svglite)
library(tidyverse)
library(waffle)

source(file = "src/figures/TCCA_palette.R")

clinical_annotation_samples <- data.table::fread(
    "results/annotation/clinical_metadata_v4_clean.tsv"
)

clinical_annotation_cells <- data.table::fread(
    "results/annotation/beyondcell_with_therapeutic_clusters.tsv"
    )


## build a doughnut chart for cll/patients
samples_cells <- clinical_annotation_samples %>%
    mutate(
        is_cell_line = study == "cell_lines_gabriella_kinker"  
    ) %>%
    group_by(is_cell_line) %>%
    summarise(
        n.samples = n()
    )

samples_cells$is_cell_line <- as.factor(samples_cells$is_cell_line)
levels(samples_cells$is_cell_line) <- c("Patients", "Cancer cell lines")
colnames(samples_cells) <- c("Model", "Samples")

samples_cells$fraction <- samples_cells$Samples / sum(samples_cells$Samples)

# Compute the cumulative percentages (top of each rectangle)
samples_cells$ymax <- cumsum(samples_cells$fraction)

# Compute the bottom of each rectangle
samples_cells$ymin = c(0, head(samples_cells$ymax, n = -1))

# Compute label position
samples_cells$label_position <- (samples_cells$ymax + samples_cells$ymin) / 2

# Compute a good label
samples_cells$label <- paste0(
    samples_cells$Model, ": ",
    samples_cells$Samples
    )


# Make the plot
patients_lines_dough <- ggplot(
    samples_cells,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = Model)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 6) +
    scale_fill_manual(
        values = c(
            "Patients" = "#F6Bd60",
            "Cancer cell lines" = "#706695")
        ) +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
        )

ggsave(
    patients_lines_dough,
    filename = "results/figures/patients_cell_lines_dough.png",
    width = 7,
    height = 7
    )

## test
pdf(
    file = "results/figures/patients_cell_lines_dough.pdf",
    width = 7,
    height = 7
    )

plot(patients_lines_dough)
dev.off()

## cancer type summary
primary_met_by_cancer <- clinical_annotation_samples %>%
    filter(study != "cell_lines_gabriella_kinker") %>%
    group_by(sample_type, treated) %>%
    reframe(
        n.patients = n()
    ) %>% 
    arrange(desc(n.patients)) %>%
    mutate(
        metagroup = paste0(sample_type, "_", treated)
    ) %>%
    select(
        metagroup, n.patients
    ) %>%
    mutate(
        metagroup = as_factor(metagroup)
    )

levels(primary_met_by_cancer$metagroup) <- c(
    "Primary tumor, untreated",
    "Primary tumor, treated", 
    "Primary tumor, unknown treatment",
    "Metastasis, untreated",
    "Metastasis, treated",
    "Metastasis, unknown treatment"
    )


metagroup_waffle <- ggplot(
    data = primary_met_by_cancer,
    aes(fill = metagroup, values = n.patients)
    ) +
    geom_waffle(n_rows = 10, size = 0.33, color = "white") +
    scale_fill_brewer(palette = "Set2", name = "") +
    coord_equal() +
    theme_void() +
    theme(
        text = element_text(family = "Arial", size = 11)
    )


ggsave(
    filename = "results/figures/waffleplot_patients.png",
    plot = metagroup_waffle,
    width = 12,
    height = 4,
    dpi = 300
    )


## Primary / metastatic proportions
primary_met_tops <- data.table::fread("results/annotation/clinical_metadata_v4_clean.tsv") %>%
    filter(
        study != "cell_lines_gabriella_kinker",
        sample_type == "p"
    ) %>%
    group_by(tumor_type, treated) %>%
    summarise(
        n.samples = n()
    ) %>%
    ungroup() %>%
    group_by(tumor_type) %>%
    mutate(
        n.total = sum(n.samples)
    ) %>%
    mutate(
        treated = case_when(
            treated == "" ~ "unknown",
            TRUE ~ treated
        )
    ) %>%
    arrange(desc(n.total))

primary_met_tops$tumor_type <- fct_reorder(
    primary_met_tops$tumor_type,
    primary_met_tops$n.total,
    .desc = FALSE
)


primary_tumors_barplot <- ggplot(
    data = primary_met_tops,
    aes(x = n.samples,
        y = tumor_type,
        fill = treated,
        label = n.samples
        )
    ) +
    geom_bar(stat = "identity") +
    geom_text(position = position_stack(vjust = 0.5)) +
    ggtitle("Top 10 cancer types by number of patients") +
    scale_fill_manual(
        name = "Treatment information",
        aesthetics = "fill",
        values = c(
            "t" = treatment_colors[["Treated"]],
            "f" = treatment_colors[["Untreated"]],
            "unknown" = "grey50"
        ),
        labels = c("Treated", "Untreated", "Unknown")
        ) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(family = "Arial")
    )

ggsave(
    plot = primary_tumors_barplot,
    filename = "results/figures/primary_tumors_barplot.png",
    height = 7,
    width = 7,
    dpi = 100
    )

mets_top <- data.table::fread("results/annotation/clinical_metadata_v4_clean.tsv") %>%
    filter(
        study != "cell_lines_gabriella_kinker",
        sample_type == "m"
    ) %>%
    group_by(refined_tumor_site, treated) %>%
    summarise(
        n.samples = n()
    ) %>%
    ungroup() %>%
    group_by(refined_tumor_site) %>%
    mutate(
        n.total = sum(n.samples)
    ) %>%
    mutate(
        treated = case_when(
            treated == "" ~ "unknown",
            TRUE ~ treated
        )
    ) %>%
    arrange(desc(n.total))

mets_top$refined_tumor_site <- fct_reorder(
    mets_top$refined_tumor_site,
    mets_top$n.total,
    .desc = FALSE
)


met_tumors_barplot <- ggplot(
    data = mets_top,
    aes(x = n.samples,
        y = refined_tumor_site,
        fill = treated,
        label = n.samples
    )
) +
    geom_bar(stat = "identity") +
    geom_text(position = position_stack(vjust = 0.5)) +
    ggtitle("Top 10 cancer types by number of patients") +
    scale_fill_manual(
        name = "Treatment information",
        aesthetics = "fill",
        values = c(
            "t" = treatment_colors[["Treated"]],
            "f" = treatment_colors[["Untreated"]],
            "unknown" = "grey50"
        ),
        labels = c("Treated", "Untreated", "Unknown")
    ) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(family = "Arial")
    )
