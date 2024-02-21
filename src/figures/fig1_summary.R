library(tidyverse)

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
        text = element_text(family = "Arial", size = 12)
        )

ggsave(
    patients_lines_dough,
    filename = "results/figures/patients_cell_lines_dough.png", 
    dpi = 300,
    height = 7,
    width = 7
    )

## cancer type summary
primary_met_by_cancer <- clinical_annotation_samples %>%
    filter(study != "cell_lines_gabriella_kinker") %>%
    group_by(tumor_type, sample_type) %>%
    mutate(
        n.cells = n()
    ) %>%
    group_by(tumor_type, sample_type) %>%
    reframe(
        n.patients = length(unique(patient)),
        n.cells = n()
    ) %>%
    filter(
        sample_type == "p"
    ) %>%
    arrange(desc(n.patients))


