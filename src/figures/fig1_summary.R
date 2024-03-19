library(RColorBrewer)
library(svglite)
library(tidyverse)
library(waffle)

source(file = "src/figures/TCCA_palette.R")

clinical_annotation_samples <- data.table::fread(
    "results/annotation/clinical_metadata_v4_clean.tsv"
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

## Doughnut plot for metagroups
primary_met_by_cancer$metagroup <- as.factor(primary_met_by_cancer$metagroup)
colnames(primary_met_by_cancer) <- c("Group", "Samples")

primary_met_by_cancer$fraction <- primary_met_by_cancer$Samples / sum(primary_met_by_cancer$Samples)

# Compute the cumulative percentages (top of each rectangle)
primary_met_by_cancer$ymax <- cumsum(primary_met_by_cancer$fraction)

# Compute the bottom of each rectangle
primary_met_by_cancer$ymin = c(0, head(primary_met_by_cancer$ymax, n = -1))

# Compute label position
primary_met_by_cancer$label_position <- (primary_met_by_cancer$ymax + primary_met_by_cancer$ymin) / 2

# Compute a good label
primary_met_by_cancer$label <- paste0(
    primary_met_by_cancer$Group, ": ",
    primary_met_by_cancer$Samples
)



## Doughnut plot for the same graph
patients_conditions <- ggplot(
    primary_met_by_cancer,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = Group)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 6, check_overlap = FALSE) +
    scale_fill_brewer(palette = "Set2", name = "") +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
    )


ggsave(
    plot = patients_conditions,
    filename = "results/figures/doughnut_patients_conditions_raw.svg",
    device = "svg"
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

## keep top 10 by grand total
primary_top <- primary_met_tops %>%
    group_by(tumor_type) %>%
    reframe(
        n.total = n.total
    ) %>%
    distinct(.keep_all = TRUE) %>%
    arrange(desc(n.total)) %>%
    head(10) %>%
    pull(tumor_type)

primary_met_tops <- primary_met_tops %>%
    filter(tumor_type %in% primary_top)

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
    ggtitle("Top 10 primary tumor types by number of patients") +
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
    filename = "results/figures/primary_tumors_barplot_raw.svg",
    device = "svg"
    )


## Additional doughnut plots

## Primary/metastasis
primary_metastatic_patients <- clinical_annotation_samples %>%
    filter(
        study != "cell_lines_gabriella_kinker"
    ) %>%
    group_by(sample_type) %>%
    summarise(
        n.samples = n()
    )

primary_metastatic_patients$sample_type <- as.factor(primary_metastatic_patients$sample_type)
levels(primary_metastatic_patients$sample_type) <- c("Metastasis", "Primary")
colnames(primary_metastatic_patients) <- c("Primary or metastasis", "Samples")

primary_metastatic_patients$fraction <- primary_metastatic_patients$Samples / sum(primary_metastatic_patients$Samples)

# Compute the cumulative percentages (top of each rectangle)
primary_metastatic_patients$ymax <- cumsum(primary_metastatic_patients$fraction)

# Compute the bottom of each rectangle
primary_metastatic_patients$ymin = c(0, head(primary_metastatic_patients$ymax, n = -1))

# Compute label position
primary_metastatic_patients$label_position <- (primary_metastatic_patients$ymax + primary_metastatic_patients$ymin) / 2

# Compute a good label
primary_metastatic_patients$label <- paste0(
    primary_metastatic_patients$`Primary or metastasis`, ": ",
    primary_metastatic_patients$Samples
)


# Make the plot
pm_dough <- ggplot(
    primary_metastatic_patients,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = `Primary or metastasis`)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 10) +
    scale_fill_manual(
        values = c(
            "Metastasis" = "#C10044",
            "Primary" ="#F0BFD0")
    ) +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
    )

ggsave(
    pm_dough,
    filename = "results/figures/pm_dough_raw.svg",
    width = 7,
    height = 7
)


## Male/Female
male_female_patients <- clinical_annotation_samples %>%
    filter(
        study != "cell_lines_gabriella_kinker"
    ) %>%
    group_by(sex) %>%
    summarise(
        n.samples = n()
    )

male_female_patients$sex <- as.factor(male_female_patients$sex)
levels(male_female_patients$sex) <- c("Unknown", "Female", "Male")
colnames(male_female_patients) <- c("Chromosomal sex", "Samples")

male_female_patients$fraction <- male_female_patients$Samples / sum(male_female_patients$Samples)

# Compute the cumulative percentages (top of each rectangle)
male_female_patients$ymax <- cumsum(male_female_patients$fraction)

# Compute the bottom of each rectangle
male_female_patients$ymin = c(0, head(male_female_patients$ymax, n = -1))

# Compute label position
male_female_patients$label_position <- (male_female_patients$ymax + male_female_patients$ymin) / 2

# Compute a good label
male_female_patients$label <- paste0(
    male_female_patients$`Chromosomal sex`, ": ",
    male_female_patients$Samples
)


# Make the plot
mf_dough <- ggplot(
    male_female_patients,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = `Chromosomal sex`)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 10) +
    scale_fill_manual(
        values = c(
            "Female" = "#1c8c8c",
            "Male" = "#ec5c44",
            "Unknown" = "gray50"
        )
    ) +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
    )

ggsave(
    mf_dough,
    filename = "results/figures/mf_dough_raw.svg",
    width = 7,
    height = 7
)


## Treated/untreated
treated_untreated_patients <- clinical_annotation_samples %>%
    filter(
        study != "cell_lines_gabriella_kinker"
    ) %>%
    group_by(treated) %>%
    summarise(
        n.samples = n()
    )

treated_untreated_patients$treated <- as.factor(treated_untreated_patients$treated)
levels(treated_untreated_patients$treated) <- c("Unknown", "Treatment naive", "Treated")
colnames(treated_untreated_patients) <- c("Treatment history", "Samples")

treated_untreated_patients$fraction <- treated_untreated_patients$Samples / sum(treated_untreated_patients$Samples)

# Compute the cumulative percentages (top of each rectangle)
treated_untreated_patients$ymax <- cumsum(treated_untreated_patients$fraction)

# Compute the bottom of each rectangle
treated_untreated_patients$ymin = c(0, head(treated_untreated_patients$ymax, n = -1))

# Compute label position
treated_untreated_patients$label_position <- (treated_untreated_patients$ymax + treated_untreated_patients$ymin) / 2

# Compute a good label
treated_untreated_patients$label <- paste0(
    treated_untreated_patients$`Treatment history`, ": ",
    treated_untreated_patients$Samples
)


# Make the plot
treatment_dough <- ggplot(
    treated_untreated_patients,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = `Treatment history`)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 10) +
    scale_fill_manual(
        values = c(
            "Treatment naive" = "#D18B6E",
            "Treated" = "#6ED1BC",
            "Unknown" = "gray50"
        )
    ) +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
    )

ggsave(
    treatment_dough,
    filename = "results/figures/treatment_dough_raw.svg",
    width = 7,
    height = 7
)

## adult/pediatric
pediatric_adult_samples <- clinical_annotation_samples %>%
    filter(
        study != "cell_lines_gabriella_kinker"
    ) %>%
    mutate(
        is_pediatric = age < 18,
    ) %>%
    group_by(is_pediatric) %>%
    summarise(
        n.samples = n()
    )

pediatric_adult_samples$is_pediatric <- as.factor(pediatric_adult_samples$is_pediatric)
levels(pediatric_adult_samples$is_pediatric) <- c("Adult", "Pediatric")
colnames(pediatric_adult_samples) <- c("Age group", "Samples")

pediatric_adult_samples$fraction <- pediatric_adult_samples$Samples / sum(pediatric_adult_samples$Samples)

# Compute the cumulative percentages (top of each rectangle)
pediatric_adult_samples$ymax <- cumsum(pediatric_adult_samples$fraction)

# Compute the bottom of each rectangle
pediatric_adult_samples$ymin = c(0, head(pediatric_adult_samples$ymax, n = -1))

# Compute label position
pediatric_adult_samples$label_position <- (pediatric_adult_samples$ymax + pediatric_adult_samples$ymin) / 2

# Compute a good label
pediatric_adult_samples$label <- paste0(
    pediatric_adult_samples$`Age group`, ": ",
    pediatric_adult_samples$Samples
)


# Make the plot
pediatric_adult_dough <- ggplot(
    pediatric_adult_samples,
    aes(
        ymax = ymax,
        ymin = ymin,
        xmax = 4,
        xmin = 3,
        fill = `Age group`)) +
    geom_rect() +
    geom_text(x = 1, aes(y = label_position, label = label), size = 10) +
    scale_fill_manual(
        values = c(
            "Adult" = "#7689DE",
            "Pediatric" = "#a9dce3",
            "Unknown" = "gray50"
        )
    ) +
    coord_polar(theta = "y") + 
    xlim(c(-1, 4)) +
    theme_void() +
    theme(
        legend.position = "none",
        text = element_text(family = "Arial", size = 24, color = "black")
    )

ggsave(
    pediatric_adult_dough,
    filename = "results/figures/pediatric_adult_raw.svg",
    width = 7,
    height = 7,
    device = "svg"
)