library(ComplexHeatmap)
library(limma)
library(tidyverse)

gsva_mat <- read.table(
    file = "results/functional/pseudo_bulk_gsva.tsv",
    check.names = FALSE
    )

gsva_mat <- as.matrix(gsva_mat)

## load metacommunity data from patients
metacommunity_data <- read_tsv(
    "results/modules/annotated/metacom_proportions_primary_wide.tsv"
    )  %>%
    mutate(
        sample_functional_mat = paste(sample, study, sep = "_"),
        sample_functional_mat = gsub(pattern = "-", replacement = "_", x = sample_functional_mat),
        sample_functional_mat = gsub(pattern = "\\.", replacement = "_", x = sample_functional_mat)
    ) %>%
    filter(
        sample_type == "p",
        treated == FALSE,
        study != "cell_lines_gabriella_kinker"
    )

metacommunity_data$best_metacom <- as.factor(metacommunity_data$best_metacom)
levels(metacommunity_data$best_metacom) <- c(
    "Metacommunity 1",
    "Metacommunity 2",
    "Metacommunity 3",
    "Metacommunity 4",
    "Metacommunity 5",
    "Metacommunity 6"
    
)

## Make sure design and matrix are in the same order
gsva_mat <- gsva_mat[, metacommunity_data$sample_functional_mat]


## Test MPs first
gsva_mat <- gsva_mat[grepl(pattern = "^MP", x = rownames(gsva_mat)), ]

## Make a 0 intercept model as it makes sense here
design <- model.matrix(~0 + best_metacom, data = metacommunity_data)
rownames(design) <- metacommunity_data$sample_functional_mat

colnames(design) <- paste0("metcom_", c(1:6))

## Initial lm fit
fit <- lmFit(object = gsva_mat, design = design)

## get all contrasts
all_contrasts <- makeContrasts(
    metcom_2 - metcom_1,
    metcom_3 - metcom_1,
    metcom_4 - metcom_1,
    metcom_5 - metcom_1,
    metcom_6 - metcom_1,
    metcom_3 - metcom_2,
    metcom_4 - metcom_2,
    metcom_5 - metcom_2,
    metcom_6 - metcom_2,
    metcom_4 - metcom_3,
    metcom_5 - metcom_3,
    metcom_6 - metcom_3,
    metcom_5 - metcom_4,
    metcom_6 - metcom_4,
    metcom_6 - metcom_5,
    levels = design
    )

fit2 <- contrasts.fit(fit, contrasts = all_contrasts)
fit2 <- eBayes(fit2)

## get all res sim.
all_tests <- decideTests(
    fit2,
    adjust.method = "BH",
    method = "global",
    p.value = 0.05
    )

## Keep significant assocs
tests <- as.data.frame(all_tests@.Data)
significant_tests <- rowSums(abs(tests)) > 0
keep <- significant_tests[significant_tests]

top_annotation <- HeatmapAnnotation(
    "Best metacommunity" = metacommunity_data$best_metacom
)

bulk_heat <- ComplexHeatmap::Heatmap(
    matrix = gsva_mat[names(keep), ],
    show_column_names = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = TRUE,
    cluster_rows = FALSE,
    cluster_row_slices = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    column_split = metacommunity_data$best_metacom,
    row_split = c(1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,5,5,6,7,8,9,9,10,11,11,11,11,11,12,12,12,13,13,13,13,13,14,14,14,14),
    top_annotation = top_annotation,
    row_names_gp = gpar(fontsize = 7),
    row_title = c(
        "Cell cycle",
        "Stress or hypoxia",
        "Protein reg.",
        "Mesenchymal",
        "Interferon and MHC II",
        "Senescence",
        "Oncogenic",
        "Respiration",
        "Secreted",
        "Cilia",
        "Lineage Neural",
        "Lineage Other",
        "Lineage Haemat.",
        "Unassigned"
        ),
    row_title_rot = 0,
    row_title_gp = (gpar(fontsize = 7)),
    column_title_gp = (gpar(fontsize = 7)),
    heatmap_legend_param = list(
        title = "Enrichment score",
        direction = "horizontal",
        title_position = "lefttop"
    )
)
