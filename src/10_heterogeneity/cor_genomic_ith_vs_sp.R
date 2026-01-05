library(dplyr)
library(tidyverse)
library(openxlsx)
library(purrr)
library(ggplot2)
library(ggpubr)

setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load genomic heterogeneity
sample_ith <- read.table(
    "single_cell/heterogeneity/ith_scores_per_sample.tsv",
    header = TRUE,
    sep = "\t"
)

sample_metadata <- read.table(
    "single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t",
    header = TRUE
) %>%
    select(study, sample, patient, refined_tumor_type, sample_type, tme_archetype_group) %>%
    distinct() %>%
    mutate(study_sample = gsub("[_. \\-]", "", paste0(study, "_", sample)))

genomic_ith <- sample_ith %>%
    mutate(match_names = gsub("[_. \\-]", "", study_sample)) %>%
    select(match_names, genomic_ith) %>%
    left_join(
        select(
            sample_metadata,
            c("study_sample", "patient", "refined_tumor_type", "sample_type")
        ),
        by = c("match_names" = "study_sample")
    ) %>%
    filter(patient != "ccl") # keep only patient samples

# Load sp per sample
full_sp_mat <- read.table(
    "beyondcell_immuno/full_sp_mat_filtered.tsv",
    sep = "\t",
    header = TRUE
)
colnames(full_sp_mat) <- gsub("[_. \\-]", "", colnames(full_sp_mat))
# Keep samples with genomic info
common_samples <- intersect(colnames(full_sp_mat), genomic_ith$match_names)
full_sp_mat <- full_sp_mat[, common_samples]
full_sp_mat <- t(full_sp_mat)
genomic_ith <- genomic_ith %>%
    filter(match_names %in% rownames(full_sp_mat))

# Drug names
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
    distinct(IDs, .keep_all = TRUE) %>%
    select(IDs, preferred.drug.names, collapsed.MoAs) %>%
    mutate(
        collapsed.MoAs = case_when(
            preferred.drug.names == "VANDETANIB" ~ "VEGFR inhibitor",
            preferred.drug.names == "DASATINIB" ~ "Kinase inhibitor",
            preferred.drug.names == "RIGOSERTIB" ~ "Other",
            preferred.drug.names == "SORAFENIB" ~ "Kinase inhibitor",
            TRUE ~ collapsed.MoAs
        )
    ) %>%
    as.data.frame()


## Compute pancancer correlation
# Scores of genomic heterogeneity
genomic_ith_scores <- "genomic_ith"

# Crear un workbook vacío
wb <- createWorkbook()

for (score in genomic_ith_scores) {
    ith_score <- genomic_ith[[score]]

    # Correlation of genomic ith with all drugs
    cor_results <- lapply(colnames(full_sp_mat), function(drug) {
        res <- cor.test(
            full_sp_mat[, drug],
            ith_score,
            method = "spearman",
            use = "pairwise.complete.obs"
        )
        tibble(
            drug = drug,
            cor = res$estimate,
            p.value = res$p.value
        )
    }) %>%
        bind_rows() %>%
        arrange(desc(cor))

    cor_results <- cor_results %>%
        left_join(drugs, by = c("drug" = "IDs"))

    # Añadir hoja al workbook
    addWorksheet(wb, sheetName = score)
    writeData(wb, sheet = score, cor_results)
}

# Guardar el archivo Excel
saveWorkbook(
    wb,
    "single_cell/heterogeneity/genomic_ith/genomic_heterogeneity_correlations.xlsx",
    overwrite = TRUE
)




# List of genomic heterogeneity scores
genomic_ith_scores <- "genomic_ith"

# List of cancer types
cancer_types <- unique(genomic_ith$refined_tumor_type)

for (score in genomic_ith_scores) {
    # Create empty workbook for this heterogeneity score
    wb <- createWorkbook()

    for (cancer in cancer_types) {
        # Filter samples for this cancer type
        idx <- genomic_ith$refined_tumor_type == cancer
        n_samples <- sum(idx)

        # Skip cancer types with fewer than 3 samples
        if (n_samples < 3) {
            message("Skipping cancer type ", cancer, " (", n_samples, " samples)")
            next
        }

        ith_score <- genomic_ith[[score]][idx]

        # Subset the SP matrix for these samples
        sp_mat_sub <- full_sp_mat[idx, , drop = FALSE]

        # Compute correlation of the heterogeneity score with all drugs
        cor_results <- lapply(colnames(sp_mat_sub), function(drug) {
            # Skip drugs with fewer than 3 non-NA observations
            if (sum(!is.na(sp_mat_sub[, drug])) < 3) {
                return(NULL)
            }

            res <- cor.test(
                sp_mat_sub[, drug],
                ith_score,
                method = "spearman",
                use = "pairwise.complete.obs"
            )
            tibble(
                drug = drug,
                cor = res$estimate,
                p.value = res$p.value
            )
        }) %>%
            bind_rows() %>%
            arrange(desc(cor))

        # Add drug information (MoA and preferred name)
        cor_results <- cor_results %>%
            left_join(drugs, by = c("drug" = "IDs"))

        # Add a worksheet for this cancer type
        addWorksheet(wb, sheetName = cancer)
        writeData(wb, sheet = cancer, cor_results)
    }

    # Save the workbook for this heterogeneity score
    saveWorkbook(
        wb,
        paste0("single_cell/heterogeneity/genomic_ith/correlations_", score, ".xlsx"),
        overwrite = TRUE
    )
}

# Compute correlations for each cancer type separately for primary and metastasis
# samples
# Select only samples with primary/metastasis annotation
genomic_ith <- genomic_ith %>% filter(sample_type %in% c("p", "m"))
genomic_ith_scores <- "genomic_ith"
full_sp_mat <- full_sp_mat[genomic_ith$match_names, ]
sample_types <- c("p", "m")
cancer_types <- unique(genomic_ith$refined_tumor_type)

for (score in genomic_ith_scores) {
    for (stype in sample_types) {
        # Filter for sample type
        idx_type <- genomic_ith$sample_type == stype
        if (sum(idx_type) < 3) next # skip if too few samples

        wb <- createWorkbook()

        for (cancer in cancer_types) {
            idx <- idx_type & genomic_ith$refined_tumor_type == cancer
            n_samples <- sum(idx)
            if (n_samples < 3) next # skip small cancer types

            ith_score <- genomic_ith[[score]][idx]
            sp_mat_sub <- full_sp_mat[idx, , drop = FALSE]

            cor_results <- lapply(colnames(sp_mat_sub), function(drug) {
                if (sum(!is.na(sp_mat_sub[, drug])) < 3) {
                    return(NULL)
                }
                res <- cor.test(
                    sp_mat_sub[, drug],
                    ith_score,
                    method = "spearman",
                    use = "pairwise.complete.obs"
                )
                tibble(
                    drug = drug,
                    cor = res$estimate,
                    p.value = res$p.value
                )
            }) %>%
                bind_rows() %>%
                arrange(desc(cor))

            # Add drug info if available
            if (exists("drugs")) {
                cor_results <- cor_results %>% left_join(drugs, by = c("drug" = "IDs"))
            }

            # Add worksheet
            addWorksheet(wb, sheetName = cancer)
            writeData(wb, sheet = cancer, cor_results)
        }

        # Save workbook
        saveWorkbook(
            wb,
            paste0("single_cell/heterogeneity/genomic_ith/correlations_", score, "_", stype, ".xlsx"),
            overwrite = TRUE
        )
    }
}


## Reformat results
path <- "single_cell/heterogeneity/genomic_ith/"

# Read all files with correlations
files <- list.files(path, pattern = "correlations_.*_(p|m)\\.xlsx$", full.names = TRUE)

results_df <- map_dfr(files, function(f) {
    score <- str_extract(basename(f), "(?<=correlations_)[^_]+")
    stype <- str_extract(basename(f), "(?<=_)[pm](?=\\.xlsx)")

    # Leer todas las hojas (cada hoja = un tipo de cáncer)
    sheets <- getSheetNames(f)

    map_dfr(sheets, function(s) {
        read.xlsx(f, sheet = s) %>%
            mutate(
                cancer_type = s,
                score = score,
                sample_type = stype
            )
    })
})

# Plot heatmap of correlations for primary samples
library(ComplexHeatmap)
source("~/bc-meta/figures/TCCA_palette.R")

corr_primary <- results_df %>%
    filter(sample_type == "p") %>%
    filter(score == "genomic") %>%
    select(cancer_type, drug, cor) %>%
    pivot_wider(names_from = cancer_type, values_from = cor, values_fill = NA) %>%
    column_to_rownames(var = "drug") %>%
    as.data.frame()

# Top annotation: cancer type
cancer_type <- as.data.frame(colnames(corr_primary))
colnames(cancer_type) <- "Cancer type"
top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = cancer_type,
    which = "column",
    col = list(
        "Cancer type" = tumor_type_colors
    ),
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = FALSE
)

drugs <- drugs[drugs$IDs %in% rownames(corr_primary), ]
rownames(drugs) <- drugs$IDs

drugs$collapsed.MoAs <- ifelse(
    drugs$collapsed.MoAs %in% names(MoAs_colors),
    drugs$collapsed.MoAs,
    "Other"
)
# Add mechanism of action to immunotherapies
drugs <- drugs[rownames(corr_primary), c("IDs", "preferred.drug.names", "collapsed.MoAs")]
rownames(drugs) <- rownames(corr_primary)
drugs$preferred.drug.names[is.na(drugs$preferred.drug.names)] <- rownames(corr_primary)[is.na(drugs$preferred.drug.names)]
MoAs <- drugs[, c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"
rownames(MoAs) <- rownames(corr_primary)

# Add mechanism of action to immunotherapies
MoAs$`Mechanism of action`[is.na(MoAs$`Mechanism of action`)] <- "Immunotherapy"
# Order drugs based on MoA
# ord <- order(MoAs[rownames(corr_primary), "Mechanism of action"])
# ordered_rownames <- rownames(corr_primary)[ord]
# ordered_drug_names <- drugs[ordered_rownames, "preferred.drug.names"]
# ordered_drug_names[is.na(ordered_drug_names)] <- ordered_rownames[is.na(ordered_drug_names)]

MoAs <- MoAs[rownames(corr_primary), , drop = FALSE]
moa_pals <- list(
    "Mechanism of action" = MoAs_colors
)

right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = MoAs,
    which = "row",
    col = moa_pals,
    show_annotation_name = FALSE,
    show_legend = TRUE
)


png(
    file = "single_cell/heterogeneity/genomic_ith/figures/heatmap_corr_primary.png",
    res = 500,
    width = 16,
    height = 22,
    units = "in"
)

mat <- as.matrix(corr_primary)
heat <- ComplexHeatmap::Heatmap(
    mat = mat,
    # col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    cluster_rows = TRUE,
    cluster_row_slices = TRUE,
    show_row_dend = TRUE,
    # row_split = MoAs$`Mechanism of action`,
    # column_order = rownames(sample_annot_df)[order(sample_annot_df$`Cancer type`)],
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    # column_split = cancer_type$`Cancer type`,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    # row_title = unique(sort(MoAs[rownames(full_sp_mat), "Mechanism of action"])),
    row_title = NULL,
    # row_order = ordered_rownames,
    row_labels = drugs[rownames(corr_primary), "preferred.drug.names"],
    show_row_names = TRUE,
    column_names_side = "top",
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 2),
    column_title = NULL,
    heatmap_legend_param = list(title = "Spearman ρ"),
    heatmap_width = unit(12, "in"),
    heatmap_height = unit(19, "in")
)

ht_opt(
    "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
    "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top")
dev.off()


# Select top correlations (positive and negative) per cancer type
top_corrs <- results_df %>%
    filter(
        score == "WSE",
        sample_type == "p",
        p.value <= 0.05,
        abs(cor) > 0.5
    ) %>%
    group_by(cancer_type) %>%
    slice_max(order_by = cor, n = 5, with_ties = FALSE) %>%
    bind_rows(
        results_df %>%
            filter(
                score == "WSE",
                sample_type == "p",
                p.value < 0.05,
                abs(cor) > 0.5
            ) %>%
            group_by(cancer_type) %>%
            slice_min(order_by = cor, n = 5, with_ties = FALSE)
    ) %>%
    arrange(cancer_type, desc(cor))

top_corrs[is.na(top_corrs$collapsed.MoAs), "collapsed.MoAs"] <- "Immunotherapy"
top_corrs[is.na(top_corrs$preferred.drug.names), "preferred.drug.names"] <- top_corrs[is.na(top_corrs$preferred.drug.names), "drug"]

# Add the number of samples in each correlation
n_samples_cancertype <- genomic_ith %>%
    filter(sample_type == "m") %>%
    group_by(refined_tumor_type) %>%
    summarise(n_samples = n())
top_corrs <- top_corrs %>%
    left_join(n_samples_cancertype, by = c("cancer_type" = "refined_tumor_type")) %>%
    as.data.frame()

library(ggplot2)

# Dotplot of top correlations
top_corrs <- top_corrs %>%
    filter(n_samples >= 10)
dotplot <- ggplot(top_corrs, aes(
    x = cancer_type,
    y = preferred.drug.names,
    color = cor,
    size = -log10(p.value)
)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        name = "Spearman ρ"
    ) +
    scale_size_continuous(name = expression(-log[10](p))) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        legend.position = "right"
    ) +
    labs(
        title = "Strong correlations between genomic ITH and Switch Point (Metastasis samples)",
        x = "Cancer type",
        y = "Drug (preferred name)"
    ) +
    coord_flip()

ggsave(
    file = "single_cell/heterogeneity/genomic_ith/figures/dotplot_metastasis_10samples.png",
    plot = dotplot,
    width = 10,
    height = 4
)


# Plot correlation for trastuzumab response in primary patient samples of breast cancer
brca_samples <- sample_metadata %>%
    filter(refined_tumor_type == "BRCA" & sample_type == "p" & patient != "ccl") %>%
    pull(study_sample)
sp <- as.data.frame(full_sp_mat[intersect(brca_samples, rownames(full_sp_mat)), "LIU_TRASTUZUMAB_RESPONSE"])
colnames(sp) <- "LIU_TRASTUZUMAB_RESPONSE"
genomic_ith_brca <- genomic_ith %>%
    filter(match_names %in% brca_samples) %>%
    select(match_names, genomic_ith) %>%
    column_to_rownames(var = "match_names")
cor_df <- sp %>%
    rownames_to_column(var = "match_names") %>%
    left_join(genomic_ith_brca %>% rownames_to_column(var = "match_names"), by = "match_names")

cor_plot <- ggscatter(
    cor_df,
    x = "genomic_ith",
    y = "LIU_TRASTUZUMAB_RESPONSE", # <- must match exactly the column name
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman",
    xlab = "WSE (Genomic Heterogeneity)",
    ylab = "LIU_TRASTUZUMAB_RESPONSE"
)
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/correlation_brca_trastuzumab.png",
    plot = cor_plot,
    width = 6, height = 6
)
cor_df <- cor_df %>%
    left_join(select(sample_metadata, c("study_sample", "tme_archetype_group")),
        by = c("match_names" = "study_sample")
    ) %>%
    mutate(predicted_response = case_when(
        LIU_TRASTUZUMAB_RESPONSE <= 0.3 ~ "Responder",
        LIU_TRASTUZUMAB_RESPONSE >= 0.7 ~ "Non-responder",
        TRUE ~ "Intermediate"
    )) %>%
    filter(predicted_response != "Intermediate")

barplot_tme <- ggplot(
    cor_df,
    aes(x = predicted_response, fill = tme_archetype_group)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = tme_group_colors) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_bw()
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/brca_trastuzumab_response_vs_tme.png",
    plot = barplot_tme,
    width = 6, height = 6
)




# Plot correlation for decitabine response in primary patient samples of pancreatic cancer
paad_samples <- sample_metadata %>%
    filter(refined_tumor_type == "PAAD" & sample_type == "p" & patient != "ccl") %>%
    pull(study_sample)
sp <- as.data.frame(full_sp_mat[intersect(paad_samples, rownames(full_sp_mat)), "LIU_TRASTUZUMAB_RESPONSE"])
colnames(sp) <- "LIU_TRASTUZUMAB_RESPONSE"
genomic_ith_brca <- genomic_ith %>%
    filter(match_names %in% brca_samples) %>%
    select(match_names, genomic_ith) %>%
    column_to_rownames(var = "match_names")
cor_df <- sp %>%
    rownames_to_column(var = "match_names") %>%
    left_join(genomic_ith_brca %>% rownames_to_column(var = "match_names"), by = "match_names")

cor_plot <- ggscatter(
    cor_df,
    x = "genomic_ith",
    y = "LIU_TRASTUZUMAB_RESPONSE", # <- must match exactly the column name
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman",
    xlab = "WSE (Genomic Heterogeneity)",
    ylab = "LIU_TRASTUZUMAB_RESPONSE"
)
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/correlation_brca_trastuzumab.png",
    plot = cor_plot,
    width = 6, height = 6
)

# Plot correlation for decitabine response in primary patient samples of pancreatic cancer
paad_samples <- sample_metadata %>%
    filter(refined_tumor_type == "PAAD" & sample_type == "p" & patient != "ccl") %>%
    pull(study_sample)
sp <- as.data.frame(full_sp_mat[intersect(paad_samples, rownames(full_sp_mat)), "LIU_TRASTUZUMAB_RESPONSE"])
colnames(sp) <- "LIU_TRASTUZUMAB_RESPONSE"
genomic_ith_brca <- genomic_ith %>%
    filter(match_names %in% brca_samples) %>%
    select(match_names, genomic_ith) %>%
    column_to_rownames(var = "match_names")
cor_df <- sp %>%
    rownames_to_column(var = "match_names") %>%
    left_join(genomic_ith_brca %>% rownames_to_column(var = "match_names"), by = "match_names")

cor_plot <- ggscatter(
    cor_df,
    x = "genomic_ith",
    y = "LIU_TRASTUZUMAB_RESPONSE", # <- must match exactly the column name
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman",
    xlab = "WSE (Genomic Heterogeneity)",
    ylab = "LIU_TRASTUZUMAB_RESPONSE"
)
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/correlation_brca_trastuzumab.png",
    plot = cor_plot,
    width = 6, height = 6
)

# Plot correlation for decitabine response in primary patient samples of pancreatic cancer
paad_samples <- sample_metadata %>%
    filter(refined_tumor_type == "PAAD" & sample_type == "p" & patient != "ccl") %>%
    pull(study_sample)
sp <- as.data.frame(full_sp_mat[intersect(paad_samples, rownames(full_sp_mat)), "sig-21342"])
colnames(sp) <- "sig.21342"
genomic_ith_brca <- genomic_ith %>%
    filter(match_names %in% brca_samples) %>%
    select(match_names, genomic_ith) %>%
    column_to_rownames(var = "match_names")
cor_df <- sp %>%
    rownames_to_column(var = "match_names") %>%
    left_join(genomic_ith_brca %>% rownames_to_column(var = "match_names"), by = "match_names")

cor_plot <- ggscatter(
    cor_df,
    x = "genomic_ith",
    y = "sig.21342", # <- must match exactly the column name
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman",
    xlab = "WSE (Genomic Heterogeneity)",
    ylab = "DECITABINE"
)
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/correlation_paad_decitabine.png",
    plot = cor_plot,
    width = 6, height = 6
)

# Plot correlation for decitabine response in primary patient samples of pancreatic cancer
brca_samples <- sample_metadata %>%
    filter(refined_tumor_type == "BRCA" & sample_type == "m" & patient != "ccl") %>%
    pull(study_sample)
sp <- as.data.frame(full_sp_mat[intersect(brca_samples, rownames(full_sp_mat)), "sig-21069"])
colnames(sp) <- "sig.21069"
genomic_ith_brca <- genomic_ith %>%
    filter(match_names %in% brca_samples) %>%
    select(match_names, genomic_ith) %>%
    column_to_rownames(var = "match_names")
cor_df <- sp %>%
    rownames_to_column(var = "match_names") %>%
    left_join(genomic_ith_brca %>% rownames_to_column(var = "match_names"), by = "match_names")

cor_plot <- ggscatter(
    cor_df,
    x = "genomic_ith",
    y = "sig.21069", # <- must match exactly the column name
    add = "reg.line",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "spearman",
    xlab = "WSE (Genomic Heterogeneity)",
    ylab = "VINCRISTINE"
)
ggsave(
    filename = "single_cell/heterogeneity/genomic_ith/figures/correlation_brca_vincristine_m.png",
    plot = cor_plot,
    width = 10, height = 6
)
