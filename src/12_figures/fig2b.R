library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)

setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/cohort_statistics/")
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")

## Plot top drugs per cluster ##
drug_pred <- read.table("drug_response_subclone_final.tsv", header = TRUE)

# Get number of subclones per cluster
cluster_sizes <- drug_pred %>%
  select(Subclone.Name, ScTherapy.Cluster) %>%
  distinct() %>%
  count(ScTherapy.Cluster, name = "n_total")

# Frequency of each drug per therapeutic cluster
drug_freq <- drug_pred %>%
    group_by(ScTherapy.Cluster, Drug.Name) %>%
    summarise(
        n_predicted = n_distinct(Subclone.Name),
        .groups = "drop"
    ) %>%
    left_join(cluster_sizes, by = "ScTherapy.Cluster") %>%
    mutate(frequency = n_predicted / n_total)

# Proportions test
drug_exclusivity_pval <- drug_freq %>%
    group_by(Drug.Name) %>%
    mutate(
        n_other = sum(n_predicted) - n_predicted,
        total_other = sum(n_total) - n_total,
        pval = pmap_dbl(
            list(n_predicted, n_total, n_other, total_other),
            function(x, n, x_other, n_other) {
                # Avoid errors when the drug is not present in other TCs
                if (n_other == 0 || n == 0) return(1)
                tryCatch(
                    prop.test(
                        x = c(x, x_other),
                        n = c(n, n_other),
                        alternative = "greater"
                    )$p.value,
                    error = function(e) 1
                )
            }
        ),
        padj = p.adjust(pval, method = "BH")
    ) %>%
    ungroup()

# Filter drugs with frequency >= 0.5 in at least one cluster
selected_drugs <- drug_exclusivity_pval %>%
    group_by(Drug.Name) %>%
    filter(max(frequency) >= 0.50) %>%
    ungroup()

# Add MoA information
drug_moas <- drug_pred %>%
  select(Drug.Name, Drug.Mechanism.Of.Action) %>%
  distinct()

selected_drugs <- selected_drugs %>%
    left_join(drug_moas, by = "Drug.Name")

tc_order <- as.character(sort(unique(as.numeric(
    selected_drugs$ScTherapy.Cluster))))

# Prepare matrices for dotplot
# Frequency matrix (color)
freq_matrix <- selected_drugs %>%
    dplyr::select(Drug.Name, ScTherapy.Cluster, frequency) %>%
    pivot_wider(
        names_from = ScTherapy.Cluster,
        values_from = frequency,
        values_fill = 0
    ) %>%
    column_to_rownames("Drug.Name") %>%
    as.matrix()

freq_matrix <- freq_matrix[, tc_order]

# Adjusted p-values matrix (asterisk)
pval_matrix <- selected_drugs %>%
    dplyr::select(Drug.Name, ScTherapy.Cluster, padj) %>%
    pivot_wider(
        names_from = ScTherapy.Cluster,
        values_from = padj,
        values_fill = 1
    ) %>%
    column_to_rownames("Drug.Name") %>%
    as.matrix()

pval_matrix <- pval_matrix[, tc_order]

# Order drugs by MoA and then by name
drug_order <- selected_drugs %>%
    distinct(Drug.Name, Drug.Mechanism.Of.Action) %>%
    mutate(Drug.Mechanism.Of.Action = factor(
        Drug.Mechanism.Of.Action,
        levels = c(
            sort(setdiff(unique(Drug.Mechanism.Of.Action), "Other")),
            "Other"  # Other siempre al final
        )
    )) %>%
    arrange(Drug.Mechanism.Of.Action, Drug.Name) %>%
    pull(Drug.Name)

drug_order <- drug_order[drug_order %in% rownames(freq_matrix)]
freq_matrix <- freq_matrix[drug_order, ]
pval_matrix <- pval_matrix[drug_order, ]

# Rename therapeutic clusters
colnames(freq_matrix) <- paste0("TC", colnames(freq_matrix))
colnames(pval_matrix) <- paste0("TC", colnames(pval_matrix))
tc_order <- paste0("TC", 1:10)

# Frequency values control bubble size
size_matrix <- freq_matrix * 2 + 0.5

# Color scale for frequency values
freq_col <- colorRamp2(
    seq(0, 1, length.out = 9),
    RColorBrewer::brewer.pal(9, "Blues")
)

# Row annotation for MoA
moa_annot <- selected_drugs %>%
    distinct(Drug.Name, Drug.Mechanism.Of.Action) %>%
    arrange(match(Drug.Name, drug_order)) %>%
    filter(Drug.Name %in% drug_order) %>%
    pull(Drug.Mechanism.Of.Action)

left_annotation <- rowAnnotation(
    MoA = moa_annot,
    col = list(MoA = MoAs_colors),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.3, "cm"),
    show_legend = FALSE
)

# Function to draw bubbles in the heatmap cells
cell_fun <- function(j, i, x, y, width, height, fill) {
    # Size and color based on frequency
    size <- size_matrix[i, j]
    col <- freq_col(freq_matrix[i, j])
    
    # Draw bubble
    grid.circle(
        x, y,
        r = unit(size, "mm"),
        gp = gpar(fill = col, col = "grey80", lwd = 0.5)
    )
    
    # Add asterisk if padj < 0.05
    if (pval_matrix[i, j] < 0.05) {
        grid.text(
            "*",
            x, y,
            vjust = 0.75,
            gp = gpar(fontsize = 8, col = "black", fontface = "bold")
        )
    }
}


# Create base heatmap
cell_size <- unit(0.6, "cm")
hm <- Heatmap(
    matrix = freq_matrix,
    name = "Frequency",
    cell_fun = cell_fun,
    rect_gp = gpar(col = "grey90", fill = NA, lwd = 0.5),
    left_annotation = left_annotation,
    
    # Filas
    cluster_rows = FALSE,
    row_order = drug_order,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    
    # Columnas
    cluster_columns = FALSE,
    column_order = tc_order,
    column_title = "Therapeutic cluster",
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    show_column_names = TRUE,
    
    show_heatmap_legend = FALSE,
    
    width = cell_size * 10,
    height = cell_size * nrow(freq_matrix)
)

# Legend for frequency (color)
freq_col_legend <- Legend(
    title = "Proportion of\n subclones",
    col_fun = freq_col,
    direction = "vertical",
    at = seq(0, 1, by = 0.25),
    labels = c("0", "0.25", "0.50", "0.75", "1")
)


# Legend for MoA
present_moas <- unique(moa_annot)
moa_levels <- c(sort(setdiff(present_moas, "Other")), "Other")
filtered_moa_colors <- MoAs_colors[moa_levels]
filtered_moa_colors <- filtered_moa_colors[!is.na(filtered_moa_colors)]

moa_legend <- Legend(
    title = "MoA",
    labels = names(filtered_moa_colors),
    legend_gp = gpar(fill = filtered_moa_colors),
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
)

# Save the plot as a png file
pdf(
    file = "/Users/mariagb/Documents/new_figures_tcca/dotplot_drugs_TC.pdf",
    # res = 300,
    width = 8,
    height = 18,
    # units = "in"
)

draw(
    hm,
    annotation_legend_list = list(freq_col_legend, moa_legend),
    annotation_legend_side = "right",
    heatmap_legend_side = "right"
)

dev.off()