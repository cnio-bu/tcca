library(dplyr)
library(readr)
library(pheatmap)

setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf/")
out_dir <- "therapeutic_from_mps/sctherapy_mps/"

## 1. MP drug predictions
# Aggregate tables
files <- list.files(out_dir, pattern = "_mps_drugs\\.tsv$", full.names = TRUE)
full_table <- bind_rows(lapply(files, read_tsv))
write.table(full_table, paste0(out_dir, "full_table_drugs_mps.tsv"), sep = "\t")

# Outside the HPC cluster
setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/therapeutic_analysis/therapeutic_from_mps")

# Match MP names
mps <- readRDS("../../functional/functional_mps/sample_wise/prefilter_cpm/metaprograms_cpm/mp_list_reordered.rds")
mps_names <- names(mps)
names(mps_names) <- paste0("MP", 1:43)

full_table <- read.table("full_table_drugs_mps.tsv")
full_table <- full_table %>%
    mutate(
        mp = names(mps_names)[match(mp, mps_names)],
        Drug_Name = toupper(Drug_Name)
    )

full_table %>%
    distinct(sample, mp) %>%
    nrow()

# Drug frequency per MP
n_samples_per_mp <- full_table %>%
    group_by(mp) %>%
    summarise(n_samples_total = n_distinct(sample))

print(n_samples_per_mp, n = 43)

drug_freq_mp <- full_table %>%
    group_by(mp, Drug_Name) %>%
    summarise(n_samples_predicted = n_distinct(sample), .groups = "drop") %>%
    left_join(n_samples_per_mp, by = "mp") %>%
    mutate(pct_samples = n_samples_predicted / n_samples_total)

# Drug list per MP at 50% threshold
drug_list_mp <- drug_freq_mp %>%
    filter(pct_samples >= 0.5) %>%
    arrange(mp, desc(pct_samples))

write.table(drug_list_mp, "drug_list_by_mp.tsv", sep = "\t")

# Summary
mp_counts <- drug_list_mp %>%
    count(mp, name = "n_drugs") %>%
    arrange(desc(n_drugs))
write.table(mp_counts, "drug_counts_mps.tsv", sep = "\t")

# 2. TC drug predictions
# Cross-validation with therapeutic cluster consensus drugs
# Load drugs per TC
dir <- "/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/cohort_statistics/"
drug_pred <- read.table(paste0(dir, "drug_response_subclone_final.tsv"), header = TRUE) %>%
    mutate(ScTherapy.Cluster = paste0("TC", ScTherapy.Cluster))
cluster_sizes <- drug_pred %>%
    select(Subclone.Name, ScTherapy.Cluster) %>%
    distinct() %>%
    count(ScTherapy.Cluster, name = "n_total")

drug_freq_tc <- drug_pred %>%
    group_by(ScTherapy.Cluster, Drug.Name) %>%
    summarise(n_predicted = n_distinct(Subclone.Name), .groups = "drop") %>%
    left_join(cluster_sizes, by = "ScTherapy.Cluster") %>%
    mutate(pct_subclones = n_predicted / n_total)

# Drug list per TC at 50% threshold
drug_list_tc <- drug_freq_tc %>%
    filter(pct_subclones >= 0.5) %>%
    arrange(ScTherapy.Cluster, desc(pct_subclones))

drug_counts <- drug_list_tc %>%
    count(ScTherapy.Cluster, name = "n_drugs") %>%
    arrange(desc(n_drugs))
write.table(drug_counts, "drug_counts_tcs.tsv", sep = "\t")

## 3. Jaccard index MPs x TCs
jaccard <- function(set1, set2) {
    u <- length(union(set1, set2))
    if (u == 0) {
        return(NA_real_)
    } # NA instead of 0 when both sets are empty
    length(intersect(set1, set2)) / u
}
# Number of TCs in which each drug is part of the consensus set (>=50% of subclones)
broad_df <- drug_list_tc %>%
    mutate(Drug.Name = toupper(Drug.Name)) %>%
    group_by(Drug.Name) %>%
    summarise(n_tcs = n_distinct(ScTherapy.Cluster), .groups = "drop") %>%
    arrange(desc(n_tcs))

broad <- broad_df %>%
    filter(n_tcs >= 6) %>%
    pull(Drug.Name)

# Drug sets per MP and TC
mp_sets <- split(toupper(drug_list_mp$Drug_Name), drug_list_mp$mp)
tc_sets <- split(toupper(drug_list_tc$Drug.Name), drug_list_tc$ScTherapy.Cluster)

mp_names <- paste0("MP", 1:43)
tc_names <- paste0("TC", 1:10)

compute_jaccard <- function(mp_sets, tc_sets) {
    jaccard_matrix <- sapply(tc_names, function(tc) {
        sapply(mp_names, function(mp) jaccard(mp_sets[[mp]], tc_sets[[tc]]))
    })
    dimnames(jaccard_matrix) <- list(mp_names, tc_names)
    jaccard_matrix
}

# Full drug sets
jac_spec <- compute_jaccard(
    lapply(mp_sets, setdiff, broad),
    lapply(tc_sets, setdiff, broad)
)
jaccard_matrix <- compute_jaccard(mp_sets, tc_sets)
write.csv(jaccard_matrix, "jaccard_mp_tc.csv")

# Verify global overlap between MP and TC collection
drugs_mp <- unique(unlist(mp_sets))
drugs_tc <- unique(unlist(tc_sets))
cat("Unique drugs in MPs:", length(drugs_mp), "\n")
cat("Unique drugs in TCs:", length(drugs_tc), "\n")
cat("Overlap:", length(intersect(drugs_mp, drugs_tc)), "\n")

# Hypergeometric test for significance
universe_size <- length(union(
    unique(toupper(full_table$Drug_Name)),
    unique(toupper(drug_pred$Drug.Name))
)) # total compounds predicted in scTherapy

phyper_overlap <- function(set1, set2, universe = universe_size) {
    k <- length(intersect(set1, set2))
    m <- length(set1)
    n <- universe - m
    q <- length(set2)
    if (k == 0) {
        return(1)
    }
    phyper(k - 1, m, n, q, lower.tail = FALSE)
}

compute_pvalue <- function(mp_sets, tc_sets) {
    pval_matrix <- sapply(tc_names, function(tc) {
        sapply(mp_names, function(mp) phyper_overlap(mp_sets[[mp]], tc_sets[[tc]]))
    })
    dimnames(pval_matrix) <- list(mp_names, tc_names)
    pval_matrix
}

pval_matrix <- compute_pvalue(mp_sets, tc_sets)

# BH correction
pval_adj_matrix <- matrix(
    p.adjust(pval_matrix, method = "BH"),
    nrow = length(mp_names),
    ncol = length(tc_names),
    dimnames = list(mp_names, tc_names)
)
sig_matrix <- ifelse(pval_adj_matrix < 0.05, "*", "")

# Save Jaccard heatmap
source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")
# Row annotation
mp_families <- data.frame(
    MP = paste0("MP", 1:43),
    family = c(
        rep("CellCycle", 6),
        "Oncogenic",
        rep("Stress", 6),
        rep("Inflammation", 3),
        rep("EMT", 5),
        rep("CellularPlasticity", 3),
        rep("ProteinRegulation", 4),
        "EpithelialSenescence",
        "MitochondrialRespiration",
        "Cilia",
        rep("LineageSpecific.Hemato", 6),
        rep("LineageSpecific.Neural", 2),
        rep("LineageSpecific.Other", 4)
    )
)

mp_families <- mp_families %>%
    column_to_rownames("MP")

top_annotation <- HeatmapAnnotation(
    Family = mp_families$family,
    col = list(Family = mp_family_colors),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.3, "cm"),
    annotation_legend_param = list(
        Family = list(
            ncol = 3,
            direction = "vertical"
        )
    )
)

# Heatmap
col_fun <- colorRamp2(
    seq(0, 0.8, length.out = 100),
    colorRampPalette(c(
        "#ffffff", "#f2f0f7", "#cbc9e2",
        "#9e9ac8", "#756bb1", "#54278f"
    ))(100)
)
hm <- Heatmap(
    matrix = t(jaccard_matrix),
    name = "Jaccard index",
    col = col_fun,
    rect_gp = gpar(col = "grey90", lwd = 0.3),
    top_annotation = top_annotation,

    # Rows
    cluster_columns = FALSE,
    column_title = "Metaprogram",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title_side = "bottom",
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_names_centered = TRUE,

    # Columns
    cluster_rows = FALSE,
    row_title = "Therapeutic cluster",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_names_rot = 0,
    row_names_centered = TRUE,

    # Legend
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
        title = "Jaccard index",
        direction = "horizontal"
    )
)

pdf(
    file = "jaccard_mps_tcs.pdf",
    # res = 300,
    width = 12,
    height = 5,
    # units = "in"
)

draw(
    hm,
    annotation_legend_side = "bottom",
    heatmap_legend_side = "bottom",
    merge_legend = TRUE
)
dev.off()






################### Bubble plot with drug frequency per MP #####################

# 1. Standarize MoAs
MoAs <- read.table(
    "/Users/mariagb/Documents/tcca/reference/final_moas - Collapsed.tsv",
    header = TRUE,
    sep = "\t"
) %>%
    select(preferred.drug.names, collapsed.MoAs) %>%
    distinct()

MoAs_sctherapy <- full_table %>%
    select(Drug_Name, MoA) %>%
    distinct() %>%
    mutate(Drug_Name = toupper(Drug_Name)) %>%
    as.data.frame()

rownames(MoAs_sctherapy) <- MoAs_sctherapy$Drug_Name

full_table <- full_table %>%
    mutate(Drug_Name = toupper(Drug_Name)) %>%
    left_join(
        select(MoAs, preferred.drug.names, collapsed.MoAs),
        by = c("Drug_Name" = "preferred.drug.names")
    )

full_table <- full_table %>%
    mutate(collapsed.MoAs = ifelse(
        is.na(collapsed.MoAs),
        MoAs_sctherapy[Drug_Name, "MoA"],
        collapsed.MoAs
    ))

full_table <- full_table %>%
    mutate(collapsed.MoAs = case_when(
        collapsed.MoAs == "thymidylate_synthase_inhibitor" ~ "DNA related agent",
        collapsed.MoAs == "tubulin_polymerization_inhibitor" ~ "Microtubule agent",
        collapsed.MoAs == "MDM_inhibitor" ~ "p53 activator/MDM2 inhibitor",
        collapsed.MoAs == "ALK_tyrosine_kinase_receptor_inhibitor" ~ "Kinase inhibitor",
        collapsed.MoAs == "Brutons_tyrosine_kinase_BTK_inhibitor" ~ "Kinase inhibitor",
        collapsed.MoAs == "inosine_monophosphate_dehydrogenase_inhibitor" ~ "DNA related agent",
        collapsed.MoAs == "RAF_inhibitor" ~ "BRAF inhibitor",
        collapsed.MoAs == "BRD-K24576554" ~ "Other",
        collapsed.MoAs == "BRD-K60237333" ~ "Other",
        collapsed.MoAs == "BRD-K95142244" ~ "Other",
        collapsed.MoAs == "MAP_kinase_inhibitor" ~ "MAPK inhibitor",
        collapsed.MoAs == "adenosine_receptor_antagonist" ~ "Other",
        collapsed.MoAs == "centromere_associated_protein_inhibitor" ~ "Microtubule agent",
        collapsed.MoAs == "PLK_inhibitor" ~ "Cell cycle arrest",
        collapsed.MoAs == "HDAC_inhibitor" ~ "Chromatin agent",
        collapsed.MoAs == "CDK_inhibitor" ~ "Cell cycle arrest",
        collapsed.MoAs == "topoisomerase_inhibitor" ~ "DNA related agent",
        collapsed.MoAs == "DNA_synthesis_inhibitor_microtubule_inhibitor" ~ "Microtubule agent",
        collapsed.MoAs == "kinesin_inhibitor_kinesinlike_spindle_protein_inhibitor" ~ "Microtubule agent",
        collapsed.MoAs == "VEGFR_inhibitor" ~ "VEGFR inhibitor",
        collapsed.MoAs == "potassium_channel_activator" ~ "Other",
        collapsed.MoAs == "KIT_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor_VEGFR_inhibitor" ~ "VEGFR inhibitor",
        collapsed.MoAs == "dehydrogenase_inhibitor_inositol_monophosphatase_inhibitor" ~ "Other",
        collapsed.MoAs == "FLT3_inhibitor_JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
        collapsed.MoAs == "Ttype_calcium_channel_blocker" ~ "Other",
        collapsed.MoAs == "BcrAbl_kinase_inhibitor_FLT3_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor" ~ "BCR-ABL inhibitor",
        collapsed.MoAs == "estrogen_receptor_antagonist_selective_estrogen_receptor_modulator_SERM" ~ "Other",
        collapsed.MoAs == "retinoid_receptor_agonist" ~ "Other",
        collapsed.MoAs == "HSP_inhibitor" ~ "HSP inhibitor",
        collapsed.MoAs == "RNA_polymerase_inhibitor" ~ "Transcription inhibitor",
        collapsed.MoAs == "PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
        collapsed.MoAs == "JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
        collapsed.MoAs == "mTOR_inhibitor_PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
        collapsed.MoAs == "Aurora_kinase_inhibitor_FLT3_inhibitor_VEGFR_inhibitor" ~ "Other",
        collapsed.MoAs == "BRAF inhibitor;VEGFR inhibitor" ~ "Kinase inhibitor",
        collapsed.MoAs == "BCR-ABL inhibitor;SRC inhibitor" ~ "BCR-ABL inhibitor",
        collapsed.MoAs == "VEGFR inhibitor;MET inhibitor" ~ "Kinase inhibitor",
        collapsed.MoAs == "CDK_inhibitor_cell_cycle_inhibitor_MCL1_inhibitor" ~ "Cell cycle arrest",
        collapsed.MoAs == "NFkB_pathway_inhibitor_proteasome_inhibitor" ~ "Ubiquitin-proteasome system inhibitor",
        collapsed.MoAs == "proteasome_inhibitor" ~ "Ubiquitin-proteasome system inhibitor",
        collapsed.MoAs == "adenylyl_cyclase_activator_phosphodiesterase_inhibitor" ~ "Other",
        collapsed.MoAs == "-" ~ "Other",
        Drug_Name == "CYT387" ~ "JAK-STAT signaling inhibitor",
        TRUE ~ collapsed.MoAs
    ))


unique(full_table$collapsed.MoAs) %in% names(MoAs_colors)


# 2. Prepare data for bubble plot
# Proportions test
drug_freq_pval <- drug_freq_mp %>%
    group_by(Drug_Name) %>%
    mutate(
        n_other = sum(n_samples_predicted) - n_samples_predicted,
        total_other = sum(n_samples_total) - n_samples_total,
        pval = pmap_dbl(
            list(n_samples_predicted, n_samples_total, n_other, total_other),
            function(x, n, x_other, n_other) {
                # Avoid errors when the drug is not present in other TCs
                if (n_other == 0 || n == 0) {
                    return(1)
                }
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
selected_drugs <- drug_freq_pval %>%
    group_by(Drug_Name) %>%
    filter(max(pct_samples) >= 0.50) %>%
    ungroup()

# Add MoA information
drug_moas <- full_table %>%
    select(Drug_Name, collapsed.MoAs) %>%
    distinct()

selected_drugs <- selected_drugs %>%
    left_join(drug_moas, by = "Drug_Name")

mp_order <- paste0("MP", 1:43)

# Prepare matrices for dotplot
# Frequency matrix (color)
freq_matrix <- selected_drugs %>%
    dplyr::select(Drug_Name, mp, pct_samples) %>%
    pivot_wider(
        names_from = mp,
        values_from = pct_samples,
        values_fill = 0
    ) %>%
    column_to_rownames("Drug_Name") %>%
    as.matrix()

freq_matrix <- freq_matrix[, mp_order]

# Adjusted p-values matrix (asterisk)
pval_matrix <- selected_drugs %>%
    dplyr::select(Drug_Name, mp, padj) %>%
    pivot_wider(
        names_from = mp,
        values_from = padj,
        values_fill = 1
    ) %>%
    column_to_rownames("Drug_Name") %>%
    as.matrix()

pval_matrix <- pval_matrix[, mp_order]

# Order drugs by MoA and then by name
drug_order <- selected_drugs %>%
    distinct(Drug_Name, collapsed.MoAs) %>%
    mutate(collapsed.MoAs = factor(
        collapsed.MoAs,
        levels = c(
            sort(setdiff(unique(collapsed.MoAs), "Other")),
            "Other" # Other at the end
        )
    )) %>%
    arrange(collapsed.MoAs, Drug_Name) %>%
    pull(Drug_Name)

drug_order <- drug_order[drug_order %in% rownames(freq_matrix)]
freq_matrix <- freq_matrix[drug_order, ]
pval_matrix <- pval_matrix[drug_order, ]

# Frequency values control bubble size
size_matrix <- freq_matrix * 2 + 0.5

# Color scale for frequency values
freq_col <- colorRamp2(
    seq(0, 1, length.out = 9),
    RColorBrewer::brewer.pal(9, "Reds")
)

# Row annotation for MoA
moa_annot <- selected_drugs %>%
    distinct(Drug_Name, collapsed.MoAs) %>%
    arrange(match(Drug_Name, drug_order)) %>%
    filter(Drug_Name %in% drug_order) %>%
    pull(collapsed.MoAs)

left_annotation <- rowAnnotation(
    MoA = moa_annot,
    col = list(MoA = MoAs_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size = unit(0.3, "cm"),
    show_legend = FALSE
)

# Top annotation
top_annotation <- HeatmapAnnotation(
    Family = mp_families$family,
    col = list(Family = mp_family_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
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
    top_annotation = top_annotation,

    # Rows
    cluster_rows = FALSE,
    row_order = drug_order,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),

    # Columnas
    cluster_columns = FALSE,
    column_order = mp_order,
    column_title = "Metaprogram",
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    column_names_rot = 45,
    show_column_names = TRUE,
    show_heatmap_legend = FALSE,
    width = cell_size * ncol(freq_matrix) / 1.2,
    height = cell_size * nrow(freq_matrix) / 1.2
)

# Legend for frequency (color)
freq_col_legend <- Legend(
    title = "Proportion of\n samples",
    col_fun = freq_col,
    direction = "vertical",
    at = seq(0, 1, by = 0.25),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
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
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
)

# Legend for MPs
mp_legend <- Legend(
    title = "Functional family",
    labels = names(mp_family_colors),
    legend_gp = gpar(fill = mp_family_colors),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
)

# Save the plot as a png file
pdf(
    file = "dotplot_drugs_MP.pdf",
    # res = 300,
    width = 12,
    height = 20,
    # units = "in"
)

draw(
    hm,
    annotation_legend_list = list(mp_legend, moa_legend, freq_col_legend),
    annotation_legend_side = "bottom",
    heatmap_legend_side = "bottom"
)

dev.off()
