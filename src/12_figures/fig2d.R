library(GenVisR)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(AnnotationHub)
library(dplyr)
library(ComplexHeatmap)

# Load CNV data
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata")
cs <- read.table("cnv_segments_clones_lvl2_cytobands.tsv", header = T, sep = "\t")

# This is a table of averages of amplification and deletions of each region
amp <- apply(cs, 1, function(row) mean(row[row > 2], na.rm = TRUE))
del <- apply(cs, 1, function(row) mean(row[row < 2], na.rm = TRUE))
amp.del <- cbind(amp, del)
amp.del <- replace(amp.del, is.nan(amp.del), 0)

# Make a granges of the cytobands of the human genome
proxy <- "mgonzalezb@cnio.es"
proxy <- httr::use_proxy(Sys.getenv("http_proxy"))
httr::set_config(proxy)
AnnotationHub::setAnnotationHubOption("PROXY", proxy)
AnnotationHub::getAnnotationHubOption("LOCAL")

hub <- AnnotationHub()

hub_hg38 <- subset(
    hub,
    (hub$species == "Homo sapiens") & (hub$genome == "hg38")
)

cytobands <- hub_hg38[[797]]
cytobands$custom_name <- paste0(seqnames(cytobands), cytobands$name)

dfA <- data.frame(
    names = cytobands$custom_name,
    chrom = as.data.frame(seqnames(cytobands))
)
dfB <- as.data.frame(ranges(cytobands))

cn.data <- cbind(dfA, dfB)
cn.data <- cn.data[, -5]

# Join cytobands to CNV data and set long format
cs <- rownames_to_column(as.data.frame(cs), var = "names")
cs <- merge(cn.data, cs, by = "names")

cs.long <- cs %>%
    pivot_longer(-c(names, value, start, end), names_to = "sampleID", values_to = "mean")

cs.long <- data.frame(
    region_name = cs.long$names,
    sample = cs.long$sampleID,
    chromosome = cs.long$value,
    start = cs.long$start,
    end = cs.long$end,
    probes = 1,
    segmean = cs.long$mean
)

cs.long <- replace(cs.long, is.na(cs.long), 2) %>%
    mutate(sample = str_replace_all(sample, "\\.", "-")) # This reverts scevan name changing

cs.long$sample <- gsub("GSM5645908_Breast_1_biol-rep", "GSM5645908_Breast_1_biol.rep", cs.long$sample) # This three fix three exceptions
cs.long$sample <- gsub("SyS11-met", "SyS11.met", cs.long$sample)
cs.long$sample <- gsub("Travaglini_Krasnow_2020_distal-1b", "Travaglini_Krasnow_2020_distal 1b", cs.long$sample)


# Load clinical and TCs metadata
metadata <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_metadata_h5ad.tsv", header = T, sep = "\t")
metadata <- metadata %>%
    mutate(study_sample = paste(study, sample, sep = "__")) %>%
    select(scevan_subclone, study, study_sample, tumor_type, sample_type, treated, therapeutic_cluster) %>%
    filter(!is.na(therapeutic_cluster)) %>%
    dplyr::rename(sample = scevan_subclone) %>%
    distinct()


# Add samples IDs and other metadata to cs.long
cs.long <- cs.long %>%
    mutate(sample = gsub("subclone", "", str_replace_all(sample, "__", "\\.")))
cs.long.merged <- left_join(cs.long, metadata, by = "sample")


# Plot heatmap of CNV data per TC
ht_data <- cs.long.merged %>%
    filter(study != "cell_lines_gabriella_kinker")

cnv_freq_list <- list()

for (tc in 1:10) {
    tc_samples <- ht_data %>% filter(therapeutic_cluster == tc)
    tc_freq <- cnFreq(tc_samples,
        CN_low_cutoff = 1.5,
        CN_high_cutoff = 2.5,
        genome = "hg38",
        out = "data"
    )
    tc_freq$data$TC <- paste0("TC", tc)
    cnv_freq_list[[paste0("TC", tc)]] <- tc_freq$data
}

# Combine frequencies for all
cnv_freq_all <- bind_rows(cnv_freq_list)


# Add region name
region_name <- ht_data %>%
    select(region_name, chromosome, start, end) %>%
    distinct()

cnv_freq_all <- cnv_freq_all %>%
    left_join(region_name, by = c("chromosome", "start", "end")) %>%
    mutate(chr_numeric = as.integer(gsub("chr", "", chromosome)))

region_order <- cnv_freq_all %>%
    distinct(region_name, chromosome, start, end, chr_numeric) %>%
    arrange(chr_numeric, start) %>%
    pull(region_name)

# Gain matrix (TCs x regions)
gain_matrix <- cnv_freq_all %>%
    dplyr::select(TC, region_name, gainProportion) %>%
    pivot_wider(
        names_from = region_name,
        values_from = gainProportion,
        values_fill = 0
    ) %>%
    column_to_rownames("TC") %>%
    as.matrix()

gain_matrix <- gain_matrix[paste0("TC", 1:10), region_order]

# Loss matrix
loss_matrix <- cnv_freq_all %>%
    dplyr::select(TC, region_name, lossProportion) %>%
    pivot_wider(
        names_from = region_name,
        values_from = lossProportion,
        values_fill = 0
    ) %>%
    column_to_rownames("TC") %>%
    as.matrix()

loss_matrix <- loss_matrix[paste0("TC", 1:10), region_order]

# Combined matrix: gains are positive, losses are negative
combined_matrix <- gain_matrix - loss_matrix

# Chromosome top annotation
region_chr <- cnv_freq_all %>%
    distinct(region_name, chr_numeric) %>%
    arrange(match(region_name, region_order)) %>%
    pull(chr_numeric)

block_colors <- ifelse(unique(region_chr) %% 2 == 0, "#bdbdbd", "#e2e2e2")
block_labels <- 1:22

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    foo = anno_block(
        gp = gpar(fill = block_colors),
        labels = block_labels,
        labels_gp = gpar(fontsize = 8)
    )
)

# Cytobands of interest to highlight
regions_interest <- c(
    "1q21.3", "12p13.2", "12q13.13", "17q21.2", "10p12.31",
    "17q21.1", "13q13.3", "13q14.11", "7p15.1", "7p21.1", "3q25.1", "8q24.3",
    "11q13.2", "19q13.2", "2q31.3", "3p11.1", "3q11.1", "5q14.1",
    "10p12.1", "13q12.11", "13q32.1", "21q22.3", "7p22.1", "7q31.1",
    "1q23.3", "6q23.2", "10p11.1", "10q11.1", "11q12.1", "17q11.1",
    "22q13.31"
)
regions_interest <- paste0("chr", regions_interest)
mark_idx <- which(region_order %in% regions_interest)
mark_labels <- region_order[mark_idx]

# Bottom annotation: mark regions of interest
bottom_annotation <- ComplexHeatmap::HeatmapAnnotation(
    mark = anno_mark(
        at = mark_idx,
        labels = mark_labels,
        labels_gp = gpar(fontsize = 7),
        side = "bottom",
        link_height = unit(3, "mm")
    ),
    which = "column"
)


# Heatmap color
cnv_colors <- circlize::colorRamp2(
    breaks = seq(-0.5, 0.5, length.out = 11),
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))
)


# Plot heatmap
hm_cnv <- Heatmap(
    matrix = combined_matrix,
    name = "CNV\nproportion",
    col = cnv_colors,
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,

    # Rows (TCs)
    cluster_rows = FALSE,
    row_order = paste0("TC", 1:10),
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),

    # Columns (regions)
    cluster_columns = FALSE,
    column_order = region_order,
    show_column_names = FALSE,
    column_split = region_chr,
    column_gap = unit(0.5, "mm"),
    column_title = NULL,

    # Legend
    heatmap_legend_param = list(
        title = "Net CNA\nproportion\n(gain - loss)",
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        direction = "vertical"
    ),
    width = unit(20, "cm"),
    height = unit(6, "cm")
)

pdf(
    file = "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/heatmap_cnv_TC.pdf",
    # res = 300,
    width = 14,
    height = 5,
    # units = "in"
)

draw(hm_cnv)

dev.off()



## COMPARE TCCA CNV profiles with TCGA CNV profiles #
# CNV freq per cancer type
cnv_freq_list <- list()

for (cancer_type in unique(ht_data$tumor_type)) {
    tc_samples <- ht_data %>% filter(tumor_type == cancer_type)
    tc_freq <- cnFreq(tc_samples,
        CN_low_cutoff = 1.5,
        CN_high_cutoff = 2.5,
        genome = "hg38",
        out = "data"
    )
    tc_freq$data$tumor_type <- cancer_type
    cnv_freq_list[[cancer_type]] <- tc_freq$data
}

# Combine frequencies for all
cnv_freq_all <- bind_rows(cnv_freq_list)
write.table(cnv_freq_all, "../seurat/compare_cnv_tcga/cnv_freq_tcca.tsv", sep = "\t")

# Download TCGA data
library(UCSCXenaTools)
library(data.table)
library(GenomicRanges)

# Download GDC PANCAN CNV dataset
tcga_cnv <- read.table("../seurat/compare_cnv_tcga/GDC-PANCAN.masked_cnv.tsv", sep = "\t", header = TRUE)
colnames(tcga_cnv) <- c("sample", "chromosome", "start", "end", "segmean")

# Extract bins from SCEVAN inferred CNVs
scevan_bins <- cnv_freq_all %>%
    distinct(chromosome, start, end)

# Proyectar TCGA a los bins de SCEVAN
project_to_scevan_bins <- function(sample_segs, bins_df) {
    seg_gr <- makeGRangesFromDataFrame(sample_segs,
        keep.extra.columns = TRUE,
        seqnames.field = "chromosome"
    )
    bins_gr <- makeGRangesFromDataFrame(bins_df,
        seqnames.field = "chromosome"
    )

    hits <- findOverlaps(bins_gr, seg_gr)

    result <- bins_df
    result$segmean <- 0 # neutral por defecto
    result$segmean[queryHits(hits)] <- sample_segs$segmean[subjectHits(hits)]
    result$sample <- unique(sample_segs$sample)
    result
}

tcga_binned <- tcga_cnv %>%
    mutate(chromosome = paste0("chr", chromosome)) %>%
    group_by(sample) %>%
    group_modify(~ project_to_scevan_bins(.x, scevan_bins)) %>%
    ungroup()

# Add cancer type
metadata <- read.table("../seurat/compare_cnv_tcga/tcga_metadata.tsv", sep = "\t", header = TRUE)
metadata <- metadata %>% select(barcode, cancer.type.abbreviation)
tcga_binned <- tcga_binned %>%
    left_join(metadata, by = c("sample" = "barcode")) %>%
    dplyr::rename(tumor_type = cancer.type.abbreviation) %>%
    filter(!is.na(tumor_type))

# Now cnFreq witht he same bins used in TCCA
cnv_freq_tcga_list <- list()

for (cancer_type in unique(tcga_binned$tumor_type)) {
    ct_samples <- tcga_binned %>% filter(tumor_type == cancer_type)

    ct_freq <- cnFreq(ct_samples,
        CN_low_cutoff = -0.2,
        CN_high_cutoff = 0.2,
        genome = "hg38",
        out = "data"
    )
    ct_freq$data$tumor_type <- cancer_type
    cnv_freq_tcga_list[[cancer_type]] <- ct_freq$data
}

cnv_freq_tcga <- bind_rows(cnv_freq_tcga_list)
write.table(cnv_freq_tcga, "../seurat/compare_cnv_tcga/cnv_freq_tcga.tsv", sep = "\t")

# Verify that bins match
identical(
    scevan_bins %>% arrange(chromosome, start),
    cnv_freq_tcga %>% distinct(chromosome, start, end) %>% arrange(chromosome, start)
)


# Create matrices with rows = bins, columns = cancer_type
cnv_freq_all <- read.table("../seurat/compare_cnv_tcga/cnv_freq_tcca.tsv", sep = "\t", header = TRUE)
cnv_freq_tcga <- read.table("../seurat/compare_cnv_tcga/cnv_freq_tcga.tsv", sep = "\t", header = TRUE)
make_cnv_matrix <- function(cnv_freq, group_col = "tumor_type") {
    cnv_freq %>%
        mutate(
            net = gainProportion - lossProportion,
            bin = paste(chromosome, start, end, sep = "_")
        ) %>%
        dplyr::select(bin, all_of(group_col), net) %>%
        pivot_wider(names_from = all_of(group_col), values_from = net) %>%
        column_to_rownames("bin") %>%
        as.matrix()
}

scevan_mat <- make_cnv_matrix(cnv_freq_all, group_col = "tumor_type")
tcga_mat <- make_cnv_matrix(cnv_freq_tcga, group_col = "tumor_type")

# Align bins (just in case)
common_bins <- intersect(rownames(scevan_mat), rownames(tcga_mat))
scevan_mat <- scevan_mat[common_bins, ]
tcga_mat <- tcga_mat[common_bins, ]

# Compute cosine similarity: CNV profiles per cancer type in TCCA vs CNV profiles
# per cancer type in TCGA
library(lsa)
common_ct <- intersect(colnames(scevan_mat), colnames(tcga_mat))
sim_matrix <- matrix(
    NA,
    nrow = length(common_ct),
    ncol = length(common_ct),
    dimnames = list(common_ct, common_ct)
)

for (ct_tcca in common_ct) {
    for (ct_tcga in common_ct) {
        sim_matrix[ct_tcca, ct_tcga] <- cosine(scevan_mat[, ct_tcca], tcga_mat[, ct_tcga])
    }
}

# Plot heatmap of similarities
library(ComplexHeatmap)
library(circlize)

sim_colors <- colorRamp2(
    seq(0, 1, length.out = 9),
    RColorBrewer::brewer.pal(9, "YlGnBu")
)

png(
    file = "../seurat/compare_cnv_tcga/heatmap_similarity_cnv_tcca_vs_tcga.png",
    res = 300,
    width = 8,
    height = 8,
    units = "in"
)
Heatmap(
    sim_matrix,
    name = "Cosine\nsimilarity",
    col = sim_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_rot = 45,
    row_title = "Your cohort (TC)",
    column_title = "TCGA cancer types",
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(round(sim_matrix[i, j], 2),
            x, y,
            gp = gpar(fontsize = 7)
        )
    }
)
dev.off()


# Function for a mirrored panel
# Table with chromosome positions
chr_order <- c(as.character(1:22))
chr_breaks <- cnv_freq_all %>%
    filter(tumor_type == unique(tumor_type)[1]) %>%
    mutate(chr_numeric = factor(gsub("chr", "", chromosome), levels = chr_order)) %>%
    arrange(chr_numeric, start) %>%
    mutate(bin_idx = row_number()) %>%
    group_by(chr_numeric) %>%
    summarise(bin_start = min(bin_idx)) %>%
    arrange(bin_start) %>%
    mutate(chromosome = as.character(chr_numeric))

chr_bands <- chr_breaks %>%
    arrange(bin_start) %>%
    mutate(
        bin_end = lead(bin_start, default = max(bin_start) + 50),
        band_fill = ifelse(row_number() %% 2 == 0, "gray95", "white")
    )

plot_cnv_mirror <- function(tcga_freq, tcca_freq, cancer_type, show_y = TRUE) {
    prep_data <- function(freq_df, source_label) {
        freq_df %>%
            filter(tumor_type == cancer_type) %>%
            mutate(chr_numeric = factor(gsub("chr", "", chromosome), levels = chr_order)) %>%
            arrange(chr_numeric, start) %>%
            mutate(bin_idx = row_number()) %>%
            dplyr::select(bin_idx, gainProportion, lossProportion) %>%
            pivot_longer(
                cols = c(gainProportion, lossProportion),
                names_to = "type",
                values_to = "freq"
            ) %>%
            mutate(
                freq   = ifelse(type == "lossProportion", -freq, freq),
                type   = ifelse(type == "gainProportion", "Gain", "Loss"),
                source = source_label
            )
    }

    tcga_data <- prep_data(tcga_freq, "TCGA")
    tcca_data <- prep_data(tcca_freq, "TCCA")

    cos_val <- round(sim_matrix[cancer_type, cancer_type], 2)

    make_panel <- function(data, title_label) {
        ggplot(data, aes(y = bin_idx)) +
            geom_rect(
                data = chr_bands,
                aes(xmin = -1, xmax = 1, ymin = bin_start, ymax = bin_end, fill = band_fill),
                inherit.aes = FALSE
            ) +
            scale_fill_identity() +
            geom_segment(
                data = filter(data, type == "Loss"),
                aes(x = 0, xend = freq, yend = bin_idx),
                color = "#4575b4", linewidth = 3
            ) +
            geom_segment(
                data = filter(data, type == "Gain"),
                aes(x = 0, xend = freq, yend = bin_idx),
                color = "#d73027", linewidth = 3
            ) +
            geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
            scale_x_continuous(
                limits = c(-1, 1),
                breaks = c(-1, 0, 1),
                labels = c("1", "0", "1")
            ) +
            scale_y_reverse(
                breaks = chr_breaks$bin_start,
                labels = chr_breaks$chromosome
            ) +
            labs(title = title_label, x = NULL, y = NULL) +
            theme_minimal(base_size = 7) +
            theme(
                panel.grid   = element_blank(),
                axis.text.y  = if (show_y) element_text(size = 5) else element_blank(),
                axis.ticks.y = element_blank(),
                plot.title   = element_text(hjust = 0.5, size = 7, color = "gray40")
            )
    }

    p_tcca <- make_panel(tcca_data, paste0(cancer_type, "\ncos=", cos_val))
    p_tcga <- make_panel(tcga_data, "")

    tcca_data / tcga_data
}

# Generate for all shared cancer types
shared_types <- intersect(
    unique(cnv_freq_tcga$tumor_type),
    unique(cnv_freq_all$tumor_type)
)

top_cancer_types <- c("LUAD", "LUSC", "BRCA", "SKCM", "ESCA", "GBM", "PAAD", "COAD", "KIRC", "OV")
plots <- lapply(seq_along(top_cancer_types), function(i) {
    plot_cnv_mirror(
        tcca_freq = cnv_freq_all,
        tcga_freq = cnv_freq_tcga,
        cancer_type = top_cancer_types[i],
        show_y = (i == 1)
    )
})

# Combinar todos en un panel
plot <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
        title = "CNV frequency profiles: TCGA (top) vs TCCA (bottom)",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
ggsave("../seurat/compare_cnv_tcga/cnv_mirror_top_comparison.pdf", width = 12, height = 4)
