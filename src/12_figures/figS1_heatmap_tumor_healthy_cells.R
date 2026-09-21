library(BPCells)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)
set.seed(1234)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

# Load sketched matrix and metadata
mat <- open_matrix_dir(dir = "cnv_cells_genes_lvl2_fullbpcellsmatrix")

metadata <- read.table(
    "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_metadata_h5ad.tsv", 
    header = T, 
    sep = "\t"
    )

metadata <- metadata %>%
    mutate(
        cell = gsub("_|\\.|-", "", cell),
        barcode_study_sample = paste(cell, study, sample, sep = "__")
  ) %>%
  filter(!(is.na(malignancy)))

# Keep only cells with CNV data
metadata <- metadata[metadata$barcode_study_sample %in% colnames(mat), ]

######## SELECT SAMPLES TO PLOT HEATMAP OF MALIGNANT AND HEALTHY CELLS ########
# Parameters
N_CELLS_PER_SUBCLONE <- 20
N_HEALTHY_PER_SAMPLE <- 30
TOTAL_SAMPLES <- 100
MIN_CELLS <- 100
MIN_SUBCLONES <- 2

# Select samples from each tumor type based on sample proportions
tumor_proportions <- metadata %>%
    filter(malignancy == "True") %>%
    distinct(sample, tumor_type) %>%
    group_by(tumor_type) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(
        proportion = n_samples / sum(n_samples),
        n_to_select = pmax(1, round(proportion * TOTAL_SAMPLES))
    )

# Select samples with enough quality (number of cells and subclones)
valid_samples <- metadata %>%
    filter(malignancy == "True") %>%
    group_by(sample, tumor_type) %>%
    summarise(
        n_cells = n(),
        n_subclones = n_distinct(scevan_subclone),
        .groups = "drop"
    ) %>%
    filter(n_cells >= MIN_CELLS,
           n_subclones >= MIN_SUBCLONES)

# Select samples stratified by tumor type
selected_samples <- valid_samples %>%
    left_join(tumor_proportions, by = "tumor_type") %>%
    group_by(tumor_type) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), unique(.x$n_to_select)))) %>%
    ungroup()

# Subsample malignant cells by subclone
malignant_cells <- metadata %>%
    filter(sample %in% selected_samples$sample,
           malignancy == "True") %>%
    group_by(sample, scevan_subclone) %>%
    slice_sample(n = N_CELLS_PER_SUBCLONE, replace = FALSE) %>%
    ungroup()

# Subsample healthy cells by sample
healthy_cells <- metadata %>%
    filter(sample %in% selected_samples$sample,
           malignancy == "False") %>%
    group_by(sample) %>%
    slice_sample(n = N_HEALTHY_PER_SAMPLE, replace = FALSE) %>%
    ungroup()

# Combine and order
metadata_filtered <- bind_rows(malignant_cells, healthy_cells) %>%
    arrange(tumor_type, sample, malignancy, scevan_subclone)

# Subset CNV matrix and metadata to selected cells
mat <- mat[, metadata_filtered$barcode_study_sample]

# Convert BPCell matrix to dense matrix
mat <- as.matrix(mat)


################# REFORMAT GENE AND CHROMOSOME NAMES ################
process_matrix <- function(mat, bin_size = 3) {
  # Parse rownames
  x <- strsplit(rownames(mat), "_")
  x <- as.data.frame(do.call(rbind, x), stringsAsFactors = FALSE)

  lrg <- filter(x, V5 == "LRG") %>%
    unite(newV4, V4, V5, sep = "_")
  colnames(lrg) <- c("chromosome", "start", "end", "name", "symbol")

  not_lrg <- filter(x, V4 != "LRG")
  colnames(not_lrg) <- c("chromosome", "start", "end", "name", "symbol")

  x <- bind_rows(lrg, not_lrg[, 1:5]) %>%
    mutate(
      chromosome_numeric = as.integer(sub("X", "", chromosome)),
      start = as.numeric(start),
      original_name = paste(chromosome, start, end, name, symbol, sep = "_")
    ) %>%
    arrange(chromosome_numeric, start)

  # Reorder and transpose
  mat <- mat[match(x$original_name, rownames(mat)), ]
  rownames(mat) <- x$symbol
  mat <- t(mat)

  stopifnot(all(colnames(mat) == x$symbol))

  # Bin genes
  grouped_list <- list()
  grouped_names <- character()

  for (chr in unique(x$chromosome_numeric)) {
    idx <- which(x$chromosome_numeric == chr)
    chr_mat <- mat[, idx, drop = FALSE]
    chr_names <- colnames(chr_mat)

    for (i in seq(1, ncol(chr_mat), by = bin_size)) {
      end_idx <- min(i + bin_size - 1, ncol(chr_mat))
      cols <- chr_mat[, i:end_idx, drop = FALSE]

      grouped_list[[length(grouped_list) + 1]] <-
        if (ncol(cols) == 1) cols[, 1] else rowMeans(cols)

      grouped_names <- c(
        grouped_names,
        paste0(chr_names[i], "_to_", chr_names[end_idx])
      )
    }
  }

  grouped_mat <- do.call(cbind, grouped_list)
  colnames(grouped_mat) <- grouped_names
  grouped_mat
}

grouped_mat <- process_matrix(mat)

# Remove some malignant cells with few CNV for plotting
row_sd <- apply(grouped_mat, 1, sd)

malignant <- metadata_filtered$malignancy == "True"
rows_to_remove <- which(row_sd < 0.075 & malignant)

grouped_mat <- grouped_mat[-rows_to_remove, ]
metadata_filtered <- metadata_filtered[-rows_to_remove, ]

## Set heatmap annotation
## Clones annotation
clones_annot <- metadata_filtered %>% dplyr::select(barcode_study_sample, malignancy, scevan_subclone, tumor_type) %>%
    mutate(malignancy = ifelse(malignancy == "True", "Malignant", "Healthy"),
           scevan_subclone = ifelse(malignancy == "Healthy", "", scevan_subclone)) %>%
    arrange(malignancy, tumor_type, scevan_subclone) %>%
    column_to_rownames("barcode_study_sample")
grouped_mat <- grouped_mat[rownames(clones_annot), ]


pals <- list(
  malignancy = c(
    Malignant = "#db4646",
    Healthy = "#0090ab"
  ),
  tumor_type = tumor_type_colors
)

subclone_levels <- unique(clones_annot$scevan_subclone[clones_annot$malignancy == "Malignant"])
subclone_idx <- as.integer(factor(clones_annot$scevan_subclone, levels = subclone_levels))

cell_colors <- case_when(
  clones_annot$malignancy == "Healthy" ~ "#64bcda",
  subclone_idx %% 2 == 0 ~ "#bdbdbd",
  TRUE ~ "#969696"
)

left_annotation <- ComplexHeatmap::HeatmapAnnotation(
  which = "row",
  malignancy = clones_annot$malignancy,
  tumor_type = clones_annot$tumor_type,
  subclone = anno_simple(
    x = cell_colors,
    col = c(
      "#64bcda" = "#64bcda",
      "#bdbdbd" = "#bdbdbd",
      "#969696" = "#969696"
    ),
    gp = gpar(col = NA)
  ),
  col = pals,
  show_legend = TRUE,
  annotation_width = unit(c(3, 3, 3), "mm"),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10),
  annotation_label = c(
    malignancy = "Malignancy",
    tumor_type = "Tumor type",
    subclone = "Subclone"
  )
)

# Gene groups annotations
grouped_symbols <- sapply(strsplit(colnames(grouped_mat), "_to_"), `[`, 1)

genes_annot <- split_names_df %>%
  filter(symbol %in% grouped_symbols) %>%
  distinct(symbol, chromosome_numeric) %>%
  column_to_rownames("symbol") %>%
  dplyr::select(chromosome_numeric)

genes_annot <- genes_annot[grouped_symbols, , drop = FALSE]
colnames(genes_annot) <- "Chromosome"
rownames(genes_annot) <- colnames(grouped_mat)

## Assign colors: dark gray for odd, light gray for even
block_labels <- 1:22
block_colors <- ifelse(block_labels %% 2 == 0, "#bdbdbd", "#e2e2e2")

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = block_colors), labels = block_labels)
)

# Heatmap color
cnv_colors <- circlize::colorRamp2(
  breaks = seq(-0.3, 0.3, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))
)

# Order malignant cells by subclone
hm <- ComplexHeatmap::Heatmap(
  mat = grouped_mat,
  col = cnv_colors,
  name = "Segm.mean",
  left_annotation = left_annotation,
  top_annotation = top_annotation,
  row_split = factor(clones_annot$malignancy, levels = c("Malignant", "Healthy")),
  
  # Rows
  cluster_rows = FALSE,
  row_title = "Cells",
  show_row_names = FALSE,
  row_title_gp = gpar(fontface = "bold", fontsize = 12),
  
  # Cols
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_title = "Chromosomes",
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  column_split = genes_annot$Chromosome,
  
  heatmap_legend_param = list(
    title = "Segm.mean\nLog₂(ratio)",
    title_gp = gpar(fontface = "bold")
  ),
  heatmap_width = unit(6, "in"),
  heatmap_height = unit(10, "in"),
)

png(
  file = "../seurat/tcca/heatmap_tumor_healthy_cells.png",
  res = 300,
  width = 10,
  height = 20,
  units = "in"
)

draw(hm)
dev.off()


##################### Zoom of the heatmap for one sample #######################
clones_annot_sample <- metadata %>%
    filter(malignancy == "True" & sample == "CID44971") %>%
    select(barcode_study_sample, scevan_subclone) %>%
    mutate(scevan_subclone = as.character(str_extract(scevan_subclone, "\\d+$"))) %>%
    column_to_rownames("barcode_study_sample")

mat <- open_matrix_dir(dir = "cnv_cells_genes_lvl2_fullbpcellsmatrix")
mat <- mat[, rownames(clones_annot_sample)]
mat <- as.matrix(mat)

grouped_mat_sample <- process_matrix(mat)

# Remove some malignant cells with few CNV for plotting
row_sd <- apply(grouped_mat_sample, 1, sd)

rows_to_remove <- which(row_sd < 0.075)

# grouped_mat_sample <- grouped_mat_sample[-rows_to_remove, ]
# clones_annot_sample <- clones_annot_sample[-rows_to_remove, ]


# Plot heatmap for the sample
left_annotation <- ComplexHeatmap::HeatmapAnnotation(
  which = "row",
  subclone = clones_annot_sample$scevan_subclone,
  col = list(
    subclone = c(
    "1" = "#a8d8ea", 
    "2" = "#ffd3b6", 
    "3" = "#a8e6cf", 
    "4" = "#ffaaa5", 
    "5" = "#d4b8e0", 
    "6" = "#fff3b0"
    )),
  show_legend = FALSE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10),
  annotation_label = c(
    scevan_subclone = "Subclone"
  )
)

# Gene groups annotations
grouped_symbols <- sapply(strsplit(colnames(grouped_mat), "_to_"), `[`, 1)

genes_annot <- split_names_df %>%
  filter(symbol %in% grouped_symbols) %>%
  distinct(symbol, chromosome_numeric) %>%
  column_to_rownames("symbol") %>%
  dplyr::select(chromosome_numeric)

genes_annot <- genes_annot[grouped_symbols, , drop = FALSE]
colnames(genes_annot) <- "Chromosome"
rownames(genes_annot) <- colnames(grouped_mat)

## Assign colors: dark gray for odd, light gray for even
block_labels <- 1:22
top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = NA, col = NA), labels = block_labels)
)

# Heatmap color
cnv_colors <- circlize::colorRamp2(
  breaks = seq(-0.3, 0.3, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))
)

# Order malignant cells by subclone
hm <- ComplexHeatmap::Heatmap(
  mat = grouped_mat_sample,
  col = cnv_colors,
  name = "Segm.mean",
  left_annotation = left_annotation,
  top_annotation = top_annotation,
  row_split = factor(clones_annot_sample$scevan_subclone, levels = as.character(seq(1:6))),
  
  # Rows
  cluster_rows = FALSE,
  row_title = "Cells",
  show_row_names = FALSE,
  row_title_gp = gpar(fontface = "bold", fontsize = 12),
  
  # Cols
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_title = "Chromosomes",
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  column_split = genes_annot$Chromosome,
  
  heatmap_legend_param = list(
    title = "Segm.mean\nLog₂(ratio)",
    title_gp = gpar(fontface = "bold", fontsize = 11)
  ),
  heatmap_width = unit(7, "in"),
  heatmap_height = unit(5, "in"),
)

png(
  file = "../seurat/tcca/subclones_sample_brca.png",
  res = 300,
  width = 9,
  height = 7,
  units = "in"
)

draw(hm)
dev.off()