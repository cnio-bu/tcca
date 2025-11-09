library(Seurat)
library(BPCells)
library(qlcMatrix)
library(dplyr)
library(ComplexHeatmap)
## Source TCCA palette
source(file = "~/bc-meta/figures/TCCA_palette.R")

setwd("/storage/scratch01/shared/projects/bc-meta/")
mat_bc <- open_matrix_dir(dir = "./beyondcell_immuno/full_mat_beyondcell")

mat_mps <- read.table(
    "functional_nmf/sample_wise/mps_ucell_scores.tsv",
    header = TRUE
)
# Transform to sparse matrix
mat_bc <- as(mat_bc, "sparseMatrix")
mat_mps <- as(mat_mps, "sparseMatrix")

# Rename cells
colnames(mat_bc) <- paste0("c", c(1:ncol(mat_bc)))
rownames(mat_mps) <- paste0("c", c(1:nrow(mat_mps)))

# Transpose mat_bc (to cells x drugs)
mat_bc <- t(mat_bc)

# Compute correlation between drug sensitivity and mp enrichment
cor_mat <- corSparse(mat_bc, mat_mps)
rownames(cor_mat) <- colnames(mat_bc)
colnames(cor_mat) <- gsub("_UCell", "", colnames(mat_mps))

# Heatmap
## get drug names
drugs <- data.table::fread("~/bc-meta/reference/final_moas - Collapsed.tsv") %>%
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
    distinct() %>%
    as.data.frame()

drugs <- drugs[drugs$IDs %in% rownames(cor_mat), ]
rownames(drugs) <- drugs$IDs

drugs$collapsed.MoAs <- ifelse(
    drugs$collapsed.MoAs %in% names(MoAs_colors),
    drugs$collapsed.MoAs,
    "Other"
)
MoAs <- drugs[rownames(cor_mat), c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"

# Add mechanism of action to immunotherapies
MoAs$`Mechanism of action`[is.na(MoAs$`Mechanism of action`)] <- "Immunotherapy"

moa_pals <- list(
    "Mechanism of action" = MoAs_colors
)


right_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = MoAs,
    which = "row",
    col = moa_pals,
    show_annotation_name = FALSE,
    show_legend = FALSE
)


pdf(
    file = "functional_nmf/correlation.pdf",
    width = 17,
    height = 19
)


heat <- ComplexHeatmap::Heatmap(
    mat = cor_mat,
    # mat = t(sketched_mat),
    right_annotation = right_annotation,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = MoAs$`Mechanism of action`,
    #row_title = NULL,
    row_title_side = "right",
    row_title_rot = 0,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = TRUE,
    row_labels = drugs[rownames(cor_mat), "preferred.drug.names"],
    show_row_names = FALSE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_legend_param = list(title = "Pearson coefficient", direction = "horizontal"),
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(16, "in")
)


draw(heat, annotation_legend_side = "top", heatmap_legend_side = "bottom")
dev.off()

