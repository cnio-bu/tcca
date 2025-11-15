library(Seurat)
library(BPCells)
library(GSVA)
library(GSEABase)
library(edgeR)
library(dplyr)
library(tidyverse)
library(uwot)
setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/")

# Load Seurat object with
seu.lvl2 <- readRDS("../single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
metadata <- read.table(
  "../single_cell/seurat/v5/tcca_metadata.tsv",
  sep = "\t",
  header = TRUE,
  row.names = NULL
) %>%
  column_to_rownames(var = "cell")

seu.lvl2@meta.data <- metadata

# Select only TCCA cells assigned to therapeutic clusters (TCs)
malignant_subset <- subset(
  seu.lvl2,
  subset = malignancy == "True" &
    scevan_subclone != "non_tumor" &
    !is.na(seu.lvl2$scTherapy_cluster)
)

# Pseudobulk at subclone level
tcca_bulk <- AggregateExpression(
  object = malignant_subset,
  slot = "counts",
  return.seurat = T,
  group.by = c("scevan_subclone")
)

# Compute Beyondcell enrichment in drug sensitivity signatures
mat <- tcca_bulk[["RNA"]]$counts
mat <- as.matrix(mat)

write.table(x = mat, file = "subclone_level/subclones_expr_mat.tsv", sep = "\t")

dge <- DGEList(counts = mat)
dge <- calcNormFactors(object = dge, method = "TMM")

## get norm. mat
mat_cpm <- cpm(
  dge,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 1
)

dim(mat_cpm)

# Read the GMT file with SSc and immunotherapy signatures
gsets <- GSEABase::getGmt(con = "../reference/bc_immuno.gmt")

## gsva parameters
gsvapar <- gsvaParam(
  exprData = mat_cpm,
  geneSets = gsets,
  kcdf = "Gaussian",
  maxDiff = TRUE
)

gsva_bc <- gsva(gsvapar)

write.table(gsva_bc, "subclone_level/subclones_bc_gsva.tsv", sep ="\t")

# For signatures with UP and DOWN gene sets compute the substraction of GSVA scores
signatures_UP <- rownames(gsva_bc)[grepl("_UP$", rownames(gsva_bc))]
signatures_DN <- rownames(gsva_bc)[grepl("_DOWN$", rownames(gsva_bc))]
base_UP <- sub("_UP$", "", signatures_UP)
base_DN <- sub("_DOWN$", "", signatures_DN)
bidirectional_sigs <- intersect(base_UP, base_DN)

final_matrix <- gsva_bc

for (base_name in bidirectional_sigs) {
  up_name <- paste0(base_name, "_UP")
  dn_name <- paste0(base_name, "_DOWN")
  
  # Compute difference
  final_matrix <- rbind(final_matrix, gsva_bc[up_name, ] - gsva_bc[dn_name, ])
  rownames(final_matrix)[nrow(final_matrix)] <- base_name
  
  # Remove UP and DOWN drug signature enrichment
  final_matrix <- final_matrix[!rownames(final_matrix) %in% c(up_name, dn_name), ]
}

write.table(final_matrix, "subclone_level/subclones_bc_gsva_aggregated.tsv", sep ="\t")


# Compute UMAP of subclones based on drug enrichment
gsva_final <- read.table("subclone_level/subclones_bc_gsva_aggregated.tsv", sep ="\t", header = TRUE)

# Perform an UWOT on the drug signatures scores
umap_transform <- uwot::umap(X = t(gsva_final),
                             n_neighbors = 20,
                             n_components = 2,
                             metric = "correlation",
                             n_threads = 4,
                             pca = NULL)


# Draw the therapeutic clusters over the umap
umap_transform <- as.data.frame(umap_transform)
umap_transform <- umap_transform %>%
    rownames_to_column(var = "subclone") %>%
    mutate(subclone = gsub("[^a-zA-Z0-9]", "", subclone)) %>%
    left_join(
        malignant_subset@meta.data %>%
            mutate(subclone = gsub("[^a-zA-Z0-9]", "", scevan_subclone), 
                   therapeutic_cluster = factor(scTherapy_cluster, levels = as.character(1:10))) %>%
            select(subclone, study, sample_type, treated, therapeutic_cluster) %>%
            distinct(),
        by = "subclone"
    ) %>%
    column_to_rownames(var = "subclone")

source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
tcs_umap <- ggplot(
    data = umap_transform,
    aes(x = V1, y = V2, color = therapeutic_cluster)) +
    geom_point(size = 3) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x.bottom = element_line(),
          axis.line.y.left   = element_line(),
          axis.line.y.right  = element_line(),
          axis.text.y.right  = element_blank(),
          axis.ticks.y.right = element_blank(),
          panel.border       = element_blank(),
          text = element_text(family = "Arial")
    ) +
    scale_color_manual(values = sctherapy_colors) +
    scale_x_continuous(name = "UMAP1") +
    scale_y_continuous(name = "UMAP2") +
    scale_color_discrete(name = "Therapeutic cluster") 

ggsave(
    filename = "subclone_level/plots/tcs_uwot.png", 
    plot = tcs_umap,
    dpi = 300,
    height = 4,
    width = 5
)


## Plot heatmap with top variance drugs
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 1. Prepare drug information
drugs <- data.table::fread("../reference/final_moas - Collapsed.tsv")
drug_info <- drugs %>%
    filter(IDs %in% rownames(gsva_final)) %>%
    group_by(IDs) %>%
    summarise(
        drug_name = paste0(first(preferred.drug.names), "_", first(studies)),
        collapsed.MoAs = paste(unique(collapsed.MoAs), collapse = ";"),  # Sin espacios
        .groups = "drop"
    ) %>%
    mutate(
        collapsed.MoAs = case_when(
            collapsed.MoAs == "VEGFR inhibitor;MET inhibitor" ~ "Kinase inhibitor",
            collapsed.MoAs == "BCR-ABL inhibitor;SRC inhibitor" ~ "BCR-ABL inhibitor",
            collapsed.MoAs == "BRAF inhibitor;VEGFR inhibitor" ~ "Kinase inhibitor",
            TRUE ~ collapsed.MoAs
        )
    ) %>%
    column_to_rownames("IDs")

# 2. Prepare drug subset and scale gsva matrix
top_rv <- matrixStats::rowVars(as.matrix(gsva_final))
top50_drugs <- names(sort(top_rv, decreasing = TRUE)[1:50])

gsva_subset <- gsva_final[top50_drugs, ]
gsva_scaled <- t(scale(t(gsva_subset), center = TRUE, scale = TRUE))

new_names <- ifelse(
    rownames(gsva_scaled) %in% rownames(drug_info),
    drug_info[rownames(gsva_scaled), "drug_name"],
    rownames(gsva_scaled)
)
rownames(gsva_scaled) <- new_names
colnames(gsva_scaled) <- gsub("[^a-zA-Z0-9]", "", colnames(gsva_scaled))


# 3. Column annotations
annotation_col <- umap_transform %>%
    rownames_to_column("subclone") %>%
    select(subclone, study, sample_type, treated, therapeutic_cluster) %>%
    mutate(
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        treated = ifelse(treated == "t", "Treated", "Untreated")
    ) %>%
    column_to_rownames("subclone")

col_ha <- HeatmapAnnotation(
    Study = annotation_col$study,
    `Sample type` = annotation_col$sample_type,
    Treated = annotation_col$treated,
    `Therapeutic cluster` = annotation_col$therapeutic_cluster,
    col = list(
        Study = study_colors[unique(annotation_col$study)],
        `Sample type` = pm_colors,
        Treated = treatment_colors,
        `Therapeutic cluster` = sctherapy_colors
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontface = "bold", fontsize = 12)
)

# 4. Row annotations
annotation_row <- drug_info %>%
    rownames_to_column("IDs") %>%
    select(IDs, drug_name, collapsed.MoAs) %>%
    mutate(drug_name = make.unique(drug_name, sep = "_")) %>%
    select(drug_name, collapsed.MoAs) %>%
    column_to_rownames("drug_name")

missing_drugs <- setdiff(rownames(gsva_scaled), rownames(annotation_row))
annotation_row <- rbind(
    annotation_row,
    data.frame(
        collapsed.MoAs = rep("Immunotherapy", 4),
        row.names = missing_drugs
    )
)[rownames(gsva_scaled), , drop = FALSE]

row_ha <- rowAnnotation(
    MoA = annotation_row[rownames(gsva_scaled), "collapsed.MoAs"],
    col = list(MoA = MoAs_colors),
    annotation_name_side = "top",
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontface = "bold", fontsize = 12),
    annotation_name_rot = 45,
    annotation_legend_param = list(
        MoA = list(direction = "vertical")
    )
)

# 5. Create and save heatmap
ht <- Heatmap(
    gsva_scaled,
    name = "Scaled\nGSVA",
    col = colorRamp2(
        seq(-3, 3, length.out = 100),
        colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
    ),
    column_order = rownames(annotation_col[order(annotation_col$therapeutic_cluster),]),
    top_annotation = col_ha,
    right_annotation = row_ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10),
    row_title = "Drugs",
    column_title = "Subclones",
    heatmap_legend_param = list(
        title = "Scaled\nGSVA",
        direction = "horizontal"
    ),
    show_heatmap_legend = TRUE
)

png(
    "subclone_level/plots/drug_signatures_subclones_heatmap.png",
    width = 15, 
    height = 12, 
    units = "in", 
    res = 300)
    
draw(ht, 
    heatmap_legend_side = "bottom",
    annotation_legend_side = "right",
    ht_gap = unit(2, "cm"),
    padding = unit(c(2, 10, 2, 20), "mm")
    )
dev.off()