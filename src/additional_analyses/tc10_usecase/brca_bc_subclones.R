library(Seurat)
library(BPCells)
library(GSVA)
library(GSEABase)
library(edgeR)
library(dplyr)
library(tidyverse)
library(uwot)
setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/brca_usecase")
set.seed(123)

# Load Seurat object with BRCA patient cells
brca <- readRDS("seu_brca_harmony.rds")

# Select only BRCA cells assigned to therapeutic clusters (TCs)
brca <- subset(brca, subset = !is.na(scTherapy_cluster))

# Pseudobulk at subclone level
brca_bulk <- AggregateExpression(
  object = brca,
  slot = "counts",
  return.seurat = T,
  group.by = c("scevan_subclone")
)

# Compute Beyondcell enrichment in drug sensitivity signatures
mat <- brca_bulk[["RNA"]]$counts
mat <- as.matrix(mat)

write.table(x = mat, file = "subclones_expr_mat.tsv", sep = "\t")

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
gsets <- GSEABase::getGmt(con = "../../reference/bc_immuno.gmt")

## gsva parameters
gsvapar <- gsvaParam(
  exprData = mat_cpm,
  geneSets = gsets,
  kcdf = "Gaussian",
  maxDiff = TRUE
)

gsva_bc <- gsva(gsvapar)

write.table(gsva_bc, "subclones_bc_gsva.tsv", sep ="\t")

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

write.table(final_matrix, "subclones_bc_gsva_aggregated.tsv", sep ="\t")


# Compute UMAP of subclones based on drug enrichment
gsva_final <- read.table("subclones_bc_gsva_aggregated.tsv", sep ="\t", header = TRUE)

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
        brca@meta.data %>%
            mutate(subclone = gsub("[^a-zA-Z0-9]", "", scevan_subclone), 
                   therapeutic_cluster = factor(scTherapy_cluster, levels = as.character(1:10))) %>%
            select(subclone, study, sample_type, treated, therapeutic_cluster) %>%
            distinct(),
        by = "subclone"
    ) %>%
    column_to_rownames(var = "subclone")

source(("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R"))
tcs_umap <- ggplot(
    data = umap_transform,
    aes(x = V1, y = V2, color = therapeutic_cluster)) +
    geom_point(size = 2, alpha = 0.7) + 
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
    scale_color_manual(
        values = sctherapy_colors,
        name = "Therapeutic cluster"
    ) +
    scale_x_continuous(name = "UMAP1") +
    scale_y_continuous(name = "UMAP2")

ggsave(
    filename = "plots/tcs_uwot.png", 
    plot = tcs_umap,
    dpi = 300,
    height = 4,
    width = 5
)


## Check the enrichment of TC10 correlated drugs
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
drugs_of_interest <- c(
            "sig-20883", "sig-21356", "sig-20886", "sig-21369", "sig-21370",
            "sig-20922", "sig-20900", "sig-21212", "sig-20889", "sig-20921",
            "sig-20995", "sig-20879", "sig-21166", "sig-20911", "sig-21041",
            "sig-20898", "sig-20925", "sig-20926", "sig-20952", "sig-21405",
            "sig-20907", "sig-20932", "sig-21323", "sig-21324", "sig-20890",
            "sig-20941", "sig-20967", "sig-21141", "sig-20980", "sig-21250",
            "sig-21238", "sig-20958", "sig-21297", "sig-21120", "sig-20881",
            "sig-21144", "sig-20949", "sig-20924", "sig-21007", "sig-20904",
            "sig-21194", "sig-20896", "sig-21193", "sig-21039", "sig-20957",
            "sig-21038", "sig-20888", "sig-21182"
        )

# 1. Prepare drug information
drugs <- data.table::fread("../../reference/final_moas - Collapsed.tsv")
drug_info <- drugs %>%
    filter(IDs %in% drugs_of_interest) %>%
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
            TRUE ~ collapsed.MoAs
        )
    ) %>%
    column_to_rownames("IDs")

# 2. Prepare drug subset and scale gsva matrix
gsva_subset <- gsva_final[drugs_of_interest, , drop = FALSE]
rownames(gsva_subset) <- drug_info[rownames(gsva_subset), "drug_name"]
colnames(gsva_subset) <- gsub("[^a-zA-Z0-9]", "", colnames(gsva_subset))
gsva_scaled <- t(scale(t(gsva_subset)))

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
    select(drug_name, collapsed.MoAs) %>%
    column_to_rownames("drug_name")

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

png("plots/drug_signatures_subclones_heatmap.png", width = 15, height = 12, units = "in", res = 300)
draw(ht, 
    heatmap_legend_side = "bottom",
    annotation_legend_side = "right",
    ht_gap = unit(2, "cm"),
    padding = unit(c(2, 10, 2, 20), "mm")
    )
dev.off()


# Epirubicin GSVA score on the top of UMAP
# Seleccionar droga de interés
drug_of_interest <- "EPIRUBICIN_GDSC"

# Añadir GSVA score al umap_transform
umap_transform_drug <- umap_transform %>%
    rownames_to_column(var = "subclone") %>%
    mutate(
        drug_score = t(gsva_subset[drug_of_interest,])
    )

# Plot con gradiente de color
drug_umap <- ggplot(
    data = umap_transform_drug,
    aes(x = V1, y = V2, color = drug_score)) +
    geom_point(size = 3) + 
    scale_color_gradient2(
        low = "#1D61F2",
        mid = "#F7F7F7",
        high = "#DA0078",
        midpoint = 0,  # Centro en 0
        limits = c(-0.8, 0.8),
        na.value = "grey50",
        name = "GSVA\nscore"
    ) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        panel.border = element_blank(),
        text = element_text(family = "Arial")
    ) +
    scale_x_continuous(name = "UMAP1") +
    scale_y_continuous(name = "UMAP2") +
    labs(title = paste("GSVA score:", drug_of_interest))

ggsave(
    filename = paste0("plots/umap_", drug_of_interest, ".png"), 
    plot = drug_umap,
    dpi = 300,
    height = 5,
    width = 5
)
