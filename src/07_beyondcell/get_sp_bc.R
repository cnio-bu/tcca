library(Seurat)
library(beyondcell)
library(dplyr)
library(tidyverse)
library(circlize)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/beyondcell_immuno")

all_studies <- list.files(
    path = ".",
    pattern = "*.rds",
    full.names = TRUE
)

samples_lvl2 <- read.table(
    "../seurat/v5/clinical_metadata_v4_clean_new.tsv",
    header = TRUE,
    sep = "\t"
) %>%
    mutate(study_sample = gsub(" ", "", paste0(study, "_", sample))) %>%
    pull(study_sample) %>%
    unique()

# Total drugs
adrenalnb_rui_chong <- readRDS("adrenalnb_rui_chong.rds")
all_drugs <- unique(unlist(
    lapply(
    adrenalnb_rui_chong,
    function(x) names(x@switch.point)
)))


study_mats_sp <- list()

for (study in all_studies){
    study_name <- basename(study)
    study_name <- gsub(x = study_name, pattern = ".rds", replacement = "")
    samples_bc <- readRDS(study)
    sample_names <- unlist(lapply(samples_bc, function(x) unique(x@meta.data$sample)))
    names(samples_bc) <- gsub(" ", "", sample_names)
    names(samples_bc) <- paste0(study_name, "_", names(samples_bc))
    print(study_name)
    # Remove samples not in level 2 (i.e., samples with fewer than 100 malignant
    # cells or missing clinical annotations)
    samples_bc <- samples_bc[names(samples_bc) %in% samples_lvl2]
    
    samples_aligned <- lapply(names(samples_bc), function(sample_name) {
        sp <- samples_bc[[sample_name]]@switch.point
        # Create vector with all possible drugs
        empty_vec <- rep(NA, length(all_drugs))
        names(empty_vec) <- all_drugs
        # Complete with existing valuees
        empty_vec[names(sp)] <- sp
        empty_vec
    })
    names(samples_aligned) <- names(samples_bc)

    # Join SPs from all samples in a single matrix
    mat_sp <- do.call(cbind, samples_aligned)
    print(dim(mat_sp))
    # Add the SP mat to the list with SP mats from all studies
    study_mats_sp[[study_name]] <- mat_sp
}

full_sp_mat <- do.call(cbind, study_mats_sp)

write.table(
    full_sp_mat, "../../beyondcell_immuno/full_sp_mat.tsv",
    sep = "\t"
)


# Filter samples for SP analysis
# Remove cells with more than 10% of zeros in the beyondcell mat and cell lines
bc <- readRDS("../../beyondcell_immuno/beyondcell_pancancer_immuno.Rds")
bc_samples <- bc@meta.data %>%
    mutate(
        study_sample = paste0(study, "_", sample),
        study_sample_match = gsub("[_. \\-]", "", paste0(study, "_", sample))
    ) %>%
    filter(patient != "ccl") %>%
    select(study_sample, study_sample_match) %>%
    distinct()

colnames(full_sp_mat) <- gsub("[_. \\-]", "", colnames(full_sp_mat))
full_sp_mat_filtered <- full_sp_mat[, bc_samples$study_sample_match]

colnames(full_sp_mat_filtered) <- bc_samples$study_sample[match(colnames(full_sp_mat_filtered), bc_samples$study_sample_match)]

write.table(
    full_sp_mat_filtered, "../../beyondcell_immuno/full_sp_mat_filtered.tsv",
    sep = "\t"
)

## Heatmap of SP per sample with clinical annotations
library(ComplexHeatmap)
source(file = "~/bc-meta/figures/TCCA_palette.R")
full_sp_mat <- read.table(
    "../../beyondcell_immuno/full_sp_mat_filtered.tsv",
    sep = "\t",
    header = TRUE
)
colnames(full_sp_mat) <- gsub("[_. \\-]", "", colnames(full_sp_mat))
# Add clinical annotations
metadata <- read.table("../seurat/v5/tcca_metadata.tsv", sep = "\t", header = TRUE)
metadata <- metadata %>%
    mutate(study_sample =  gsub("[_. \\-]", "", paste0(study, "_", sample))) %>%
    select(
        study_sample,
        study,
        sample,
        patient,
        sex,
        age,
        refined_tumor_site,
        sample_type,
        treated,
        refined_tumor_type,
        tme_archetype,
        tme_archetype_group
    ) %>%
    distinct()

translat_human_sites <- c(
    "adrenal_gland" = "Adrenal gland",
    "bladder" = "Bladder",
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "breast" = "Breast",
    "colon" = "Colon",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "kidney" = "Kidney",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue"
)


clinical_features <- metadata %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(refined_tumor_type %in% c("ALL", "CLL", "LAML", "MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        sample_origin = ifelse(patient == "ccl", "Cell line", "Patient")
    )


sample_annot_df <- clinical_features %>%
    select(
        study_sample,
        sex,
        adult_pediatric,
        is_blood,
        sample_type,
        summarised_tumor_site,
        treated,
        refined_tumor_type,
       #sample_origin,
        tme_archetype
    ) %>%
    filter(study_sample %in% colnames(full_sp_mat)) %>%
    as.data.frame() 

sample_annot_df$summarised_tumor_site <- translat_human_sites[sample_annot_df$summarised_tumor_site]

rownames(sample_annot_df) <- sample_annot_df$study_sample
sample_annot_df$study_sample <- NULL
full_sp_mat <- full_sp_mat[, rownames(sample_annot_df)]

colnames(sample_annot_df) <- c(
    "Chromosomal sex",
    "Age group",
    "Solid/Liquid",
    "Sample type",
    "Sample site",
    "Treatment",
    "Cancer type",
    #"Sample origin",
    "TME archetype"
)

pals <- list(
    "Chromosomal sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Sample type" = pm_colors,
    "Sample site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    #"Sample origin" = patient_ccl_colors,
    "Cancer type" = tumor_type_colors,
    "TME archetype" = tme_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = sample_annot_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_legend = FALSE
)


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

drugs <- drugs[drugs$IDs %in% rownames(full_sp_mat), ]
rownames(drugs) <- drugs$IDs

drugs$collapsed.MoAs <- ifelse(
    drugs$collapsed.MoAs %in% names(MoAs_colors),
    drugs$collapsed.MoAs,
    "Other"
)
MoAs <- drugs[rownames(full_sp_mat), c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"
rownames(MoAs) <- rownames(full_sp_mat)

# Add mechanism of action to immunotherapies
MoAs$`Mechanism of action`[is.na(MoAs$`Mechanism of action`)] <- "Immunotherapy"
# Order drugs based on MoA
ord <- order(MoAs[rownames(full_sp_mat), "Mechanism of action"])
ordered_rownames <- rownames(full_sp_mat)[ord]
ordered_drug_names <- drugs[ordered_rownames, "preferred.drug.names"]
ordered_drug_names[is.na(ordered_drug_names)] <- ordered_rownames[is.na(ordered_drug_names)]

MoAs <- MoAs[ordered_rownames, , drop = FALSE]
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
    file = "../../beyondcell_immuno/results/heatmap_sp_per_sample_filtered.png",
    res = 500,
    width = 19,
    height = 19,
    units = "in"
)


mat <- as.matrix(full_sp_mat)
mat <- scale(x = mat, center = TRUE, scale = TRUE)
apply(head(apply(mat, 2, function(x) is.na(x))), 2, sum)[apply(head(apply(mat, 2, function(x) is.na(x))), 2, sum) != 0]
heat <- ComplexHeatmap::Heatmap(
    mat = mat[ordered_rownames,],
    #col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    cluster_rows = FALSE,
    cluster_row_slices = TRUE,
    row_split = MoAs$`Mechanism of action`,
    column_order = rownames(sample_annot_df)[order(sample_annot_df$`Cancer type`)],
    cluster_columns = FALSE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split = sample_annot_df$`Cancer type`,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    #row_title = unique(sort(MoAs[rownames(full_sp_mat), "Mechanism of action"])),
    row_title = NULL,
    row_order = ordered_rownames,
    row_labels  = ordered_drug_names,
    show_row_names = TRUE,
    column_names_side = "top",
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 2),
    column_title = NULL,
    heatmap_legend_param = list(title = "Switch Point"),
    heatmap_width = unit(14, "in"),
    heatmap_height = unit(17, "in")
)

ht_opt(
    "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"), "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
    "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top")
dev.off()


# Trastuzumab response in primary patient samples of breast cancer
brca <- metadata %>%
    filter(refined_tumor_type == "BRCA" & sample_type == "p" & patient != "ccl") %>%
    mutate(study_sample = gsub("[_. \\-]", "", paste0(study, "_", sample))) %>%
    select(study_sample, tme_archetype_group)

brca$sp <- t(full_sp_mat["LIU_TRASTUZUMAB_RESPONSE", brca$study_sample])
brca <- brca %>%
    mutate(sp_group = case_when(
        sp < 0.5 ~ "Responder",
        sp >= 0.5 ~ "Non-responder"
    ))



barplot <- ggplot(
    brca,
    aes(x = sp_group, fill = tme_archetype_group)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = tme_group_colors) +
    labs(x = "Predicted response (SP)", y = "Proportion of samples", fill = "TME archeatype") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 16, margin = margin(t = 6), face = "bold"),
        axis.title.y = element_text(size = 16, margin = margin(r = 6), face = "bold"),
        axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)
    ) +
    guides(fill = guide_legend(ncol = 2))
ggsave(
    filename = "../../beyondcell_immuno/results/barplot_trastuzumab_response_brca.png",
    plot = barplot,
    width = 7, height = 6)

# Plot SPs of Wang Melanoma Hot-vs-Cold signature per TME archetype group
sp_hot_vs_cold <- as.data.frame(t(full_sp_mat["WANG_MELANOMA_HOT_VS_COLD", ]))

sp_hot_vs_cold <- sp_hot_vs_cold %>%
    rownames_to_column(var = "study_sample") %>%
    left_join(
        clinical_features %>%
            select(study_sample, tme_archetype_group),
        by = "study_sample"
    ) %>%
    mutate(tme_archetype_group = factor(tme_archetype_group,
        levels = c(
            "Tcell_centric",
            "Myeloid_centric",
            "Immune_rich",
            "Immune_stromal",
            "Immune_desert",
            "none"
        )
    )) %>%
    filter(tme_archetype_group != "none") %>%
    # Compute the inverse of SP (1 - SP) to improve interpretability
    mutate(sensitivity = 1 - WANG_MELANOMA_HOT_VS_COLD)


boxplot <- ggplot(
    sp_hot_vs_cold,
    aes(
        x = tme_archetype_group,
        y = sensitivity,
        fill = tme_archetype_group
    )
) +
    geom_boxplot() +
    scale_fill_manual(values = tme_group_colors) +
    labs(
        x = "TME archetype group",
        y = "Predicted response to ICB (anti-PD1 & anti-CTLA4)"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave(
    "../../beyondcell_immuno/results/boxplot_wang_melanoma_hot_vs_cold.pdf",
    plot = boxplot, 
    width = 8, 
    height = 4
)

# Statistical test
sp_hot_vs_cold <- sp_hot_vs_cold %>%
    filter(!is.na(sensitivity)) %>%
    mutate(group = ifelse(tme_archetype_group %in% c("Tcell_centric", "Myeloid_centric", "Immune_rich"), "Immune_hot", "Immune_cold"))

wilcox.test(sensitivity ~ group, data = sp_hot_vs_cold)