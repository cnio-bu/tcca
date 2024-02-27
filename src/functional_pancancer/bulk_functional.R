library(ComplexHeatmap)
library(edgeR)
library(fabia)
library(GSVA)
library(tidyverse)

mat <- read_tsv(
    file = "results/functional/pancancer_pseudobulk.tsv",
    ) %>%
    as.data.frame()

rownames(mat) <- mat$gene
mat$gene <- NULL

mat <- as.matrix(mat)

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

## gsva
gsets <- GSEABase::getGmt(con = "reference/combined_gsets_functional.gmt")

gsets_fix <- list()
for(gset in gsets){
    gset@geneIds <- gset@geneIds[gset@geneIds != ""]
    gsets_fix <- c(gsets_fix, gset)
}

gsets_fixed <- GSEABase::GeneSetCollection(gsets_fix)

## gsva parameters
gsvapar <- gsvaParam(
    exprData = mat_cpm,
    geneSets = gsets_fixed,
    kcdf = "Gaussian",
    maxDiff = TRUE
    )

gsva_es <- gsva(gsvapar)

## fix gsva sample names
colnames(gsva_es) <- gsub(pattern = "-", replacement = "_", x = colnames(gsva_es))
colnames(gsva_es) <- gsub(pattern = "\\.", replacement = "_", x = colnames(gsva_es))
colnames(gsva_es)[colnames(gsva_es) == "T19_1_adrenalnb_rui_chong"] <- "T19_adrenalnb_rui_chong"

n_hidden_factors = 20

fabias <- fabia::fabias(
    X = gsva_es,
    p = n_hidden_factors,
    center = 2,
    norm = 0,
    nL = 1,
    non_negative = 1
    )

res <- fabia::extractBic(fact = fabias, thresZ = 0.7)

all_biclusters <- list()
all_samples <- list()
for(i in c(1:n_hidden_factors)){
    this_biclust <- res$bic[i, ]
    ## check length of the bicluster rowwise (aka drugs)
    if(length(this_biclust$bixn) <= 0) {
        next
    }else {
        named_list_sigs <- this_biclust$bixn
        sig_contributions <- this_biclust$bixv
        named_list_samples <- this_biclust$biypn
        sample_contributions <- this_biclust$biypv
        names(sig_contributions) <- named_list_sigs
        names(sample_contributions) <- named_list_samples
        all_biclusters[[i]] <- sig_contributions
        all_samples[[i]] <- sample_contributions
    }
}

bicluster_table <- enframe(all_biclusters) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "signature" = value_id
    )

## Try to remove drugs whose sign is opposite from the bicluster consensus
bicluster_table <- bicluster_table %>%
    mutate(
        sig_dir = sign(cluster_contribution)
    ) %>%
    group_by(bicluster) %>%
    mutate(
        dom_sig =  case_when(
            sum(sig_dir) >= 0 ~ 1,
            TRUE ~ -1
        )
    ) %>%
    filter(
        dom_sig == sig_dir
    )


## sample table
sample_table <- enframe(all_samples) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "sample" = value_id
    ) %>%
    mutate(
        sig_dir = sign(cluster_contribution),
        dom_sig =  case_when(
            sum(sig_dir) >= 0 ~ 1,
            TRUE ~ -1
            )
    ) %>%
    filter(
        dom_sig == sig_dir
    )

sample_table <- sample_table %>%
    mutate(
        abs_contribution = abs(cluster_contribution)
    ) %>%
    group_by(sample) %>%
    arrange(desc(abs_contribution)) 

sample_table <- sample_table[!duplicated(sample_table$sample), ]

sample_table <- sample_table %>%
    mutate(
        sample_study = sample
    ) %>%
    mutate(
        sample_study = str_replace_all(string = sample_study, pattern = "\\.", replacement = "_"),
        sample_study = str_replace_all(string = sample_study, pattern = "-", replacement = "_"),
        sample_study = str_replace_all(string = sample_study, pattern = "\\.", replacement = "_")
    )

## load clinical database
clinical_metadata <- data.table::fread(
    "results/annotation/clinical_metadata_v4_clean.tsv"
    )


## human readable origins
translat_human_sites <- c(
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "adrenal_gland" = "Adrenal gland",
    "breast" = "Breast",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue"
)

clinical_features <- clinical_metadata %>%
    mutate(
        metagroup = case_when(
            treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_untreated",
            treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_untreated",
            treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_treated",
            treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_treated",
            study == "cell_lines_gabriella_kinker" ~ "cell_line",
            TRUE ~ "other"
        )
    ) %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary")
    ) %>%
    mutate(
        sample_functional_mat = paste(sample, study, sep = "_"),
        sample_functional_mat = gsub(pattern = "-", replacement = "_", x = sample_functional_mat),
        sample_functional_mat = gsub(pattern = "\\.", replacement = "_", x = sample_functional_mat)
    ) %>%
    filter(
        sample_functional_mat %in% unique(sample_table$sample_study)
    )


top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Primary or metastasis" = clinical_features$sample_type,
    "Age" = clinical_features$adult_pediatric,
    "Treatment" = clinical_features$treated,
    which = "column"
)


row_labels <- rownames(gsva_es[bicluster_table$signature, ])
row_labels <- gsub(pattern = "_UP", replacement = "", x = row_labels)

row_labels <- stringr::str_to_title(string = row_labels)

## Detect hallmarks and edit them
row_labels[
    stringr::str_starts(pattern = "Hallmark_", string = row_labels)
    ] <- paste0(row_labels[
        stringr::str_starts(pattern = "Hallmark_", string = row_labels)
    ], " (Hallmark)")

row_labels <- stringr::str_remove(string = row_labels, pattern = "^Hallmark_")

# Remove CancerSEA labels and append them to the end
row_labels[
    stringr::str_ends(pattern = "cancersea", string = row_labels)
] <- paste0(row_labels[
    stringr::str_ends(pattern = "cancersea", string = row_labels)
], " (CancerSEA)")
 
row_labels <- stringr::str_remove(string = row_labels, pattern = "_cancersea")
row_labels <- stringr::str_replace_all(string = row_labels, pattern = "_", replacement = " ")
row_labels <- stringr::str_replace(string = row_labels, pattern = "^Mp", replacement = "MP")

## Final edits can be performed in post-processing
   
## bulk heat by fabia biclust
bulk_heat <- ComplexHeatmap::Heatmap(
    matrix = gsva_es[bicluster_table$signature, sample_table$sample_study],
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    row_split = bicluster_table$bicluster,
    row_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    column_split = sample_table$bicluster,
    column_order = sample_table[order(sample_table$bicluster), ]$sample_study,
    row_names_gp = gpar(fontsize = 7),
    top_annotation = top_annotation,
    row_labels = row_labels
    )


