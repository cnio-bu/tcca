library(ComplexHeatmap)
library(circlize)
library(clustree)
library(edgeR)
library(factoextra)
library(fabia)
library(GSVA)
library(NbClust)
library(tidyverse)

source(file = "src/figures/TCCA_palette.R")


## function fix for factoextra
my_fviz_nbclust <- function(x, print.summary = TRUE, barfill = "steelblue", barcolor = "steelblue"){
    best_nc <- x$Best.nc
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    
    ss <- summary(best_nc$Number_clusters)
    cat("Among all indices: \n===================\n")
    for (i in 1:length(ss)) {
        cat("*", ss[i], "proposed ", names(ss)[i], "as the best number of clusters\n")
    }
    cat("\nConclusion\n=========================\n")
    cat("* According to the majority rule, the best number of clusters is ", 
        names(which.max(ss)), ".\n\n")
    
    df <- data.frame(Number_clusters = names(ss), freq = ss, 
                     stringsAsFactors = TRUE)
    p <- ggpubr::ggbarplot(df, x = "Number_clusters", y = "freq", 
                           fill = "steelblue", color = "steelblue") +
        ggplot2::labs(x = "Number of clusters k", 
                      y = "Frequency among all indices",
                      title = paste0("Optimal number of clusters - k = ", 
                                     names(which.max(ss))))
    p
}

mat <- read.table(
    file = "results/functional/pancancer_pseudobulk.tsv",
    check.names = FALSE
    )

mat <- as.matrix(mat)
#mat <- mat[, !grepl("cell.lines.gabriella.kinker", x = colnames(mat))]
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
    # 
    # if(grepl(pattern = "MP", x = gset@setName)){
    #     gsets_fix <- c(gsets_fix, gset)
    # }
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

gsva_es_centered <- scale(x = t(gsva_es), center = TRUE, scale = TRUE)

write.table(x = gsva_es, file = "results/functional/pseudo_bulk_gsva.tsv")

## Get hidden factors optimal
feature_cors <- corrplot::corrplot(
    cor(gsva_es_centered),
    type = "upper",
    method = "ellipse",
    tl.cex = 0.9
    )

res.pca <- prcomp(t(gsva_es_centered))

# Visualize eigenvalues/variances
pca_screen <- fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
ggsave(
    plot = pca_screen,
    filename = "results/functional/screeplot_clustering_pseudobulk.png",
    dpi = 300
    )

set.seed(120394)

# function to compute total within-cluster sum of squares
elbow_scree <- fviz_nbclust(
    t(gsva_es_centered),
    FUN = hcut,
    k.max = 24,
    method = "wss"
    ) + 
    theme_minimal()  + 
    ggtitle("Optimal elbow for functional matrix")

ggsave(
    plot = elbow_scree,
    filename = "results/functional/elbow_clustering_pseudobulk.png",
    dpi = 300
)

# nbclust
res.nbclust <- NbClust(t(gsva_es), 
                       distance = "euclidean",
                       min.nc = 2, 
                       max.nc = 24, 
                       method = "complete", 
                       index ="all"
                       )


nbclust_plot <- my_fviz_nbclust(res.nbclust) + 
    theme_minimal() + 
    ggtitle("NbClust's optimal number of clusters")

ggsave(
    plot = nbclust_plot,
    filename = "results/functional/consensus_index_clusters.png",
    dpi = 300
    )

## clustree now
tmp <- NULL
for (k in 1:10){
    tmp[k] <- kmeans(scale(t(gsva_es)), k, nstart = 30)
}
df <- data.frame(tmp)
# add a prefix to the column names
colnames(df) <- seq(1:10)
colnames(df) <- paste0("k",colnames(df))
# get individual PCA
df.pca <- prcomp(df, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df <- bind_cols(as.data.frame(df), as.data.frame(ind.coord))
clustree(df, prefix = "k")


n_hidden_factors <- 7

fabias <- fabia::fabias(
    X = t(gsva_es_centered),
    p = n_hidden_factors,
    center = 0,
    norm = 0,
    nL = 1,
    non_negative = 1
    )

res <- fabia::extractBic(fact = fabias)

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
            refined_tumor_site %in% names(translat_human_sites) ~ translat_human_sites[refined_tumor_site],
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


pals = list(
    "Age" = age_colors,
    "Chromosomal sex" = sex_colors,
    "Solid/Liquid" = sl_colors,
    "Primary or metastasis" = pm_colors,
    "Tumor site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "Therapeutic Cluster" = tcs_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    "Age" = clinical_features$adult_pediatric,
    "Chromosomal sex" = clinical_features$sex,
    "Treatment" = clinical_features$treated,
    "Primary or metastasis" = clinical_features$sample_type,
    "Tumor site" = clinical_features$summarised_tumor_site,
    col = pals,
    which = "column"
)

## Use TCCA color palette

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
    cluster_column_slices = TRUE,
    cluster_row_slices = TRUE,
    row_split = bicluster_table$bicluster,
    row_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    column_split = sample_table$bicluster,
    column_order = sample_table[order(sample_table$bicluster), ]$sample_study,
    row_names_gp = gpar(fontsize = 7),
    top_annotation = top_annotation,
    row_labels = row_labels,
    heatmap_legend_param = list(
        title = "Enrichment score",
        direction = "horizontal",
        title_position = "lefttop"
        )
    )

svg(
    filename = "results/figures/heatmap_all_samples_all_gsets.svg",
    width = 14,
    height = 14
)
draw(bulk_heat)
dev.off()

png(
    filename = "results/figures/pseudobulk_heatmap.png",
    res = 300,
    bg = "white",
    units = "in",
    width = 14,
    height = 14
    )

draw(bulk_heat, heatmap_legend_side = "bottom")
dev.off()

### test

## bulk heat by fabia biclust
bulk_heat <- ComplexHeatmap::Heatmap(
    matrix = gsva_es,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    cluster_column_slices = TRUE,
    cluster_row_slices = TRUE,
    row_split = c(1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,5,5,6,7,8,9,9,10,11,11,11,11,11,12,12,12,13,13,13,13,13,14,14,14,14),
#    row_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    column_split = 7,
#    column_order = sample_table[order(sample_table$bicluster), ]$sample_study,
    row_names_gp = gpar(fontsize = 7),
    top_annotation = top_annotation,
#    row_labels = row_labels,
    heatmap_legend_param = list(
        title = "Enrichment score",
        direction = "horizontal",
        title_position = "lefttop"
    )
)


