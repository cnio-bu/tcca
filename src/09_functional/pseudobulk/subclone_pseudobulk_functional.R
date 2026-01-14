library(Seurat)
library(BPCells)
library(GSVA)
library(GSEABase)
library(edgeR)
library(dplyr)
library(tidyverse)
library(factoextra)
library(NbClust)
library(clustree)
library(ComplexHeatmap)

setwd("/storage/scratch01/shared/projects/bc-meta/")
set.seed(123)

# Load Seurat object with
seu.lvl2 <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
metadata <- read.table(
  "single_cell/seurat/v5/tcca_metadata.tsv",
  sep = "\t",
  header = TRUE,
  row.names = NULL
) %>%
  column_to_rownames(var = "cell") %>%
  mutate(study_sample = paste0(study, "_", sample))

seu.lvl2@meta.data <- metadata

seu.lvl2 <- JoinLayers(seu.lvl2)

# Select only TCCA cells assigned to therapeutic clusters (TCs)
malignant_subset <- subset(
  seu.lvl2,
  subset = malignancy == "True" &
   !(scevan_subclone %in% c("","non_tumor"))
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

# Read the GMT file with published gene signatures (MSigDB, Gavish et al, Barkley et al, etc)
gsets <- GSEABase::getGmt(con = "reference/combined_gsets_functional.gmt")
gsets_fix <- list()
for(gset in gsets){
    gset@geneIds <- gset@geneIds[gset@geneIds != ""]
    gset@geneIds <- unlist(strsplit(gset@geneIds, split = "\\s+"))

    gsets_fix <- c(gsets_fix, gset)
}

gsets_fixed <- GSEABase::GeneSetCollection(gsets_fix)

# Add our own MPs
mp_list <- readRDS("functional_nmf/sample_wise/metaprograms_cpm/mp_list_reordered.rds")
names(mp_list) <- paste0("MP_", 1:length(mp_list), "_TCCA_UP")

# Join all functional signatures
combined_gsets <- GeneSetCollection(
    c(gsets_fix, 
      lapply(names(mp_list), function(mp_name) {
          GeneSet(
              geneIds = mp_list[[mp_name]],
              setName = mp_name,
              shortDescription = paste("Metaprogram", mp_name)
          )
      })
    )
)


## gsva parameters
gsvapar <- gsvaParam(
  exprData = mat_cpm,
  geneSets = combined_gsets,
  kcdf = "Gaussian",
  maxDiff = TRUE
)

gsva_fc <- gsva(gsvapar)
colnames(gsva_fc) <- gsub("-", "_", colnames(gsva_fc))

# For signatures with UP and DOWN gene sets compute the substraction of GSVA scores
signatures_UP <- rownames(gsva_fc)[grepl("_UP$", rownames(gsva_fc))]
signatures_DN <- rownames(gsva_fc)[grepl("_DOWN|DN$", rownames(gsva_fc))]
base_UP <- sub("_UP$", "", signatures_UP)
base_DN <- sub("_DOWN|_DN$", "", signatures_DN)
bidirectional_sigs <- intersect(base_UP, base_DN)

final_mat <- gsva_fc

for (base_name in bidirectional_sigs) {
  up_name <- paste0(base_name, "_UP")
  dn_name <- paste0(base_name, "_DOWN")
  dn_name <- ifelse(dn_name %in% rownames(gsva_fc), 
                    dn_name, 
                    paste0(base_name, "_DN")
                    )

  # Compute difference
  final_mat <- rbind(final_mat, gsva_fc[up_name, ] - gsva_fc[dn_name, ])
  rownames(final_mat)[nrow(final_mat)] <- base_name
  
  # Remove UP and DOWN drug signature enrichment
  final_mat <- final_mat[!rownames(final_mat) %in% c(up_name, dn_name), ]
}

write.table(final_mat, "functional/subclones_pseudobulk_gsva.tsv", sep ="\t")




#### Move to local running of biclustering and visualization ######
setwd("/Users/mariagb/OneDrive-CNIO/2nd_year/bc-meta/functional_pseudobulk/subclone")

# Santi's function to extract optimal number of clusters from NbClust object
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

# Run biclustering
final_mat <- read.table( "subclones_pseudobulk_gsva.tsv", sep = "\t", header = TRUE)
head(final_mat)
final_mat_centered <- scale(x = t(final_mat), center = TRUE, scale = TRUE)

# Compute correlation between functional signatures
png("correlation_pathways.png", units = "in", width = 12, height = 12, res = 300)
feature_cors <- corrplot::corrplot(
    cor(final_mat_centered),
    type = "upper",
    method = "ellipse",
    tl.cex = 0.8,
    tl.col = "black"
)
dev.off()

# Perform PCA on the functional signature matrix
res.pca <- prcomp(t(final_mat_centered))
## Visualize eigenvalues/variances
pca_screen <- fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
ggsave(
    plot = pca_screen,
    filename = "screeplot_clustering_pseudobulk.png",
    dpi = 300
    )

set.seed(120394)


# function to compute total within-cluster sum of squares
elbow_screen <- fviz_nbclust(
    t(final_mat_centered),
    FUN = hcut,
    k.max = 24,
    method = "wss"
    ) + 
    theme_minimal()  + 
    ggtitle("Optimal elbow for functional matrix")

ggsave(
    plot = elbow_screen,
    filename = "elbow_clustering_pseudobulk.png",
    dpi = 300
)


# nbclust
res.nbclust <- NbClust(t(final_mat), 
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
    filename = "consensus_index_clusters.png",
    dpi = 300
    )

# Clustree
tmp <- list()
for (k in 1:10){
    tmp[[k]] <- kmeans(scale(t(final_mat)), k, nstart = 30)
}
df <- do.call(cbind, lapply(tmp, function(x) x$cluster))

# add a prefix to the column names
colnames(df) <- seq(1:10)
colnames(df) <- paste0("k",colnames(df))
# get individual PCA
df.pca <- prcomp(df, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df <- bind_cols(as.data.frame(df), as.data.frame(ind.coord))
clustree(df, prefix = "k")



# Run biclustering with FABIA
n_hidden_factors <- 8

fabias <- fabia::fabias(
    X = t(final_mat_centered),
    p = n_hidden_factors,
    center = 0,
    norm = 0,
    nL = 1,
    non_negative = 1
    )

res <- fabia::extractBic(fact = fabias)


# Extract subclones and signatures per bicluster
all_biclusters <- list()
all_subclones <- list()
for(i in c(1:n_hidden_factors)){
    this_biclust <- res$bic[i, ]
    ## check length of the bicluster rowwise (aka drugs)
    if(length(this_biclust$bixn) <= 0) {
        next
    }else {
        named_list_sigs <- this_biclust$bixn
        sig_contributions <- this_biclust$bixv
        named_list_subclones <- this_biclust$biypn
        subclone_contributions <- this_biclust$biypv
        names(sig_contributions) <- named_list_sigs
        names(subclone_contributions) <- named_list_subclones
        all_biclusters[[i]] <- sig_contributions
        all_subclones[[i]] <- subclone_contributions
    }
}


bicluster_table <- enframe(all_biclusters) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "signature" = value_id
    )

## Try to remove signatures whose sign is opposite from the bicluster consensus
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



## subclone table
subclone_table <- enframe(all_subclones) %>%
    unnest_longer(col = "value") %>%
    rename(
        "bicluster" = name,
        "cluster_contribution" = value,
        "subclone" = value_id
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

subclone_table <- subclone_table %>%
    mutate(
        abs_contribution = abs(cluster_contribution)
    ) %>%
    group_by(subclone) %>%
    arrange(desc(abs_contribution)) 

subclone_table <- subclone_table[!duplicated(subclone_table$subclone), ]
subclone_table <- subclone_table %>%
    mutate(subclone = gsub("[^A-Za-z0-9]", "", subclone))


# Heatmap of biclusters
tcca_metadata <- read.table("../../cohort_statistics/tcca_metadata.tsv", 
                                sep = "\t", header = TRUE)

translat_human_sites <- c(
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "adrenal_gland" = "Adrenal gland",
    "breast" = "Breast",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue",
    "bladder" = "Bladder",
    "colon" = "Colon",
    "kidney" = "Kidney"
)

subclone_metadata <- tcca_metadata %>%
    filter(malignancy == "True" &
            !(scevan_subclone %in% c("","non_tumor"))
        ) %>%
    mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "Other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(refined_tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated != "", ifelse(treated == "t", "Treated", "Untreated"), NA),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        therapeutic_cluster = ifelse(is.na(scTherapy_cluster), "Unclassified", scTherapy_cluster)
    ) %>%
    mutate(scevan_subclone = gsub("[^A-Za-z0-9]", "", scevan_subclone)) %>%
    filter(scevan_subclone %in% subclone_table$subclone)

clinical_annot <- subclone_metadata %>%
    select(
        scevan_subclone,
        sex,
        adult_pediatric,
        is_blood,
        sample_type,
        summarised_tumor_site,
        treated,
        refined_tumor_type,
        tme_archetype,
        therapeutic_cluster) %>%
    as.data.frame() %>%
    distinct()
rownames(clinical_annot) <- NULL
clinical_annot <- clinical_annot %>% column_to_rownames(var = "scevan_subclone")
clinical_annot$summarised_tumor_site <-  translat_human_sites[clinical_annot$summarised_tumor_site]
clinical_annot <- clinical_annot[subclone_table$subclone, ]

colnames(clinical_annot) <- c(
    "Chromosomal sex",
    "Age group",
    "Solid/Liquid",
    "Sample type",
    "Sample site",
    "Treatment",
    "Cancer type",
    "TME archetype",
    "Therapeutic cluster"
)

source("/Users/mariagb/Documents/tcca/src/12_figures/TCCA_palette.R")
sctherapy_colors <- c(sctherapy_colors, "Unclassified" = "#BBB9B7")
pals <- list(
    "Chromosomal sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Sample type" = pm_colors,
    "Sample site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "Cancer type" = tumor_type_colors,
    "TME archetype" = tme_colors,
    "Therapeutic cluster" = sctherapy_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df =  clinical_annot,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  show_legend = FALSE
)

colnames(final_mat) <- gsub("[^A-Za-z0-9]", "", colnames(final_mat))
bulk_heat <- ComplexHeatmap::Heatmap(
    matrix = final_mat[bicluster_table$signature, subclone_table$subclone],
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    cluster_column_slices = TRUE,
    cluster_row_slices = TRUE,
    row_split = bicluster_table$bicluster,
    row_order = bicluster_table[order(bicluster_table$bicluster), ]$signature,
    column_split = subclone_table$bicluster,
    column_order = subclone_table[order(subclone_table$bicluster), ]$subclone,
    row_labels = gsub("_UP", "", bicluster_table$signature),
    row_names_gp = gpar(fontsize = 7),
    top_annotation = top_annotation,
    heatmap_legend_param = list(
        title = "Enrichment score",
        direction = "horizontal",
        title_position = "lefttop"
        )
    )
png("bicluster_heatmap_tcs8.png", width = 12, height = 15, units = "in", res = 300)
draw(bulk_heat, heatmap_legend_side = "top")
dev.off()

pdf("bicluster_heatmap_tcs8.pdf", width = 12, height = 15)
draw(bulk_heat, heatmap_legend_side = "top")
dev.off()
