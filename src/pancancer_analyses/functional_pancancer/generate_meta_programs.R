# # ----------------------------------------------------------------------------
# # Load packages and seurat objects
# # ----------------------------------------------------------------------------
library(Seurat)
library(BPCells)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ComplexHeatmap)
source("~/bc-meta/figures/TCCA_palette.R")
source("~/bc-meta/functional_pancancer/robust_nmf_programs.R")


## Input:
# nmf_programs is a list in which each entry contains NMF gene-scores of a
# single sample and k value. In our study  # we ran NMF using ranks 4-9 on the
# top 7000 genes in each sample. Hence each entry in nmf_programs is a
# matrix with 7000 rows (genes) X k columns (NMF programs). We transform this
# list of matrices into a list where each element corrresponds to all the NMF
# gene-scores of a single sample (7000 genes x 39 programs).
# We define MPs in 2 steps:
# 1) The function robust_nmf_programs.R performs filtering, so that programs
# selected for defining MPs are:
#    i) Robust            - recur in more that one rank within the sample
#    ii) Non-redundant    - once a NMF program is selected, other programs
#                         within the sample that are similar to it are removed
#    iii) Not sample specific - has similarity to NMF programs in other samples
# ** Please see https://github.com/gabrielakinker/CCLE_heterogeneity for more
# details on how to define robust NMF programs
# 2) Selected NMFs are then clustered iteratively. At the end of the process,
# each cluster generates a list of the 50 genes (i.e. the MP) that represent
# the NMF programs that contributed to the cluster. Notably, not all initially
# selected NMFs end up participating in a cluster.


# Load NMF matrices
setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf/")

nmf_programs <- readRDS(
  "geneNMFprograms_allsamples.rds"
)

# Select programs from k = 4..10,15, 20, 25, 30
# pattern <- paste(paste0("k", c(4:10, 20)),collapse = "|")
# nmf_programs <- nmf_programs[grep(pattern, names(nmf_programs), value = TRUE)]

#-------------------------------------------------------------------------------
# Adapt RcppML::nmf output to define MPs
# ------------------------------------------------------------------------------
nmf_programs <- lapply(nmf_programs, function(sample_k) sample_k$w)
sample_names <- gsub("\\.k\\d+$", "", names(nmf_programs))
sample_split <- split(seq_along(nmf_programs), sample_names)

# Get a list where each entry is a matrix with 2000 rows (genes) x 39 columns
# (NMF programs)
sample_nmf <- lapply(names(sample_split), function(sample) {
  sample_programs <- nmf_programs[sample_split[[sample]]]
  names(sample_programs) <- paste0(sample, "_", "k", c(4:9))
  sample_programs <- lapply(names(sample_programs), function(n) {
    colnames(sample_programs[[n]]) <- paste(n,
                                            seq(1, ncol(sample_programs[[n]])),
                                            sep = ".")
    return(sample_programs[[n]])
  })
  sample_programs <- Reduce(f = cbind, x = sample_programs)
  return(sample_programs)
})

names(sample_nmf) <- names(sample_split)

#-------------------------------------------------------------------------------
# Select robust NMF programs using code from Gabriela Kinker et al, 2020
# ------------------------------------------------------------------------------

## Parameters
intra_min_parameter <- 35
intra_max_parameter <- 10
inter_min_parameter <- 10


# get top 50 genes for each NMF program
nmf_programs <- lapply(sample_nmf, function(x) {
  apply(x, 2, function(y) names(sort(y, decreasing = TRUE))[1:50])
})

# For each sample, select robust NMF programs (i.e. observed using different 
# ranks in the same sample), remove redundancy due to multiple ranks, and apply 
# a filter based on the similarity to programs from other samples.
nmf_filter_ccle <- robust_nmf_programs(
  nmf_programs,
  intra_min = intra_min_parameter,
  intra_max = intra_max_parameter,
  inter_filter = TRUE,
  inter_min = inter_min_parameter
)
print(paste("Robust_programs_selected!", format(Sys.time(), "%a %b %d %X %Y")))

nmf_programs <- lapply(nmf_programs, function(x) {
  x[, is.element(colnames(x), nmf_filter_ccle), drop = FALSE]
})
nmf_programs <- do.call(cbind, nmf_programs)


# Calculate similarity between programs
nmf_intersect <- apply(nmf_programs, 2, function(x) {
  apply(nmf_programs, 2, function(y) length(intersect(x, y)))
})

# Hierarchical clustering of the similarity matrix
nmf_intersect_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(
  as.dendrogram(nmf_intersect_hc),
  colMeans(nmf_intersect)
)
nmf_intersect <- nmf_intersect[
  order.dendrogram(nmf_intersect_hc),
  order.dendrogram(nmf_intersect_hc)
]

print(paste(
  "Similarity between programs computed by intersections",
  format(Sys.time(), "%a %b %d %X %Y")
))

saveRDS(nmf_programs, "metaprograms/nmf_programs_mat.rds")
saveRDS(nmf_intersect_hc, "metaprograms/nmf_intersect_hc.rds")
saveRDS(nmf_intersect, "metaprograms/nmf_intersect.rds")


nmf_programs <- readRDS("metaprograms/nmf_programs_mat.rds")
nmf_intersect_hc <- readRDS("metaprograms/nmf_intersect_hc.rds")
nmf_intersect <- readRDS("metaprograms/nmf_intersect.rds")


# ------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ------------------------------------------------------------------------------

### Parameters for clustering
# The minimal intersection cutoff for defining the first NMF program in a cluster
min_intersect_initial <- 15
# The minimal intersection cutoff for adding a new NMF to the forming cluster
min_intersect_cluster <- 15
# The minimal group size to consider for defining the first NMF program in a cluster
min_group_size <- 10

sorted_intersection <- sort(apply(nmf_intersect, 2, function(x) {
  (length(which(x >= min_intersect_initial)) - 1)
}), decreasing = TRUE)

# Every entry contains the NMFs of a chosen cluster
cluster_list <- list()
# Every entry contains the genes of the metaprogram
mp_list <- list()
# Every entry contains the genes of the metaprogram and their NMF scores
mp_nmf_scores <- list()

k <- 1
curr_cluster <- c()
nmf_intersect_original <- nmf_intersect

while (sorted_intersection[1] > min_group_size) {
  curr_cluster <- c(curr_cluster, names(sorted_intersection[1]))
  # Intersection between all remaining NMFs and genes in MP
  # Genes in the forming MP are first chosen to be those in the first NMF,
  # 'genes_mp' always has only 50 genes and evolves during cluster's formation
  genes_mp <- nmf_programs[, names(sorted_intersection[1])]
  # remove selected NMF
  nmf_programs <- nmf_programs[, -match(
    names(sorted_intersection[1]),
    colnames(nmf_programs)
  )]
  # Intersection between all other NMFs and 'genes_mp'
  intersection_with_genes_mp <- sort(apply(nmf_programs, 2, function(x) {
    length(intersect(genes_mp, x))
  }), decreasing = TRUE)
  # 'nmf_history' contains genes from all NMFs in the current cluster,
  # for redefining 'genes_mp' after adding a new NMF
  nmf_history <- genes_mp
  
  # Create gene list composed of intersecting genes (in descending order by
  # frequency). When the number of genes with a given frequency span beyond the
  # 50th genes, they are sorted according to their NMF score.
  while (intersection_with_genes_mp[1] >= min_intersect_cluster) {
    curr_cluster <- c(curr_cluster, names(intersection_with_genes_mp)[1])
    
    ## 'genes_mp' is newly defined each time according to all NMFs in the 
    # current cluster
    genes_mp_temp <- sort(
      table(c(
        nmf_history,
        nmf_programs[, names(intersection_with_genes_mp)[1]]
      )),
      decreasing = TRUE
    )
    # Genes with overlap equal to the 50th gene
    genes_at_border <- genes_mp_temp[which(genes_mp_temp == genes_mp_temp[50])]
    
    if (length(genes_at_border) > 1) {
      # Sort last genes in genes_at_border according to maximal NMF gene scores
      # Run across all NMF programs in curr_cluster and extract NMF scores for
      # each gene
      genes_curr_nmf_score <- c()
      
      for (i in curr_cluster) {
        curr_sample <- gsub("_k\\d\\.\\d", "", i)
        matched_indices <- match(names(genes_at_border),
                                 rownames(sample_nmf[[curr_sample]]))
        matched_genes <- matched_indices[!is.na(matched_indices)]
        q <- sample_nmf[[curr_sample]][matched_genes, i]
        # Sometimes when adding genes the names do not appear
        names(q) <- rownames(sample_nmf[[curr_sample]])[matched_genes]
        genes_curr_nmf_score <- c(genes_curr_nmf_score, q)
      }
      
      genes_curr_nmf_score_sort <- sort(genes_curr_nmf_score, decreasing = TRUE)
      genes_curr_nmf_score_sort <- genes_curr_nmf_score_sort[unique(
        names(genes_curr_nmf_score_sort)
      )]
      
      genes_mp_temp <- c(
        names(genes_mp_temp[which(genes_mp_temp > genes_mp_temp[50])]),
        names(genes_curr_nmf_score_sort)
      )
    } else {
      genes_mp_temp <- names(genes_mp_temp)[1:50]
    }
    nmf_history <- c(
      nmf_history,
      nmf_programs[, names(intersection_with_genes_mp)[1]]
    )
    genes_mp  <- genes_mp_temp[1:50]
    
    # Remove selected NMF
    nmf_programs <- nmf_programs[, -match(names(intersection_with_genes_mp)[1],
                                          colnames(nmf_programs))]
    # Intersection between all other NMFs and 'genes_mp'
    intersection_with_genes_mp <- sort(apply(nmf_programs, 2, function(x) {
      length(intersect(genes_mp, x))
    }), decreasing = TRUE)
  }
  
  # Store the final list of genes and their aggregated NMF scores for the cluster
  genes_mp_nmf_scores <- c()
  
  # Aggregate NMF scores for the final genes in the meta-program
  for (i in curr_cluster) {
    curr_sample <- gsub("_k\\d\\.\\d", "", i)
    matched_indices <- match(genes_mp, rownames(sample_nmf[[curr_sample]]))
    matched_genes <- matched_indices[!is.na(matched_indices)]
    q <- sample_nmf[[curr_sample]][matched_genes, i]
    names(q) <- rownames(sample_nmf[[curr_sample]])[matched_genes]
    genes_mp_nmf_scores <- c(genes_mp_nmf_scores, q)
  }
  
  # Aggregate scores to handle duplicate gene entries
  genes_mp_nmf_scores <- tapply(
    genes_mp_nmf_scores, 
    names(genes_mp_nmf_scores), 
    max
  )
  
  # Store the cluster, meta-program genes, and their scores
  cluster_list[[paste0("MP", k)]] <- curr_cluster
  mp_list[[paste0("MP", k)]] <- genes_mp
  mp_nmf_scores[[paste0("MP", k)]] <- genes_mp_nmf_scores[genes_mp]
  k <- k + 1
  
  # Remove current chosen cluster
  nmf_intersect <- nmf_intersect[
    -match(curr_cluster, rownames(nmf_intersect)),
    -match(curr_cluster, colnames(nmf_intersect))
  ]
  # Sort intersection of remaining NMFs not included in any of the previous 
  # clusters
  sorted_intersection <- sort(apply(nmf_intersect, 2, function(x) {
    (length(which(x >= min_intersect_initial)) - 1)
  }), decreasing = TRUE)
  
  curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}
print(paste("MPs computed!", format(Sys.time(), "%a %b %d %X %Y")))
saveRDS(cluster_list, "metaprograms/cluster_list.rds")
saveRDS(mp_list, "metaprograms/mp_list.rds")

cluster_list <- readRDS("metaprograms/cluster_list.rds")
mp_list <- readRDS("metaprograms/mp_list.rds")


# ------------------------------------------------------------------------------
# Heatmap of Jaccard similaity scores for NMF programs clustered into MPs
# ------------------------------------------------------------------------------

# Reorder MPs to group related ones together and remove MP11, linked only to brain metastases of 
# melanoma (with skin pigmentation enrichment)
mps <- c("MP1", "MP8", "MP14", "MP3", "MP13", "MP2", "MP9", "MP4", "MP10", "MP12", 
         "MP5", "MP6",  "MP7", "MP11")
cluster_list <- cluster_list[mps]
names(cluster_list) <- paste0("MP", 1:14)


#  Sort Jaccard similarity plot according to new clusters:
inds_sorted <- c()
for (j in seq_along(cluster_list)) {
  inds_sorted <- c(inds_sorted, match(
    cluster_list[[j]],
    colnames(nmf_intersect_original)
  ))
}

nmf_intersect_sort <- nmf_intersect_original[inds_sorted, inds_sorted]
rownames(nmf_intersect_sort) <- 1: ncol(nmf_intersect_sort)
colnames(nmf_intersect_sort) <- 1: ncol(nmf_intersect_sort)
labels <- ifelse(1:ncol(nmf_intersect_sort) %in% seq(0, 3000, by = 250), 
                 rownames(nmf_intersect_sort), "")

# Transform cluster list into a vector
mp_members <- unlist(lapply(seq_along(cluster_list), function(i) {
  group <- names(cluster_list)[i] # Get the name of the current group
  setNames(rep(group, length(cluster_list[[i]])), cluster_list[[i]])
}))

# Add tumor sample annotation to the heatmap
cluster_annotation <- data.frame(
  cluster = mp_members,
  stringsAsFactors = FALSE
)

programs_names <- rownames(cluster_annotation)
cluster_annotation$study_sample <- gsub(
  "_k\\d+\\.\\d+", "",
  rownames(cluster_annotation)
)
clinical <- data.table::fread(
  "../../clinical_metadata_v4_clean.tsv"
)
clinical$study_sample <- paste0(clinical$study, "__", clinical$sample)

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

# Add inferred sex
seu <- readRDS("../../seu_lvl2_sex_inferred.rds")
sex_inferred <- seu@meta.data %>%
  mutate(study_sample = paste0(study, "__", sample)) %>%
  select(sex, study_sample) %>%
  distinct()

cluster_annotation <- cluster_annotation %>%
  left_join(clinical, by = "study_sample") %>%
  left_join(sex_inferred, by = "study_sample") %>%
  mutate(
    summarised_tumor_site = case_when(
      refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
      TRUE ~ "other"
    ),
    adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
    is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML", "MM"),
                      "Liquid", "Solid"
    ),
    treated = ifelse(treated == "t", "Treated", "Untreated"),
    sex = ifelse(sex.y == "f", "Female", "Male"),
    sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
    cluster = factor(mp_members, levels = names(cluster_list))
  )

cluster_annotation$summarised_tumor_site <- translat_human_sites[
  cluster_annotation$summarised_tumor_site
]

cluster_annotation <- cluster_annotation %>%
  select(
    sex,
    adult_pediatric,
    is_blood,
    sample_type,
    summarised_tumor_site,
    treated,
    cluster
  )

rownames(cluster_annotation) <- 1: nrow(cluster_annotation)
colnames(cluster_annotation) <- c(
  "Sex",
  "Age group",
  "Solid/Liquid",
  "Sample type",
  "Sample site",
  "Treatment",
  "Meta-program"
)

# Reorder NMF programs to match the names in nmf_intersect
pals <- list(
  "Sex" = sex_colors,
  "Age group" = age_colors,
  "Solid/Liquid" = sl_colors,
  "Sample type" = pm_colors,
  "Sample site" = tumor_sites_colors,
  "Treatment" = treatment_colors
)


top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = cluster_annotation[names(pals)],
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                                 labels_gp = gpar(fontsize = 12),
                                 title_gap = unit(10, "mm")),
  show_legend = c(
    "Sample site" = FALSE
  )
)

# Customize legends for tumor site
tumor_site_legend <- Legend(
  at = names(pals$`Sample site`),
  legend_gp = gpar(fill = pals$`Sample site`),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 12),
  ncol = 3,  # Split Group legend into 2 columns
  gap = unit(5, "mm"),
  title = "Sample site"
)

# Create the color ramp
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(20), rev(magma(323, begin = 0.18)))
n_total <- length(custom_magma)
n_white <- 20
n_magma <- n_total - n_white
breaks <- c(
  seq(0, 4, length.out = n_white), 
  seq(4.0001, 25, length.out = n_magma)
  )

custom_col <- colorRamp2(breaks, custom_magma)

custom_col <- colorRamp2(seq(0, 25, length.out = length(custom_magma)), custom_magma)
heat <- ComplexHeatmap::Heatmap(
  mat = nmf_intersect_sort,
  col = custom_col,
  top_annotation = top_annotation,
  cluster_rows = FALSE,
  cluster_row_slices = TRUE,
  row_split = cluster_annotation$`Meta-program`,
  cluster_columns = FALSE,
  cluster_column_slices = TRUE,
  show_column_dend = FALSE,
  column_split = cluster_annotation$`Meta-program`,
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson",
  row_labels = labels,
  column_labels = labels,
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  row_names_side = "left",
  column_names_side = "bottom",
  row_title = "Robust NMF programs",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "Robust NMF programs",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title_side = "left",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    at = c(0, 5, 10, 15, 20, 25),
    title = "Similarity\n(Jaccard index)",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    title_gap = unit(10, "mm")
  ),
  heatmap_width = unit(10, "in"),
  heatmap_height = unit(10, "in")
)
png(
  file = "./metaprograms/heatmap_mps.png",
  res = 500,
  width = 16,
  height = 14,
  units = "in"
)
draw(heat, 
     annotation_legend_side = "top", 
     heatmap_legend_side = "right", 
     annotation_legend_list = list(tumor_site_legend))
dev.off()



# Annotate MPs
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(openxlsx)

mp_list <- mp_list[mps]
names(mp_list) <- paste0("MP", 1:14)
nmf_programs <- readRDS("metaprograms/nmf_programs_mat.rds")
universe <- rownames(sample_nmf[[1]])
msig_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
msig_list <- split(x = msig_df$gene_symbol, f = msig_df$gs_name)

func_annot <- lapply(mp_list, function(metaprogram) {
  fgRes <- fgsea::fora(pathways = msig_list,
                       genes = metaprogram,
                       universe = universe,
                       maxSize = 200)
  
  fgRes <- fgRes[fgRes$pval <= 0.05,]
  return(fgRes)
})


# Save results as a excel file
wb <- createWorkbook()
# Write each data frame to a separate sheet
for (i in seq_along(func_annot)) {
  addWorksheet(wb, sheetName = names(func_annot)[i]) # Use the list names as sheet names
  writeData(wb, sheet = names(func_annot)[i], x = func_annot[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, file = "fgsea_hallmarks_ora.xlsx", overwrite = TRUE)

# Check MPs from Kinker et al, Gavish et al., and Barkey et al and other functional gene sets
gmt_df <- read.gmt("../fork/bc-meta/reference/combined_gsets_functional.gmt")
gmt_list <- split(x = gmt_df$gene, f = gmt_df$term)
func_annot <- lapply(mp_list, function(program) {
  fgRes <- fgsea::fora(pathways = gmt_list,
                       genes = program,
                       universe = universe)
  
  fgRes <- fgRes[fgRes$pval <= 0.05,]
  return(fgRes)
})

# Save results as a excel file
wb <- createWorkbook()
# Write each data frame to a separate sheet
for (i in seq_along(func_annot)) {
  addWorksheet(wb, sheetName = names(func_annot)[i]) # Use the list names as sheet names
  writeData(wb, sheet = names(func_annot)[i], x = func_annot[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, file = "fgsea_custom_gsets_ora.xlsx", overwrite = TRUE)



# Check other customized functional gene sets from different papers
gmt_df <- read.gmt("~/Documents/VAR_score/functional_signatures.gmt") 
gmt_list <- split(x = gmt_df$gene, f = gmt_df$term)
func_annot <- lapply(mp_list, function(program) {
  fgRes <- fgsea::fora(pathways = gmt_list,
                       genes = program,
                       universe = universe)
  
  fgRes <- fgRes[fgRes$pval <= 0.05,]
  return(fgRes)
})

# Save results as a excel file
wb <- createWorkbook()
# Write each data frame to a separate sheet
for (i in seq_along(func_annot)) {
  addWorksheet(wb, sheetName = names(func_annot)[i]) # Use the list names as sheet names
  writeData(wb, sheet = names(func_annot)[i], x = func_annot[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, file = "fgsea_var_gsets_ora.xlsx", overwrite = TRUE)


# Check overlap between metaprograms
pairwise_intersections <- sapply(names(mp_list), function(x) {
  sapply(names(mp_list), function(y) {
    length(intersect(mp_list[[x]],mp_list[[y]]))
  })
})

# Convert to a matrix for readability
pairwise_matrix <- matrix(pairwise_intersections, nrow = length(mp_list),
                          dimnames = list(names(mp_list), names(mp_list)))
print(pairwise_matrix)
