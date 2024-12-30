# # ----------------------------------------------------------------------------------------------------
# # Load packages and seurat objects
# # ----------------------------------------------------------------------------------------------------
library(Seurat)
library(BPCells)
library(reshape2)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ComplexHeatmap)
setwd("/home/lmgonzalezb/Documents/bc-meta/")
# source("fork/bc-meta/src/figures/TCCA_palette.R")
source(
  "fork/bc-meta/src/pancancer_analyses/functional_pancancer/robust_nmf_programs.R"
)

##---------------------------------FUNCTIONS ---------------------------------##
read.gmt <- function(gmt_file) {
  sigs_list <- list()
  sigs <- scan(gmt_file, what = character(), sep = "\n")
  for (sig in sigs) {
    sig <- unlist(strsplit(sig, "\t"))
    sig <- unique(sig[nzchar(sig)])
    sigs_list[[sig[1]]] <- sig[3:length(sig)]
  }
  return(sigs_list)
}

## Input:
# geneNMF.programs is a list in which each entry contains NMF gene-scores of a single sample and k value. In our study  # we ran NMF using ranks 4-9 on the top 7000 genes in each sample. Hence each entry in geneNMF.programs is a matrix
# with 7000 rows (genes) X k columns (NMF programs). We transform this list of matrices into a list where each element
# corrresponds to all the NMF gene-scores of a single sample (7000 genes x 39 programs).
# We define MPs in 2 steps:
# 1) The function robust_nmf_programs.R performs filtering, so that programs selected for defining MPs are:
#    i) Robust                - recur in more that one rank within the sample
#    ii) Non-redundant        - once a NMF program is selected, other programs within the sample that are similar to it are removed
#    iii) Not sample specific - has similarity to NMF programs in other samples
# ** Please see https://github.com/gabrielakinker/CCLE_heterogeneity for more details on how to define robust NMF programs
# 2) Selected NMFs are then clustered iteratively. At the end of the process, each cluster generates a list of the 50
# genes (i.e. the MP) that represent the NMF programs that contributed to the cluster. Notably, not all initially selected NMFs end up participating in a cluster.


# Load NMF matrices
geneNMF.programs <- readRDS("./functional_mps/geneNMFprograms_allsamples.rds")

#-----------------------------------------------------------------------------------------------------
# Adapt GeneNMF::multiNMF output to define MPs
# ----------------------------------------------------------------------------------------------------
geneNMF.programs <- lapply(geneNMF.programs, function(sample_k)
  sample_k$w)
sample.names <- gsub("\\.k\\d+$", "", names(geneNMF.programs))
sample.split <- split(seq_along(geneNMF.programs), sample.names)

# Get a list where each entry is a matrix with 2000 rows (genes) x 39 columns (NMF programs)
sample.NMFprograms <- lapply(names(sample.split), function(sample) {
  sample.programs <- geneNMF.programs[sample.split[[sample]]]
  names(sample.programs) <- paste0(sample, "_", "k", c(4:9))
  sample.programs <- lapply(names(sample.programs), function(n) {
    colnames(sample.programs[[n]]) <- paste(n, seq(1, ncol(sample.programs[[n]])), sep = ".")
    return(sample.programs[[n]])
  })
  sample.programs <- Reduce(f = cbind, x = sample.programs)
  return(sample.programs)
})

names(sample.NMFprograms) <- names(sample.split)


#-----------------------------------------------------------------------------------------------------
# Select NMF programs
# ----------------------------------------------------------------------------------------------------

## Parameters
intra_min_parameter <- 35
intra_max_parameter <- 10
inter_min_parameter <- 10


# get top 50 genes for each NMF program
nmf_programs <- lapply(sample.NMFprograms, function(x)
  apply(x, 2, function(y)
    names(sort(y, decreasing = T))[1:50]))

# for each sample, select robust NMF programs (i.e. observed using different
# ranks in the same sample), remove redundancy due to multiple ranks, and apply
# a filter based on the similarity to programs from other samples.
nmf_filter_ccle <- robust_nmf_programs(
  nmf_programs,
  intra_min = intra_min_parameter,
  intra_max = intra_max_parameter,
  inter_filter = T,
  inter_min = inter_min_parameter
)
print(paste("Robust_programs_selected!", format(Sys.time(), "%a %b %d %X %Y")))

nmf_programs <- lapply(nmf_programs, function(x)
  x[, is.element(colnames(x), nmf_filter_ccle), drop = F])
nmf_programs <- do.call(cbind, nmf_programs)


# calculate similarity between programs
nmf_intersect <- apply(nmf_programs , 2, function(x)
  apply(nmf_programs , 2, function(y)
    length(intersect(x, y))))

# hierarchical clustering of the similarity matrix
nmf_intersect_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), 
                            colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), 
                               order.dendrogram(nmf_intersect_hc)]

print(paste(
  "Similarity between programs computed by intersections",
  format(Sys.time(), "%a %b %d %X %Y")
))

saveRDS(nmf_programs, "./functional_mps/nmf_programs_mat_allsamples.rds")
saveRDS(nmf_intersect_hc, "./functional_mps/nmf_intersect_hc_allsamples.rds")
saveRDS(nmf_intersect, "./functional_mps/nmf_intersect_allsamples.rds")


nmf_programs <- readRDS("./functional_mps/nmf_programs_mat_allsamples.rds")
nmf_intersect_hc <- readRDS("./functional_mps/nmf_intersect_hc_allsamples.rds")
nmf_intersect <- readRDS("./functional_mps/nmf_intersect_allsamples.rds")


# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------

### Parameters for clustering
min_intersect_initial <- 15   # the minimal intersection cutoff for defining the first NMF program in a cluster
min_intersect_cluster <- 15   # the minimal intersection cutoff for adding a new NMF to the forming cluster
min_group_size <- 10     # the minimal group size to consider for defining the first NMF program in a cluster

sorted_intersection <- sort(apply(nmf_intersect , 2, function(x)
  (length( which(x >= min_intersect_initial)) - 1)) , decreasing = TRUE)

cluster_list <- list()   # every entry contains the NMFs of a chosen cluster
mp_list <- list()
k <- 1
curr_cluster <- c()
nmf_intersect_original <- nmf_intersect

while (sorted_intersection[1] > min_group_size) {
  curr_cluster <- c(curr_cluster , names(sorted_intersection[1]))
  
  # intersection between all remaining NMFs and genes in MP
  # Genes in the forming MP are first chosen to be those in the first NMF. 
  # Genes_MP always has only 50 genes and evolves during the formation of the cluster
  genes_mp <- nmf_programs[, names(sorted_intersection[1])]
  # remove selected NMF
  nmf_programs <- nmf_programs[, -match(names(sorted_intersection[1]) , colnames(nmf_programs))]
  # intersection between all other NMFs and genes_mp
  intersection_with_genes_mp  <- sort(apply(nmf_programs, 2, function(x)
    length(intersect(genes_mp, x))) , decreasing = TRUE)
  # has genes in all NMFs in the current cluster, for redefining genes_mp after adding a new NMF
  nmf_history <- genes_mp
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). 
  # When the number of genes with a given frequency span beyond the 50th genes, they are sorted according to their NMF score.
  while (intersection_with_genes_mp[1] >= min_intersect_cluster) {
    curr_cluster <- c(curr_cluster , names(intersection_with_genes_mp)[1])
    
    ## Genes_MP is newly defined each time according to all NMFs in the current cluster
    genes_mp_temp <- sort(table(c(nmf_history , 
                                  nmf_programs[, names(intersection_with_genes_mp)[1]])), 
                          decreasing = TRUE)
    
    ### genes with overlap equal to the 50th gene
    genes_at_border <- genes_mp_temp[which(genes_mp_temp == genes_mp_temp[50])]
    
    if (length(genes_at_border) > 1) {
      ### Sort last genes in genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in curr_cluster and extract NMF scores for each gene
      genes_curr_nmf_score <- c()

      for (i in curr_cluster) {
        curr_sample <- gsub("_k\\d\\.\\d", "", i)
        matched_indices <- match(names(genes_at_border),
                                 rownames(sample.NMFprograms[[curr_sample]]))
        matched_genes <- matched_indices[!is.na(matched_indices)]
        q <- sample.NMFprograms[[curr_sample]][matched_genes, i]
        ### sometimes when adding genes the names do not appear
        names(q) <- rownames(sample.NMFprograms[[curr_sample]])[matched_genes]
        genes_curr_nmf_score <- c(genes_curr_nmf_score, q)
      }

      genes_curr_nmf_score_sort <- sort(genes_curr_nmf_score, decreasing = TRUE)
      genes_curr_nmf_score_sort <- genes_curr_nmf_score_sort[unique(names(genes_curr_nmf_score_sort))]

      genes_mp_temp <- c(names(genes_mp_temp[which(genes_mp_temp > genes_mp_temp[50])]) ,
                         names(genes_curr_nmf_score_sort))
      
    } else {
      genes_mp_temp <- names(genes_mp_temp)[1:50]
    }
    
    nmf_history  <- c(nmf_history , nmf_programs[, names(intersection_with_genes_mp)[1]])
    genes_mp  <- genes_mp_temp[1:50]

    # remove selected NMF
    nmf_programs <- nmf_programs[, -match(names(intersection_with_genes_mp)[1] ,
                                          colnames(nmf_programs))]
    
    # intersection between all other NMFs and genes_mp
    intersection_with_genes_mp <- sort(apply(nmf_programs, 2, function(x)
      length(intersect(genes_mp, x))) , decreasing = TRUE) 

  }
  
  cluster_list[[paste0("MP", k)]] <- curr_cluster
  mp_list[[paste0("MP", k)]] <- genes_mp
  k <- k + 1

  # Remove current chosen cluster
  nmf_intersect  <- nmf_intersect[-match(curr_cluster, rownames(nmf_intersect)) , 
                                  -match(curr_cluster, colnames(nmf_intersect))]
  
  # Sort intersection of remaining NMFs not included in any of the previous clusters
  sorted_intersection <-  sort(apply(nmf_intersect , 2, function(x)
    (length(
      which(x >= min_intersect_initial)
    ) - 1)) , decreasing = TRUE)
  
  curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}
print(paste("MPs computed!", format(Sys.time(), "%a %b %d %X %Y")))
saveRDS(cluster_list, "functional_mps/cluster_list_allsamples.rds")
saveRDS(mp_list, "functional_mps/mp_list_allsamples.rds")

cluster_list <- readRDS("./functional_mps/cluster_list_allsamples.rds")
mp_list <- readRDS("./functional_mps/mp_list_allsamples.rds")

# ----------------------------------------------------------------------------------------------------
# Heatmap of Jaccard similaity scores for NMF programs clustered into MPs
# ----------------------------------------------------------------------------------------------------
# Reorder MPs to group related ones together and remove MP11, linked only to brain metastases of 
# melanoma (with skin pigmentation enrichment)
mps <- c("MP1", "MP8", "MP14", "MP3", "MP13", "MP2", "MP9", "MP4", "MP10", "MP12", 
         "MP5", "MP6",  "MP7")
cluster_list <- cluster_list[mps]
names(cluster_list) <- paste0("MP", 1:13)
# Transform cluster list into a vector
cl_members <- unlist(lapply(seq_along(cluster_list), function(i) {
  group <- names(cluster_list)[i]  # Get the name of the current group
  setNames(rep(group, length(cluster_list[[i]])), cluster_list[[i]])
}))

# Add tumor sample annotation to the heatmap
cluster_annotation <- data.frame(
  cluster = cl_members,
  stringsAsFactors = FALSE
)
programs_names <- rownames(cluster_annotation)
cluster_annotation$study_sample <- gsub("_k\\d+\\.\\d+", "", rownames(cluster_annotation))
clinical <- data.table::fread(
    "clinical_metadata_v4_clean.tsv"
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

cluster_annotation <- cluster_annotation %>% 
left_join(clinical, by = "study_sample") %>% 
mutate(
        summarised_tumor_site = case_when(
            refined_tumor_site %in% names(translat_human_sites) ~ refined_tumor_site,
            TRUE ~ "other"
        ),
        adult_pediatric = ifelse(age >= 16, "Adult", "Pediatric"),
        is_blood = ifelse(tumor_type %in% c("ALL", "CLL", "LAML","MM"), "Liquid", "Solid"),
        treated = ifelse(treated == "t", "Treated", "Untreated"),
        sex = ifelse(sex == "f", "Female", "Male"),
        sample_type = ifelse(sample_type == "m", "Metastasis", "Primary"),
        cluster = factor(cl_members, levels = names(cluster_list) )
)

cluster_annotation$summarised_tumor_site <-  translat_human_sites[cluster_annotation$summarised_tumor_site]

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

rownames(cluster_annotation) <- programs_names

colnames(cluster_annotation) <- c(
    "Sex",
    "Age group",
    "Solid/Liquid",
    "Sample type",
    "Sample site",
    "Treatment",
    "Meta-program"
)

pals <- list(
    "Sex" = sex_colors,
    "Age group" = age_colors,
    "Solid/Liquid" = sl_colors,
    "Sample type" = pm_colors,
    "Sample site" = tumor_sites_colors,
    "Treatment" = treatment_colors,
    "Meta-program" = mps_colors
)


top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df =  cluster_annotation,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 16),
    annotation_legend_param = list( "Meta-program" = list(title_gp = gpar(fontsize = 18, fontface = "bold"),  # Title font size
                                                          labels_gp = gpar(fontsize = 18))),
    show_legend = c("Sex" = FALSE,    
                    "Age group" = FALSE,
                    "Solid/Liquid" = FALSE,
                    "Sample type" = FALSE,
                    "Sample site" = FALSE,
                    "Treatment" = FALSE,
                    "Meta-program" = TRUE))

#  Sort Jaccard similarity plot according to new clusters:
inds_sorted <- c()
for (j in 1:length(cluster_list)) {
  inds_sorted <- c(inds_sorted , match(cluster_list[[j]], colnames(nmf_intersect_original)))
  
}

heat <- ComplexHeatmap::Heatmap(
    mat = nmf_intersect_original[inds_sorted, inds_sorted],
    col = colorRamp2(c(2, 25), hcl_palette = "Inferno", reverse = TRUE),
    # right_annotation = right_annotation,
    top_annotation = top_annotation,
    row_order = rownames(cluster_annotation[order(cluster_annotation$`Meta-program`), , drop = FALSE]),
    cluster_rows = FALSE,
    cluster_row_slices = TRUE,
    row_split = cluster_annotation$`Meta-program`,
    row_title = NULL,
    column_order = rownames(cluster_annotation[order(cluster_annotation$`Meta-program`), , drop = FALSE]),
    cluster_columns = FALSE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split =  cluster_annotation$`Meta-program`, 
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_side = "top",
    column_title = NULL,
    heatmap_legend_param = list(title = "Similarity\n(Jaccard index)", 
                                title_gp = gpar(fontsize = 18, fontface = "bold"),    # Title font size for the heatmap legend
                                labels_gp = gpar(fontsize = 18)),
    heatmap_width = unit(10, "in"),
    heatmap_height = unit(10, "in")
)
png(
    file = "figures/functional/all_samples_nmf/heatmap_mps_15.png",
    res = 500,
    width = 16, 
    height = 14,
    units = "in"
)
draw(heat, annotation_legend_side = "right")
dev.off()



# Annotate MPs
library(msigdbr)
library(fgsea)
library(openxlsx)

seu.lvl2 <- readRDS("seu_lvl2_sex_inferred.rds")
msig_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
msig_list <- split(x = msig_df$gene_symbol, f = msig_df$gs_name)

func_annot <- lapply(mp_list, function(program) {
  names <- program
  program <- seq(1, length(program))
  names(program) <- program
  fgRes <- fgsea::fora(pathways = msig_list,
                       genes = program,
                       universe = rownames(seu.lvl2),
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
saveWorkbook(wb, file = "functional_mps/fgsea_hallmarks_ora.xlsx", overwrite = TRUE)

# Check MPs from Kinker et al, Gavish et al., and Barkey et al and other functional gene sets
gmt_list <- read.gmt("fork/bc-meta/reference/combined_gsets_functional.gmt") 
func_annot <- lapply(mp_list, function(program) {
  fgRes <- fgsea::fora(pathways = gmt_list,
                       genes = program,
                       universe = rownames(seu.lvl2))
  
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
saveWorkbook(wb, file = "functional_mps/fgsea_custom_gsets_ora.xlsx", overwrite = TRUE)



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