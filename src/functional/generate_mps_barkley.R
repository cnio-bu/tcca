library(igraph)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(scales)
# Load NMF matrices
setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf")

nmf_programs <- readRDS("geneNMFprograms_allsamples_4-30.rds")

# Select programs from k = 4..10,15, 20, 25, 30
pattern <- paste(paste0("k", c(4:10, 15, 20, 25, 30)),collapse = "|")
nmf_programs <- nmf_programs[grep(pattern, names(nmf_programs), value = TRUE)]

#-------------------------------------------------------------------------------
# Define functions
# ------------------------------------------------------------------------------
NMFToModules <- function(res, gmin = 5) {
    # Extract scores and coefficients
    scores <- res$w
    coefs <- res$h

    # Remove if fewer than gmin genes
    ranks_x <- t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
    ranks_y <- apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)

    for (i in 1:ncol(scores)) {
        ranks_y[ranks_x[, i] > 1, i] <- Inf
    }

    modules <- apply(ranks_y, 2, function(m) {
        a <- sort(m[is.finite(m)])
        a <- a[a == 1:length(a)]
        names(a)
    })

    l <- sapply(modules, length)
    keep <- (l >= 5)

    if (!any(keep)) {
        warning(
            "No modules passed the gene filter (gmin = ",
            gmin, "). Returning an empty list."
        )
        return(list())
    }

    scores <- scores[, keep, drop = FALSE]
    coefs <- coefs[keep, , drop = FALSE]

    # Find modules again with filtered data
    ranks_x <- t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
    ranks_y <- apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)

    cat("Dimensions of ranks_x:", dim(ranks_x), "\n")
    
    for (i in 1:ncol(scores)) {
        ranks_y[ranks_x[, i] > 1, i] <- Inf
    }

    modules <- apply(ranks_y, 2, function(m) {
        a <- sort(m[is.finite(m)])
        a <- a[a == 1:length(a)]
        names(a)
    })

    # Handle cases where only one module passes the filter
    if (is.matrix(modules) && ncol(modules) == 1) {
        modules <- list(as.vector(modules))
    }
    
    # Name and format the modules
    # names(modules) <- sapply(modules, "[", 1)
    # names(modules) <- paste("m", names(modules), sep = "_")
    # names(modules) <- gsub("-", "_", names(modules))

    return(modules)
}


#-------------------------------------------------------------------------------
# Adapt RcppML::nmf output to define MPs
# ------------------------------------------------------------------------------
sample_names <- gsub("\\.k\\d+$", "", names(nmf_programs))
sample_split <- split(seq_along(nmf_programs), sample_names)

# Get a list where each entry is a matrix with 2000 rows (genes) x 39 columns
# (NMF programs)
sample_nmf <- lapply(names(sample_split), function(sample) {
    sample_programs <- nmf_programs[sample_split[[sample]]]
    sample_programs <- lapply(names(sample_programs), function(k) {
        k_programs <- sample_programs[[k]]
        colnames(k_programs$w) <- paste(k, colnames(k_programs$w), sep = "_")
        rownames(k_programs$h) <- paste(k, rownames(k_programs$h), sep = "_")
        return(k_programs)
    })
    names(sample_programs) <- paste0(c(4:10, 15, 20, 25, 30))
    return(sample_programs)
})

names(sample_nmf) <- names(sample_split)


# Select rank
select.rank <- lapply(sample_nmf, function(sample.ranks) {
    modules.list <- lapply(sample.ranks, NMFToModules, gmin = 5)
    print(sapply(modules.list, length))
    comp <- as.numeric(names(modules.list)) - sapply(modules.list, length)
    mi <- min(comp)
    r <- names(which(comp == mi))
    r <- r[length(r)]
    print(r)
    rank <- sample.ranks[[r]]
    return(rank)
})

modules.list <- lapply(select.rank, NMFToModules)

## Modules from graph ##
seu.lvl2 <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
genes.all <- rownames(seu.lvl2)
all <- unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) <- unlist(sapply(modules.list, names))
ta <- table(unlist(all))
genes.use <- names(ta)[ta > 1]

# Filter non-overlapping modules
for (i in 1:5) {
    all <- unlist(modules.list, recursive = FALSE, use.names = TRUE)
    all <- lapply(all, intersect, genes.all)
    sim <- sapply(all, function(x) {
        sapply(all, function(y) {
            length(intersect(x, y)) / length(union(x, y))
        })
    })
    keep <- rownames(sim)[apply(sim, 1, function(x) {
        sum(x > 0.05) >= 3
    })]
    all <- all[keep]
    modules.list <- lapply(names(modules.list), function(x) {
        li <- modules.list[[x]]
        li[names(li)[paste(x, names(li), sep = ".") %in% keep]]
    })
    names(modules.list) <- names(select.rank)
    ta <- table(unlist(all))
    genes.use <- names(ta)[ta > 1]
    print(length(all))
}

saveRDS(modules.list, "barkley_ksubset/modules_list_barkley.rds")
saveRDS(sim, "barkley_ksubset/similarity_matrix_barkley.rds")
saveRDS(all, "barkley_ksubset/all_barkley.rds")

# Adjacency matrix, list by cancer
cancer_types <- seu.lvl2@meta.data %>%
    mutate(study_sample = paste0(study, "__", sample), tumor_type = factor(tumor_type)) %>%
    select(study_sample, tumor_type) %>%
    unique() %>%
    mutate(order = match(study_sample, names(modules.list))) %>%
    arrange(order) %>%
    select(-order) %>%
    pull(tumor_type)

adj <- matrix(0, nrow = length(genes.use), ncol = length(genes.use))
adj.list <- list()
for (can in levels(cancer_types)) {
    sub <- matrix(0, nrow = length(genes.use), ncol = length(genes.use))
    rownames(sub) <- genes.use
    colnames(sub) <- genes.use
    for (s in names(modules.list)[cancer_types == can]) {
        for (mod in modules.list[[s]]) {
            mod <- intersect(mod, genes.use)
            for (x in mod) {
                for (y in mod) {
                    sub[x, y] <- sub[x, y] + 1
                }
            }
        }
    }
    diag(sub) <- 0
    adj.list[[can]] <- sub
    # adj = adj + (sub > 0)
    adj <- adj + sub
}
adj_keep <- adj

# Remove low connections
adj[] <- (adj >= 2)
# adj[adj <= 1] = 0
for (i in 1:5) {
    keep <- names(which(rowSums(adj) >= 2))
    adj <- adj[keep, keep]
    print(dim(adj))
}


# Cluster
g <- graph_from_adjacency_matrix(adj, diag = FALSE, mode = "undirected", weighted = TRUE)
modules <- communities(cluster_infomap(g, nb.trials = 200))
names(modules) <- paste0("m_", sapply(modules, "[", 1))

saveRDS(modules, "barkley_ksubset/modules_barkley.rds")







modules <- modules[unlist(lapply(modules, function(x) length(x) >= 9))]

## Plot oiverlap between modules
ovlp_consensus <- sapply(all, function(x){
    sapply(modules, function(y){
        pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})

ovlp_consensus[is.infinite(ovlp_consensus)] <- max(ovlp_consensus[is.finite(ovlp_consensus)])
ovlp_consensus <- ovlp_consensus[,apply(ovlp_consensus, 2, function(x){any(x > 3)})]
df <- data.frame('sample' = sapply(colnames(ovlp_consensus), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))

df$top_consensus <- apply(ovlp_consensus, 2, which.max)
df$top <- factor(names(modules)[df$top], levels = names(modules))

top_ann <- HeatmapAnnotation(sample = df[, c('sample')], which = 'column', show_legend = FALSE)
side_ann <- HeatmapAnnotation(onsensus_module = df[, c('top')], which = 'row')

ovlp <- sapply(all, function(x){
  sapply(all, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
ovlp[is.infinite(ovlp)] <- max(ovlp[is.finite(ovlp)])
ovlp <- ovlp[colnames(ovlp_consensus),colnames(ovlp_consensus)]

colors_heatmap <- colorRamp2(seq(0, 6, length = 7), brewer_pal(palette = "RdBu", direction = -1)(7))
h <- Heatmap(name = 'Module overlap p-value', sim, 
            top_annotation = top_ann,
            right_annotation = side_ann,
            row_order = order(df$top),
            column_order = order(df$top),
            show_column_names = FALSE,
            show_row_names = FALSE)
            #show_row_names = TRUE, row_names_gp = gpar(cex = 0.4)

png(
  file = "barkley_ksubset/heatmap_mps.png",
  res = 500,
  width = 16,
  height = 14,
  units = "in"
)
draw(h, annotation_legend_side = "right")
dev.off() 