library("beyondcell")
library("fabia")
library("tidyverse")

### SNAKEMAKE I/O ###
bc_list <- snakemake@input[["bc_list"]]
moas_table <- snakemake@input[["moas_table"]]
where_to_save <- snakemake@output[["module_dir"]]

all_samples <- readRDS(bc_list)
moas <- read.table(moas_table, header = T, sep ="\t")
moas <- filter(moas,  studies != "CTRP")
sigs_to_keep <- moas$IDs

for(sample in all_samples){
    
    ## get norm mat.
    normalized_bc <- sample@normalized
    normalized_bc <- subset(normalized_bc, rownames(normalized_bc) %in% sigs_to_keep)

    ## C optimized code, let it run in single thread
    res <- fabia::fabias(
        X = normalized_bc,
        p = 50, ## Hidden factors = biclusters. Max. to 50 (581 / 50 ~ 11 drugs)
        cyc = 500, ## iterations, keep it at 500
        spz = 0.5, ## minimum sparseness, Laplace.
        non_negative = 0, ## Allow negative factors a.k.a. negative vectors
        random = 1.0, ## allow random initialization of loadings,
        center = 2, ## median centering
        lap = 1, ## minimal value of the variational param
        nL = 1, ## do not allow drugs to be in > 1 biclust, it will report mirrowed biclusts otherwise,
        lL = 100, ## do not allow biclusters > 100 drugs. It tends to aggregate spurious clusters
    )
    
    ## save the objects before table extraction, just in case
    dir.create(where_to_save, showWarnings = FALSE)

    sample_name <- unique(sample@meta.data$sample)
    save_dir <- paste0(where_to_save, "/", sample_name, "_biclusters.rds")
    saveRDS(object = res, file = save_dir)
    
    ## expand the module table and save in tab format.
    biclusters <- extractBic(fact = res, thresZ = 0.5)
    
    all_biclusters <- list()
    all_cells <- list()
    for(i in c(1:50)){
        this_biclust <- biclusters$bic[i, ]
        ## check length of the bicluster rowwise (aka drugs)
        if(length(this_biclust$bixn) <= 5) {
            next
        }else {
            named_list_sigs <- this_biclust$bixn
            sig_contributions <- this_biclust$bixv
            named_list_cells <- this_biclust$biypn
            cell_contributions <- this_biclust$biypv
            names(sig_contributions) <- named_list_sigs
            names(cell_contributions) <- named_list_cells
            all_biclusters[[i]] <- sig_contributions
            all_cells[[i]] <- cell_contributions
        }
    }

    bicluster_table <- enframe(all_biclusters) %>%
        unnest_longer(col = "value") %>%
        rename(
            "bicluster" = name,
            "cluster_contribution" = value,
            "signature" = value_id
        )

    cell_table <- enframe(all_cells) %>%
        unnest_longer(col = "value", simplify = TRUE) %>%
        rename(
        "bicluster" = name,
        "cluster_cell_contribution" = value,
        "cell_name" = value_id
    )
    
    cell_table <- cell_table %>%
        right_join(
            bicluster_table,
            by = "bicluster",
            relationship = "many-to-many"
            )
    
    cell_table$information_content <- res@avini[cell_table$bicluster]
    
    save_table <- paste0(where_to_save, "/", sample_name, "_clusters.tsv")
    ## save file
    write.table(
        x = cell_table,
        file = save_table,
        sep="\t",
        row.names = FALSE
        )
}
