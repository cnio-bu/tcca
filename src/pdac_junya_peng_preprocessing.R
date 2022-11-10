library("tidyverse")
library("Seurat")

## load common funs

source("src/sc_functions.R")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

#said the working directory  
data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/metanalisis_bc")

#load data
load(paste0(data_folder,
    "/single_cell/raw/pdac_junya_peng/pdac_junya_peng.RData")
)

#check version of seurat object
Version(sc_Peng_PDAC)
#version is 4.0.1


#Add anotation to cells
anot <- read.table(paste0(data_folder,
    "/single_cell/raw/pdac_junya_peng/all_celltype.txt"),
    sep = "\t",
    header = T,
    row.names = 1)


sc <- sc_Peng_PDAC


sc_Peng_PDAC@meta.data <- merge(sc_Peng_PDAC@meta.data,
    anot,
    by = "row.names",
    all.x = T)



row.names(sc_Peng_PDAC@meta.data) <- sc_Peng_PDAC@meta.data$Row.names
sc_Peng_PDAC@meta.data <- sc_Peng_PDAC@meta.data[,-1]

# Change default ident for later split and QC
sc_Peng_PDAC@meta.data$orig.ident <- sc_Peng_PDAC@meta.data$Sample
Idents(sc_Peng_PDAC) <- "orig.ident"


sc_Peng_PDAC@meta.data <- sc_Peng_PDAC@meta.data[rownames(sc@meta.data),]

## Split the merged obj

sc_Peng_PDAC.samples <- Seurat::SplitObject(
    object = sc_Peng_PDAC,
    split.by = "orig.ident"
    )


unique(sc_Peng_PDAC.samples$T2@meta.data$Sample)



#Filtering and normalize data
sc_Peng_PDAC_filtered <- lapply(X = sc_Peng_PDAC.samples,
                            FUN = filter_sc,    
                            res_dir = paste0(data_folder,
                            "/single_cell/qc/pdac_junya_peng/")
                            )

sc_Peng_PDAC_filtered <- lapply(X = sc_Peng_PDAC_filtered,
                            FUN = normalize_and_scale
)

# Get rid of "NULL" samples, there are no cells left in these
sc_Peng_PDAC_filtered[sapply(sc_Peng_PDAC_filtered, is.null)] <- NULL


## Filter and subset tumor samples

sc_Peng_PDAC_filtered <- sc_Peng_PDAC_filtered[names(sc_Peng_PDAC_filtered) 
%in% paste0('T',1:24)]


keep_malignant <- function(sc) {
    sc_filtered <- subset(
        x = sc,
        subset = cluster == "Ductal cell type 2"
        )
    return(sc_filtered)
}

sc_Peng_PDAC_malignant <- lapply(sc_Peng_PDAC_filtered, keep_malignant)    

# Get rid of "NULL" samples, if there were no malignants left
sc_Peng_PDAC_malignant[sapply(sc_Peng_PDAC_malignant, is.null)] <- NULL

saveRDS(
    object = sc_Peng_PDAC_malignant, 
    file = paste0(data_folder, "/single_cell/obj/pdac_junya_peng/all_malignant.rds")
)
