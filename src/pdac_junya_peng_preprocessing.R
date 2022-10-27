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
data_folder <- paste0(data_folder, "/single_cell")

#load data
load(paste0(data_folder,
                               "/raw/pdac_junya_peng/pdac_junya_peng.RData")
)






data_folder
"/home/ccarreterop/Documents/metanalisis_bc_todos/metanalisis_bc/single_cell/raw/"