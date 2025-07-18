# BACKEND: It's a good practice to have every component related with BACKEND separated from the FRONTEND scripts.

# Installing packages:
# install.packages("dplyr")
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("DT")
# install.packages("Seurat")
# install.packages("SeuratObject")
# install.packages("SeuratDisk")
# install.packages("ggplot2")
# install.packages("hdf5r")
# install.packages("callr")

# Install and load packages automatically not included in the .yaml environment:
source("setup_shiny_app.R")

# Upload libraries:
library(remotes)
library(dplyr)
library(shiny)
library(shinydashboard)
library(DT)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(hdf5r)
library(ShinyCell)
library(callr)


# File size uploading threshold configuration (We'll use 60GB, our .h5ad file size is ~50GB):
options(shiny.maxRequestSize = 60 * 1024^2 * 1024) # ~60 GB


#Upload always available files for download
clinical_md      <- read.delim("www/clinical_metadata_v4_clean.tsv")
annotations_raw  <- read.delim("www/tcca_annotation_raw.tsv")
subclone_tsv     <- read.delim("www/subclone_level_annotated.tsv")

# .h5ad file for representation path:
h5ad_raw_path    <- "www/tcca_raw_mat.h5ad"
path <- "/home/lmgonzalezb/Documents/bc-meta/bc-meta_repo/bc-meta/src"
shinycell_subapp_path <- paste0(path, "Shiny_TCCA_App/shinyAppH5ad")
main_app_path <- paste0(path, "Shiny_TCCA_App")



