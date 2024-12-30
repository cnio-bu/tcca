library(remotes)
#remotes::install_github("carmonalab/GeneNMF") #from Github
library(GeneNMF)
library(Seurat)
library(BPCells)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)

setwd("/storage/scratch01/shared/projects/bc-meta/")

# Load Seurat object
seu.lvl2 <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")

# Select malignant cells
malignant <- subset(seu.lvl2, subset = malignancy == TRUE)
colnames(malignant) <- paste0("c", 1:ncol(malignant))

# Subset malignant cells with computed Beyondcell Scores
# bc <- readRDS("beyondcell/results/beyondcell_pancancer_final.Rds")
# malignant <- subset(seu.lvl2, cells = colnames(bc))

# Normalize data from each sample to run NMF later
malignant <- Seurat::NormalizeData(malignant,
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000)
malignant <- Seurat::FindVariableFeatures(malignant, selection.method = "vst", nfeatures = 7000)              
malignant$study_sample <- paste0(malignant$study, "__", malignant$sample) 

# Create a list of Seurat objects (one for each sample)
malignant <- JoinLayers(malignant)
print("Layers joined")
malignant <- SplitObject(malignant, split.by = "study_sample")
print("Object splitted into samples")

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf")
saveRDS(malignant, "malignant_split_allsamples.rds")

# Compute NMF for each sample
geneNMF.programs <- GeneNMF::multiNMF(malignant, assay = "RNA", slot = "data", k = 4:9, nfeatures = 7000)
print("NMF computed for each sample")

saveRDS(geneNMF.programs, "geneNMFprograms_allsamples.rds")


# geneNMF.programs <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf/geneNMFprograms.rds")

# geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                         nMP = 7,
#                                         weight.explained = 0.7,
#                                         max.genes = 100)

# saveRDS(geneNMF.metaprograms, "/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf/geneNMFmetaprograms.rds")

# geneNMF.metaprograms <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf/geneNMFmetaprograms_jaccard.rds")

# setwd("/storage/scratch01/users/mgonzalezb/bc-meta/functional_nmf/")
# ph <- plotMetaPrograms(geneNMF.metaprograms)
# ggsave("plots/plot_mps_10_jccard.png", plot = ph, dpi = 300, height = 10, width = 10)