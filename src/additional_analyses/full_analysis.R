library("beyondcell")
library("Seurat")
library("tidyverse")

gs <- GetCollection(SSc, n.genes = 250)

merged_sunny_dalia <- readRDS("raw/merged_dalia_sunny_hugo.rds")

only_dalia <- unique(merged_sunny_dalia$orig.ident)
only_dalia <- only_dalia[1:62]

dalia_merged <- subset(x = merged_sunny_dalia, orig.ident %in% only_dalia)
dalia_malignant <- subset(x = dalia_merged, type == "malignant")

rm(merged_sunny_dalia)
gc()

dalia_malignant <- FindVariableFeatures(object = dalia_malignant)
dalia_malignant <- ScaleData(object = dalia_malignant, features = rownames(dalia_malignant))
dalia_malignant <- RunPCA(object = dalia_malignant)
dalia_malignant <- FindNeighbors(object = dalia_malignant)
dalia_malignant <- FindClusters(object = dalia_malignant, resolution = 0.5)
dalia_malignant <- RunUMAP(dalia_malignant, dims = 1:10)


seu <- DimPlot(object = dalia_malignant, reduction = "umap", group.by = "sample")

erlotinib <- gs@genelist$"sig-21305"$up

dalia_malignant <- AddModuleScore(
    object = dalia_malignant,
    features = list(erlotinib),
    ctrl = 25,
    assay = "RNA",
    name = "erlotinib_bc"
)

erlo_umap <- FeaturePlot(object = dalia_malignant, reduction = "umap", features = "erlotinib_bc1")

ggsave(plot = erlo_umap, filename = "full_seurat_umap_by_sample.png")

bc <- bcScore(
    sc = dalia_malignant,
    gs = GetCollection(SSc, include.pathways = FALSE),
    expr.thres = 0.1
)

saveRDS(object = bc, file = "dalia_merged_bc.rds")

bc@normalized[is.nan(bc@normalized)] <- 0
bc@data[is.nan(bc@data)] <- 0

## test
library("uwot")

test <- uwot::umap(
    X = t(bc@normalized),
    n_neighbors = 20,
    n_components = 2,
    metric = "correlation",
    scale = TRUE,
    spread = 2,
    min_dist = 0.2,
)

test <- as.data.frame(test)

bc@meta.data$my_umap_1 <- test$V1
bc@meta.data$my_umap_2 <- test$V2

a <- ggplot(data = bc@meta.data, aes(x = my_umap_1, y = my_umap_2)) + 
    geom_point(aes(color = cancer), size = .5, alpha = .5) + 
    theme_bw()



bc <- bcUMAP(bc = bc, pc = 5, k.neighbors = 10, npcs = 50, seed = 120394)

my_clusters <- bcClusters(
    bc = bc,
    idents = "bc_clusters_res.0.2",
    UMAP = "beyondcell"
)


my_features <- bcClusters(bc = bc, idents = "nFeature_RNA", factor.col = FALSE)
my_studies <- bcClusters(bc = bc, idents = "cancer")

bc <- bcRegressOut(bc = bc, vars.to.regress = "nFeature_RNA")

bc <- bcUMAP(bc = bc, pc = 5, k.neighbors = 20, res = 0.1)


my_clusters_rgrss <- bcClusters(
    bc = bc,
    idents = "bc_clusters_res.0.1",
    UMAP = "beyondcell"
)


test2 <- uwot::umap(
    X = t(bc@normalized),
    n_neighbors = 20,
    n_components = 2,
    metric = "correlation",
    scale = TRUE,
    spread = 2,
    min_dist = 0.2,
)


test2 <- as.data.frame(test2)

bc@meta.data$my_umap_1 <- test2$V1
bc@meta.data$my_umap_2 <- test2$V2

a <- ggplot(data = bc@meta.data, aes(x = my_umap_1, y = my_umap_2)) + 
    geom_point(aes(color = bc_clusters_res.0.1), size = .5, alpha = .5) + 
    theme_bw()


mat <- bc@normalized

dim(mat)

mat <- scale(x = mat, center = TRUE, scale = TRUE)
mat[is.na(mat)] <- 0
sv <- svd(x = t(mat))

U <- sv$u
V <- sv$v
D <- sv$d

## U are un-scaled PC, Z is scaled PC
Z <- t(mat) %*% V

## this is basically a glorified PCA

# X is now = U (pattern) * D(pattern) * V (amplitude matrix)
# Z = U*D

dim(V)

feature_loadings <- V
rownames(feature_loadings) <- rownames(mat)
feature_loadings <- feature_loadings[, 1:5]

library("tidyverse")

features <- feature_loadings %>%
    as_tibble(rownames = "drug") %>%
    pivot_longer(cols = V1:V5, names_to = "component", values_to = "loading") %>%
    mutate(abs_loading = abs(loading)) %>%
    arrange(desc(abs_loading))


write.csv(x = features, file = "features.csv")

feat_v1 <- features[features$component == "V1", "drug"]
feat_v1 <- head(feat_v1, n = 10) 
feat_v1 <- pull(feat_v1, "drug")

all_sigs <- bcSignatures(bc = bc, UMAP = "beyondcell", signatures = list(values = c(feat_v1)))

this_sig = 0
for (sig in all_sigs) {
    
    ggsave(plot = sig, filename = paste0(this_sig, ".png"))
    this_sig <- this_sig + 1
}


feat_v2 <- head(features[features$component == "V2", "drug"], n = 10) %>%
    pull("drug")


all_sigs2 <- bcSignatures(bc = bc, UMAP = "beyondcell", signatures = list(values = c(feat_v2)))

this_sig = 0
for (sig in all_sigs2) {
    
    ggsave(plot = sig, filename = paste0(this_sig, ".png"))
    this_sig <- this_sig + 1
    
}

bc_rgrss_cancer <- bcClusters(bc = bc, idents = "cancer", UMAP = "beyondcell")

ggsave(plot = bc_rgrss_cancer, filename = "bc_rgrss_by_cancer_type.png")

full_metadata <- bc@meta.data

write.csv(x = full_metadata, file = "full_metadata.csv")

bc <- bcRanks(bc = bc, idents = "bc_clusters_res.0.1")

squares <- bc4Squares(bc = bc, idents = "bc_clusters_res.0.1", levels(c("1", "2")))

ggsave(plot = squares, filename = "bc_squares.png")


# V1 top 50
top_50_f1 <- features %>%
    filter(component == "V1") %>%
    arrange(desc(abs_loading)) %>%
    head(50) %>%
    pull(drug)


top_for_contingency <- features %>%
    group_by(component) %>%
    arrange(desc(abs_loading)) %>%
    slice_head(n = 50) %>%
    ungroup()


## get moas
all_moas <- FindDrugs(bc = bc, x = top_for_contingency$drug)

top_for_contingency_annotated <- top_for_contingency %>%
    left_join(y = all_moas[, c("bc.names", "MoAs")], by = c("drug" = "bc.names")) %>%
    mutate(main_moa = stringr::str_remove_all(string = MoAs, pattern = ",.*"))


moa_arrest_v1 <- top_for_contingency_annotated %>%
    filter(component == "V1" & main_moa == "DNA replication inhibitor") %>%
    distinct()

moa_arrest_not_v1 <- top_for_contingency_annotated %>%
    filter(component != "V1" & main_moa == "DNA replication inhibitor") %>%
    distinct()


not_moa_arrest_v1 <- top_for_contingency_annotated %>%
    filter(component == "V1" & main_moa != "DNA replication inhibitor") %>%
    distinct()

not_moa_arrest_not_v1 <- top_for_contingency_annotated %>%
    filter(component != "V1" & main_moa != "DNA replication inhibitor") %>%
    distinct()


arrest_cont_table <- matrix(
    c(nrow(moa_arrest_v1),
      nrow(moa_arrest_not_v1),
      nrow(not_moa_arrest_v1),
      nrow(not_moa_arrest_not_v1)
      ),
    nrow = 2,
    dimnames = list(c("V1", "Others"),
                    c("Arrest", "Others")
                    )
)

my_test <- fisher.test(arrest_cont_table, alternative)


vd <- uwot::umap(X = b)


## new mat
new_mat <- bc@normalized
new_mat <- t(new_mat)
new_mat <- scale(new_mat, center = TRUE, scale = TRUE)

new_sv <- svd(x = new_mat)

new_loadings <- new_sv$v
rownames(new_loadings) <- rownames(bc@normalized)



features2 <- new_loadings %>%
    as_tibble(rownames = "drug") %>%
    pivot_longer(cols = V1:V555, names_to = "component", values_to = "loading") %>%
    mutate(abs_loading = abs(loading)) %>%
    arrange(desc(abs_loading))






### TUESDAY POST REU 
bc <- readRDS(file = "/local/sagarcia/bc-meta/seurat/dalia_merged_bc.rds")

bc@data[is.nan(bc@data)] <- 0
bc@normalized[is.nan(bc@normalized)] <- 0
bc@scaled[is.nan(bc@scaled)] <- 0

mat <- bc@normalized

mat <- t(mat)
mat <- scale(x = mat, center = TRUE, scale = TRUE)

sv <- svd(x = mat, nu nu = 50, nv = 50)