import scanpy as sc
import pandas as pd
import scvi
import torch
import pymde
import matplotlib.pyplot as plt
import os
import seaborn as sns

scvi.settings.seed = 120394
torch.set_float32_matmul_precision("high")

# Set working directory
os.chdir("/storage/scratch01/shared/projects/bc-meta/")

# Subset malignant cells
adata = sc.read_h5ad(
    "single_cell/seurat/tcca/tcca_annotated_clustered.h5ad"
)
malignants = adata[adata.obs["malignancy"] == True]
malignants.obs_names = [f"c{i}" for i in range(1, malignants.n_obs + 1)]

# Load UCell scores of metaprograms (MPs)
ucell_scores_mps = pd.read_csv("functional_nmf/mps_ucell_scores.tsv", sep=" ")
ucell_scores_mps.index = [f"c{i}" for i in range(1, len(ucell_scores_mps) + 1)]

# Add UCell scores as a new embedding
malignants.obsm["X_MPsignatures"] = ucell_scores_mps.values

# Compute neighbors and UMAP using MP signatures
sc.pp.neighbors(malignants, use_rep="X_MPsignatures", n_pcs=ucell_scores_mps.shape[1])
sc.tl.umap(malignants)

# Save AnnData object and plot UMAP by study
malignants.write("functional_nmf/malignant_umap_mps.h5ad")
sc.pl.umap(malignants, color="study")
plt.savefig("functional_nmf/figures/malignant_umap_mps.png", dpi=600,
            bbox_inches='tight')

# Combine UCell scores with AnnData metadata
malignants.obs = pd.concat([malignants.obs, ucell_scores_mps], axis=1)

sc.pl.umap(malignants, color=ucell_scores_mps.columns, ncols=5)
plt.savefig("functional_nmf/figures/malignant_umap_mps_scores.png",
            dpi=600, bbox_inches='tight')

# Define general cancer types
cancer_type = {
    "Brain Cancer": ["GBM", "MB", "OGD"],
    "Neuroblastic Tumors": ["GNB", "NB"],
    "Blood Cancer": ["ALL", "LAML", "CLL", "MM"],
    "Skin Cancer": ["BCC", "SKCM", "SKSC", "SKAM", "UVM"],
    "Sarcoma/Soft Tissue Cancer": ["SARC", "GIST", "MESO"],
    "Breast Cancer": ["BRCA"],
    "Lung Cancer": ["SCLC", "NSCLC", "LUAD", "LUSC", "LCLC", "PLEU"],
    "Ovarian Cancer": ["OV"],
    "Colon/Colorectal Cancer": ["COAD", "READ"],
    "Endometrial/Uterine Cancer": ["CESC", "UCEC", "UCS"],
    "Liver/Biliary Cancer": ["LIHC", "CHOL"],
    "Bladder Cancer": ["BLCA"],
    "Head and Neck Cancer": ["HNSC"],
    "Prostate Cancer": ["PRAD"],
    "Kidney Cancer": ["KRCC", "KTCC", "KIRC", "KIRCH"],
    "Esophageal Cancer": ["ESCA", "ESCC"],
    "Pancreatic Cancer": ["PAAD"],
    "Other": ["THCA", "STAD", "MISC"]
}

# Map tumor codes to cancer types
inverse_cancer_type = {
    code: cancer for cancer, codes in cancer_type.items() for code in codes
}
malignants.obs["broad_cancer_type"] = malignants.obs["tumor_type"].map(
    inverse_cancer_type
)
print(malignants.obs)

# Define colors for broad cancer types
broad_cancer_type_palette = {
    "Blood Cancer": "#A3181B", "Brain Cancer": "#B2509E",
    "Neuroblastic Tumors": "#F06616", "Breast Cancer": "#db447a",
    "Skin Cancer": "#5E2D2C", "Lung Cancer": "#158A88",
    "Sarcoma/Soft Tissue Cancer": "#e8d52c", "Esophageal Cancer": "#007EB5",
    "Bladder Cancer": "#367040", "Liver/Biliary Cancer": "#03543C",
    "Pancreatic Cancer": "#694E85", "Ovarian Cancer": "#E834EB",
    "Prostate Cancer": "#005D95", "Colon/Colorectal Cancer": "#a7495a",
    "Endometrial/Uterine Cancer": "#FAD2D9", "Head and Neck Cancer": "#97D1A9",
    "Kidney Cancer": "#918050", "Other": "#6b6363"
}

sc.pl.umap(malignants, color="broad_cancer_type",
           palette=broad_cancer_type_palette)
plt.savefig("functional_nmf/figures/umap_mps_scores_cancertype.png",
            dpi=600, bbox_inches='tight')

# Subset metastatic skin cancer cells and plot violin plots
malignants_brca = malignants[
    (malignants.obs["broad_cancer_type"] == "Skin Cancer") & 
    (malignants.obs["sample_type"] == "m")
]
malignants_brca.obs.columns = malignants_brca.obs.columns.str.replace(
    r'_UCell$', '', regex=True
)
sc.pl.violin(malignants_brca, keys=[f"MP{i}" for i in range(1, 15)],
             ylabel="UCell score")
fig = plt.gcf()
fig.set_size_inches(12, 5)
plt.savefig("functional_nmf/figures/vln_mps_scores_skincancer_metastasis.png",
            dpi=600, bbox_inches='tight')

# Plot therapeutic clusters
tcs = pd.read_csv("/storage/scratch01/shared/projects/bc-meta/beyondcell/results/tcs.tsv",
                  sep="\t")
tcs.set_index("new_cell_id", inplace=True)
tcs = tcs.loc[:, "therapeutic_clusters_k.300.res.0.5"]

malignants.obs.columns = malignants.obs.columns.str.strip()
malignants.obs = malignants.obs.rename(
    columns={"therapeutic_clusters_k.300.res.0.5": "therapeutic_clusters"}
)
malignants.obs["therapeutic_clusters"] = malignants.obs["therapeutic_clusters"]\
    .astype("int").astype("category").astype("str")
malignants = malignants[malignants.obs["therapeutic_clusters"].notna()].copy()

tcs_colors = {
    "5": "#A3A500", "4": "#F8766D", "3": "#FFD64C",
    "2": "#00BF7D", "1": "#00B0F6", "0": "#E76BF3"
}

sc.pl.umap(malignants, color="therapeutic_clusters", palette=tcs_colors)
plt.savefig("functional_nmf/figures/umap_mps_scores_therapeutic_clusters.png",
            dpi=600, bbox_inches='tight')
