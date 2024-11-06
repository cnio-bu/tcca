import pandas as pd
import numpy as np
from pycirclize import Circos
import matplotlib.pyplot as plt
import os

os.chdir("/home/lmgonzalezb/Documents/bc-meta/")

# Load data
metadata = pd.read_csv("tcca_annotation_raw.tsv", sep="\t", index_col=0)

# Create new columns: study_patient and study_sample
metadata["study_patient"] = metadata["study"] + "_" + metadata["patient"]
metadata["study_sample"] = metadata["study"] + "_" + metadata["sample"]

# Correct "KIRCH" by "KIRC" in tumor type
metadata["tumor_type"] = np.where(
    metadata["tumor_type"] == "KIRCH", "KIRC", metadata["tumor_type"]
)

#  One neck skin cancer sample is a lymph_node metastasis.
metadata["refined_tumor_site"] = np.where(
    metadata["tumor_site"] == "Neck", "lymph_node", metadata["refined_tumor_site"]
)

# One ovarian cancer sample is an endometrial cancer sample that metastasized from endometrium
# (primary cancer is endometrial cancer).
metadata["tumor_type"] = np.where(
    (metadata["tumor_type"] == "OV")
    & (metadata["sample_type"] == "m")
    & (metadata["tumor_site"] == "ovary"),
    "UCEC",
    metadata["tumor_type"],
)


# Create a new column with tumor types names for patients and 'ccl' label for cell lines.
metadata["tumor_type2"] = np.where(
    metadata["patient"] == "ccl", "ccl", metadata["tumor_type"]
)


# Create cancer type names for the plot
cancer_type = {
    "Brain Cancer": ["GBM", "MB", "OGD"],
    "Neuroblastic Tumors": ["GNB", "NB"],
    "Blood Cancer": ["ALL", "LAML", "CLL", "MM"],
    "Skin Cancer": ["BCC", "SKCM", "SKSC", "SKAM", "UVM"],
    "Sarcoma/Soft Tissue Cancer": ["SARC", "GIST"],
    "Breast Cancer": ["BRCA"],
    "Lung Cancer": ["SCLC", "NSCLC", "LUAD", "LUSC", "LCLC", "PLEU"],
    "Ovarian Cancer": ["OV"],
    "Colon/Colorectal Cancer": ["COAD", "READ"],
    "Endometrial/Uterine Cancer": ["CESC", "UCEC", "UCS"],
    "Liver/Biliary Cancer": ["LIHC", "CHOL"],
    "Bladder Cancer": ["BLCA"],
    "Head and Neck Cancer": ["HNSC"],
    "Prostate Cancer": ["PRAD"],
    "Kidney Cancer": ["KRCC", "KTCC", "KIRC"],
    "Esophageal Cancer": ["ESCA", "ESCC"],
    "Pancreatic Cancer": ["PAAD"],
    "Miscellaneous Cancer": ["MISC"],
    "Cell lines": ["ccl"],
}

# Convert cancer_type_broad dictionary to a dataframe
cancer_type_broad_df = pd.DataFrame(
    [(k, v) for k, lst in cancer_type.items() for v in lst],
    columns=["cancer_type_broad", "tumor_type2"],
)

# Merge metadata with cancer_type_broad_df
metadata = metadata.merge(cancer_type_broad_df, on="tumor_type2", how="left")

# Group by cancer_type_broad and summarize the data
circlize_data = metadata.groupby("cancer_type_broad").agg(
    cell_count=("cancer_type_broad", "size"),
    patient_count=(
        "study_sample",
        lambda x: x[metadata["patient"] != "ccl"].nunique(dropna=False),
    ),
    cell_line_count=(
        "study_sample",
        lambda x: x[metadata["patient"] == "ccl"].nunique(dropna=False),
    ),
    primary_count=(
        "study_sample",
        lambda x: x[metadata["sample_type"] == "p"].nunique(dropna=False),
    ),
    metastasis_count=(
        "study_sample",
        lambda x: x[metadata["sample_type"] == "m"].nunique(dropna=False),
    ),
    unknown_count=(
        "study_sample",
        lambda x: x[~metadata["sample_type"].isin(["p", "m"])].nunique(dropna=False),
    ),
    treated_count=(
        "study_sample",
        lambda x: x[metadata["treated"] == "t"].nunique(dropna=False),
    ),
    untreated_count=(
        "study_sample",
        lambda x: x[metadata["treated"] == "f"].nunique(dropna=False),
    ),
    unknown_treat_count=(
        "study_sample",
        lambda x: x[~metadata["treated"].isin(["t", "f"])].nunique(dropna=False),
    ),
    sample_count=("study_sample", "nunique"),
    total_patient_count=("study_patient", "nunique"),
)

# Compute primary and metastatic number of samples as percentage.
circlize_data[["primary_count", "metastasis_count", "unknown_count"]] = (
    circlize_data[["primary_count", "metastasis_count", "unknown_count"]].div(
        circlize_data["sample_count"], axis=0
    )
    * 100
).round(2)

# Compute primary and metastatic number of samples as percentage.
circlize_data[["treated_count", "untreated_count", "unknown_treat_count"]] = (
    circlize_data[["treated_count", "untreated_count", "unknown_treat_count"]].div(
        circlize_data["sample_count"], axis=0
    )
    * 100
).round(2)

# Rename tumor sites.
translat_human_sites = {
    "adrenal_gland": "Adrenal gland",
    "bladder": "Bladder",
    "bone_marrow": "Bone marrow",
    "brain": "Brain",
    "breast": "Breast",
    "colon": "Colon",
    "skin": "Skin",
    "esophagus": "Esophagus",
    "oesophagus": "Esophagus",
    "kidney": "Kidney",
    "liver": "Liver",
    "lung": "Lung",
    "lymph_node": "Lymph node",
    "other": "Other",
    "ovary": "Ovary",
    "pancreas": "Pancreas",
    "prostate": "Prostate",
    "soft_tissue": "Soft tissue",
    "endometrium": "Endometrium",
    "rectum": "Rectum",
    "cervix": "Cervix",
    "bone": "Bone",
    "oral_cavity": "Oral cavity",
}

metadata["tumor_site_broad"] = np.where(
    metadata["refined_tumor_site"].isin(translat_human_sites.keys()),
    metadata["refined_tumor_site"],
    "other",
)

metadata["tumor_site_broad"] = metadata["tumor_site_broad"].map(translat_human_sites)

# Get number of samples from each of the tumor sites.
# Group by tumor_site_broad and cancer_type_broad to calculate number of samples per site
tumor_site_counts = (
    metadata.groupby(["cancer_type_broad", "tumor_site_broad"])
    .agg(sample_count_per_site=("study_sample", "nunique"))
    .reset_index()
)

# Pivot the tumor_site_broad into separate columns
tumor_site_counts_pivot = tumor_site_counts.pivot_table(
    index="cancer_type_broad",
    columns="tumor_site_broad",
    values="sample_count_per_site",
    fill_value=0,
)

# Transform tumor site counts to percentage.
tumor_site_counts_pivot = round(
    tumor_site_counts_pivot.div(tumor_site_counts_pivot.sum(axis=1), axis=0) * 100, 2
)

# Join tumor site counts with samples counts in circlize_data
circlize_data = circlize_data.join(tumor_site_counts_pivot)

# Add a total_count column
# circlize_data['total_count'] = circlize_data['patient_count'] + circlize_data['cell_line_count']

# Sort by total_count in descending order
circlize_data = circlize_data.sort_values(
    by="sample_count", ascending=False
).reset_index()

# Place the cell lines at the end of the circular plot
circlize_data = pd.concat(
    [circlize_data.iloc[1:], circlize_data.iloc[[0]]], ignore_index=True
)

circlize_data = circlize_data.set_index("cancer_type_broad")

# Create a dataframe only with patient samples or ccl.
# circlize_data_pat = circlize_data.iloc[0:len(circlize_data) - 1]
# circlize_data_ccl = circlize_data.iloc[[len(circlize_data) - 1]]

# Save the table
circlize_data.to_csv("input_circlize.tsv", sep="\t", index=False, header=True)


# Initialize circos sectors
sectors = circlize_data.index.tolist()
sectors = dict(zip(sectors, [1] * len(sectors)))
circos = Circos(sectors, space=2, start=0, end=180)


for i in range(0, len(circlize_data.index)):
    sector = circos.sectors[i]
    if circlize_data.index[i] == "Cell lines":
        sample_max = circlize_data.loc["Cell lines", "sample_count"].tolist()
        # cell_max = circlize_data.loc["Cell lines", "cell_count"].tolist()
        y_ticks_sample = [0, sample_max]
        # y_ticks_cell = [0, cell_max]
    else:
        sample_max = max(circlize_data.drop("Cell lines")["sample_count"])
        # cell_max = max(circlize_data.drop("Cell lines")["sample_count"])
        y_ticks_sample = [0, sample_max]
        # y_ticks_cell = [0, cell_max]

    # Plot cell lines and patient number per cancer type.
    # track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
    # track1.stacked_bar(
    #     circlize_data.loc[[sector.name], ["patient_count", "cell_line_count"]],
    #     width=0.6,
    #     cmap={"patient_count": "#427394", "cell_line_count": "#F78C1F"},
    #     vmax=sample_max,
    #     show_label=False,
    # )
    # track1.axis()
    # y_ticks = np.linspace(0, 150, 6, dtype=int)
    # y_labels = list(map(str, y_ticks))
    # track1.yticks(
    #     y_ticks, labels=None, vmax=sample_max, side="left"
    # )
    # Plot primary, metastatic or unknown sample number per cancer type.
    track1 = sector.add_track((85, 100), r_pad_ratio=0.1)
    x = np.arange(sector.start, sector.end) + 0.5
    y = [circlize_data.loc[sector.name, "sample_count"]]
    track1.bar(
        x, y, color="#740090", width=0.6, vmax=max(circlize_data["sample_count"])
    )
    track1.axis()
    y_ticks = np.linspace(0, max(circlize_data["sample_count"]), 8, dtype=int)
    y_labels = list(map(str, y_ticks))
    track1.yticks(
        y_ticks, labels=None, vmax=max(circlize_data["sample_count"]), side="left"
    )

    # Plot primary, metastatic or unknown sample number per cancer type.
    track2 = sector.add_track((69, 84), r_pad_ratio=0.1)
    track2.stacked_bar(
        circlize_data.loc[
            [sector.name], ["primary_count", "metastasis_count", "unknown_count"]
        ],
        width=0.6,
        cmap={
            "primary_count": "#F0BFD0",
            "metastasis_count": "#C10044",
            "unknown_count": "#808080",
        },
        bar_kws=dict(ec="black", lw=0),
        vmax=100.01,
        show_label=False,
    )
    track2.axis()
    y_ticks = np.linspace(0, 100, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track2.yticks(y_ticks, labels=None, vmax=100.01, side="left")

    # Plot untreated, treated or unknown sample number per cancer type.
    track3 = sector.add_track((53, 68), r_pad_ratio=0.1)
    track3.stacked_bar(
        circlize_data.loc[
            [sector.name], ["untreated_count", "treated_count", "unknown_treat_count"]
        ],
        width=0.6,
        cmap={
            "untreated_count": "#D18B6E",
            "treated_count": "#6ED1BC",
            "unknown_treat_count": "#808080",
        },
        bar_kws=dict(ec="black", lw=0),
        vmax=100.01,
        show_label=False,
    )
    track3.axis()
    y_ticks = np.linspace(0, 100, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track3.yticks(y_ticks, labels=None, vmax=100.01, side="left")

    # Plot untreated, treated or unknown sample number per cancer type.
    track4 = sector.add_track((37, 52), r_pad_ratio=0.1)
    track4.stacked_bar(
        circlize_data.loc[
            [sector.name], metadata["tumor_site_broad"].unique().tolist()
        ],
        width=0.6,
        cmap={
            "Bone marrow": "#FFE072",
            "Brain": "#CA9A8C",
            "Adrenal gland": "#FF9A4A",
            "Breast": "#FF8FAB",
            "Skin": "#5E2D2C",
            "Lung": "#D3C3E0",
            "Soft tissue": "#FB6467",
            "Esophagus": "#5780FE",
            "Bladder": "#58B368",
            "Lymph node": "#B6F884",
            "Liver": "#309975",
            "Pancreas": "#B47EB3",
            "Ovary": "#9C0D38",
            "Prostate": "#2997CF",
            "Colon": "#005D95",
            "Kidney": "#918050",
            "Endometrium": "#F6A38B",
            "Rectum": "#9EDDF9",
            "Cervix": "#E7F20A",
            "Bone": "#FFED42",
            "Oral cavity": "#97D1A9",
            "Other": "#BBB9B7",
        },
        bar_kws=dict(ec="black", lw=0),
        vmax=100.01,
        show_label=False,
    )
    track4.axis()
    y_ticks = np.linspace(0, 100, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track4.yticks(y_ticks, labels=None, vmax=100.01, side="left")

    # Plot cell number per cancer type.
    track5 = sector.add_track((21, 36), r_pad_ratio=0.1)
    x = np.arange(sector.start, sector.end) + 0.5
    y = [circlize_data.loc[sector.name, "cell_count"]]



    track5.bar(x, y, color="#427394", width=0.6, vmax=max(circlize_data["cell_count"]))
    track5.axis()
    y_ticks = np.linspace(0, 350000, 8, dtype=int)
    y_labels = list(map(str, y_ticks))
    track5.yticks(y_ticks, labels=None, vmax=350000, side="left")

fig = circos.plotfig()
plt.savefig("cohort_circular.pdf", format="pdf", bbox_inches="tight")

plt.show()
