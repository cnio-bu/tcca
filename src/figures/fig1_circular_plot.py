import pandas as pd
import numpy as np
from pycirclize import Circos
import matplotlib.pyplot as plt
import os

os.chdir("/home/lmgonzalezb/Documents/bc-meta/")

# Load data
metadata = pd.read_csv("tcca_annotation_raw.tsv", sep="\t", index_col=0)

# Correct "KIRCH" by "KIRC" in tumor type
metadata["tumor_type"] = np.where(
    metadata["tumor_type"] == "KIRCH", "KIRC", metadata["tumor_type"]
)

# Create tumor site names for the plot
cancer_type = {
    "Brain Cancer": ["GBM", "MISC", "MB", "OGD"],
    "Neuroblastic Tumors": ["GNB", "NB"],
    "Blood Cancer": ["ALL", "LAML", "CLL", "MM"],
    "Skin Cancer": ["BCC", "SKCM", "SKSC", "SKAM"],
    "Sarcoma": ["SARC"],
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
    "Other": ["MISC", "UVM", "THCA", "STAD", "GIST"],
}


# Convert cancer_type_broad dictionary to a dataframe
cancer_type_broad_df = pd.DataFrame(
    [(k, v) for k, lst in cancer_type.items() for v in lst],
    columns=["cancer_type_broad", "tumor_type"],
)

# Merge metadata with cancer_type_broad_df
metadata = metadata.merge(cancer_type_broad_df, on="tumor_type", how="left")

# Create new columns: study_patient and study_sample
metadata["study_patient"] = metadata["study"] + "_" + metadata["patient"]
metadata["study_sample"] = metadata["study"] + "_" + metadata["sample"]

# Group by cancer_type_broad and summarize the data
circlize_data = (
    metadata.groupby("cancer_type_broad")
    .agg(
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
            lambda x: x[~metadata["sample_type"].isin(["p", "m"])].nunique(
                dropna=False
            ),
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
    )
    .reset_index()
)

# Add a total_count column
# circlize_data['total_count'] = circlize_data['patient_count'] + circlize_data['cell_line_count']

# Sort by total_count in descending order
circlize_data = circlize_data.sort_values(by="sample_count", ascending=False)

# Save the table
circlize_data.to_csv("input_circlize.tsv", sep="\t", index=False, header=True)

# Set cancer_type_broad as row index.
circlize_data.set_index("cancer_type_broad", inplace=True)


# Initialize circos sectors
sectors = circlize_data.index.tolist()
sectors = dict(zip(sectors, [1] * len(sectors)))
circos = Circos(sectors, space=2, start=0, end=180)



for sector in circos.sectors:
    # Plot cell lines and patient number per cancer type.
    track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
    track1.stacked_bar(
        circlize_data.loc[[sector.name], ["patient_count", "cell_line_count"]],
        width=0.6,
        cmap={"patient_count": "#427394", "cell_line_count": "#F78C1F"},
        vmax=max(circlize_data["sample_count"]),
        show_label=False,
    )
    track1.axis()
    y_ticks = np.linspace(0, 150, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track1.yticks(
        y_ticks, labels=None, vmax=max(circlize_data["sample_count"]), side="left"
    )

    # Plot primary, metastatic or unknown sample number per cancer type.
    track2 = sector.add_track((60, 79), r_pad_ratio=0.1)
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
        vmax=max(circlize_data["sample_count"]),
        show_label=False,
    )
    track2.axis()
    y_ticks = np.linspace(0, 150, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track2.yticks(
        y_ticks, labels=None, vmax=max(circlize_data["sample_count"]), side="left"
    )

    # Plot untreated, treated or unknown sample number per cancer type.
    track3 = sector.add_track((40, 59), r_pad_ratio=0.1)
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
        vmax=max(circlize_data["sample_count"]),
        show_label=False,
    )
    track3.axis()
    y_ticks = np.linspace(0, 150, 6, dtype=int)
    y_labels = list(map(str, y_ticks))
    track3.yticks(
        y_ticks, labels=None, vmax=max(circlize_data["sample_count"]), side="left"
    )

    # Plot cell number per cancer type.
    track4 = sector.add_track((20, 39), r_pad_ratio=0.1)
    x = np.arange(sector.start, sector.end) + 0.5
    y = [circlize_data.loc[sector.name, "cell_count"]]
    track4.bar(x, y, color="#3faa6d", width=0.6, vmax=max(circlize_data["cell_count"]))
    track4.axis()
    y_ticks = np.linspace(0, 350000, 8, dtype=int)
    y_labels = list(map(str, y_ticks))
    track4.yticks(y_ticks, labels=None, vmax=350000, side="left")

fig = circos.plotfig()
plt.savefig("cohort_circular.pdf", format="pdf", bbox_inches="tight")

plt.show()

