#!/bin/bash

# List of study names (36 in total)
studies=($(find /storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl1 -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sed "s/_v5//"))

# Loop through studies and submit jobs with specific memory allocation
for study in "${studies[@]}"; do
    if [[ "$study" == "nsclc_stefan_salcher"  || "$study" == "cell_lines_gabriella_kinker" ]]; then
        sbatch --job-name=nmf_subclone \
               --output=logs/nmf_%A_%a.out \
               --error=logs/nmf_%A_%a.err \
               --time=05:00:00 \
               --mem=150G \
               --cpus-per-task=1 \
               --wrap="Rscript ./compute_nmf_subclone.R $study"
    else
        sbatch --job-name=nmf_subclone \
               --output=logs/nmf_%A_%a.out \
               --error=logs/nmf_%A_%a.err \
               --time=05:00:00 \
               --mem=20G \
               --cpus-per-task=1 \
               --wrap="Rscript ./compute_nmf_subclone.R $study"
    fi
done