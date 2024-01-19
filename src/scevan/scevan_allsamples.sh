#!/bin/bash

## Set directorires
rscript_dir="scevan_allsamples.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw"

## Get list of files' full paths
files=($(find "${data_dir}" -type f))

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file}"

    #Run command
    echo "Running SCEVAN for ${file}..."
    sbatch -c 20 -o log.txt -e error.txt --mem=300G -t1440 --wrap "$r_command"

done