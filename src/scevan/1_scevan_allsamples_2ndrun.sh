#!/bin/bash

## Set directorires
rscript_dir="1_scevan_allsamples.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/"

## Get list of files' full paths
files=("urothelial_chen_v5")

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file} 5"

    #Run command
    echo "Running SCEVAN for ${file}..."
    sbatch -c 20 -o log.txt -e error.txt --mem=500G -t1440 --wrap "$r_command"

done