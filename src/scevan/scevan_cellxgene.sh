#!/bin/bash

## Set directorires
rscript_dir="scevan_cellxgene.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/cna/"

## Get list of directories' full paths
dirs=($(find "${data_dir}" -maxdepth 1 -mindepth 1 -type d))

## Loop through directories
for dir in "${dirs[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${dir}"

    #Run command
    echo "Running for ${dir}..."
    sbatch -c 10 -o log.txt -e error.txt --mem=300G -t1440 --wrap "$r_command"

done
