#!/bin/bash

## Set directorires
rscript_dir="filter_cellxgene_tables.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes"

## Get list of files' full paths
files=($(find "${data_dir}" -type f))

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file}"

    #Run command
    sbatch -c 10 -o log.txt -e error.txt --mem=80G -t120 --wrap "$r_command"

done