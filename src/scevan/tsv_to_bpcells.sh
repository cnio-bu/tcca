#!/bin/bash

## Set directorires
rscript_dir="tsv_to_bpcells.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes_lvl3"

## Get list of files' full paths
files=($(find "${data_dir}" -type f))

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file}"

    #Run command
    sbatch -c 10 -o log.txt -e error.txt --mem=100G -t200 --wrap "$r_command"

done