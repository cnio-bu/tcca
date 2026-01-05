#!/bin/bash

## Set directorires
rscript_dir="6_tsv_to_bpcells_lvl2.R"
data_dir="/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/cnv_cells_genes"

## Get list of files' full paths
files=($(find "${data_dir}" -type f))

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file}"

    #Run command
    sbatch -c 10 -o log.txt -e error.txt --mem=100G -t200 --wrap "$r_command"

done