#!/bin/bash

## Set directorires
rscript_dir="scevan_allsamples_2ndrun.R"

## Get list of files' full paths
files=("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/adrenalnb_rui_chong.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/bcc_catherine_dyao.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/brca_bhupinder_pal.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/brmets_jana_biermann.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/cell_lines_gabriella_kinker.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/nsclc_stefan_salcher.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/pancancer_sunny_wu.rds" "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw/urothelial_chen.rds")

## Loop through files
for file in "${files[@]}"; do
    #Set command
    r_command="Rscript ${rscript_dir} ${file}"

    #Run command
    echo "Running SCEVAN for ${file}..."
    sbatch -c 6 -o log.txt -e error.txt --mem=300G -plong -t10000 --wrap "$r_command"

done