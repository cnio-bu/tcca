#!/bin/bash
#SBATCH --job-name=sctherapy_mp
#SBATCH --array=1-611%100
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=logs_sctherapy/mp_%a.out
#SBATCH --error=logs_sctherapy/mp_%a.err
#SBATCH --exclude=hm001


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

echo "Processing sample: $SAMPLE"
Rscript run_sctherapy_mps.R "$SAMPLE"