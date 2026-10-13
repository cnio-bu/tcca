#!/bin/bash
#SBATCH --job-name=mp_scoring
#SBATCH --array=1-611%100
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/mp_%a.out
#SBATCH --error=logs/mp_%a.err
#SBATCH --exclude=hm001

mkdir -p logs

# Read sample name based on array indexes
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_list.txt)

echo "Processing sample: $SAMPLE"
Rscript run_high_mp_cells.R "$SAMPLE"
