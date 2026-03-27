#!/bin/bash
#SBATCH --job-name=skiome_fig
#SBATCH --output=skiome_fig_%j.out
#SBATCH --error=skiome_fig_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

export R_LIBS_USER=~/R:${R_LIBS_USER:-}

cd $SLURM_SUBMIT_DIR
Rscript figure_skiome_round3.R
