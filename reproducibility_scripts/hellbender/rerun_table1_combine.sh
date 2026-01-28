#!/bin/bash
#SBATCH --job-name=combine_table1
#SBATCH --output=combine_table1_%j.out
#SBATCH --error=combine_table1_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Load R module
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null

# Set R library path
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Run combine script
cd $SLURM_SUBMIT_DIR
Rscript combine_table1_results.R
