#!/bin/bash
#SBATCH --job-name=skiome_melsi
#SBATCH --output=skiome_validation_%j.out
#SBATCH --error=skiome_validation_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Load R module
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Run the validation script
cd $SLURM_SUBMIT_DIR
Rscript skiome_validation.R
