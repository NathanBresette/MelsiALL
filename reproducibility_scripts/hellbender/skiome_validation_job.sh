#!/bin/bash
#SBATCH --job-name=skiome_melsi
#SBATCH --output=skiome_validation_%j.out
#SBATCH --error=skiome_validation_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Load R module
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Set number of threads for R (if packages support it)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

# Run the validation script
cd $SLURM_SUBMIT_DIR
Rscript skiome_validation.R
