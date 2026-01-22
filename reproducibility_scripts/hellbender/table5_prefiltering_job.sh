#!/bin/bash
#SBATCH --job-name=melsi_prefilter
#SBATCH --output=table5_prefilter_%A_%a.out
#SBATCH --error=table5_prefilter_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-150

# Load R module (Hellbender uses lowercase 'r')
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Calculate total number of simulations
# 3 effect sizes × 50 simulations = 150 total
TOTAL_SIMS=150

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}

# Run the R script in parallel mode
cd $SLURM_SUBMIT_DIR
Rscript table5_prefiltering.R $SIM_INDEX $TOTAL_SIMS
