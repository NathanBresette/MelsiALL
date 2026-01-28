#!/bin/bash
#SBATCH --job-name=melsi_scalability
#SBATCH --output=table3_scalability_%A_%a.out
#SBATCH --error=table3_scalability_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-100

# Load R module (Hellbender uses lowercase 'r') - exit if it fails
module load r/4.4.0 2>/dev/null || module load r 2>/dev/null || {
    echo "ERROR: Failed to load R module"
    exit 1
}

# Verify Rscript is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found after module load"
    exit 1
fi

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Calculate total number of simulations
# 10 conditions (5 sample sizes + 5 taxa sizes) × 10 simulations = 100 total
TOTAL_SIMS=100

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}

# Run the R script in parallel mode
cd $SLURM_SUBMIT_DIR
Rscript table3_scalability.R $SIM_INDEX $TOTAL_SIMS
