#!/bin/bash
#SBATCH --job-name=melsi_corr_missing
#SBATCH --output=table6_corr_missing_%A_%a.out
#SBATCH --error=table6_corr_missing_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=151,153,154,155,157,158,159,161,162,163,165,166,167,169,170,171,173,174,175,177,178,179,181,182,183,185,186,187,189,190,191,193,194,195,197,198,199

# Load R module - exit if it fails
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Module load failed."
    exit 1
fi

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Calculate total number of simulations
# 4 correlation levels × 50 simulations = 200 total
TOTAL_SIMS=200

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}

# Run the R script in parallel mode
cd $SLURM_SUBMIT_DIR
Rscript table6_feature_correlation.R $SIM_INDEX $TOTAL_SIMS
