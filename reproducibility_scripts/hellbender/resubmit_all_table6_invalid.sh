#!/bin/bash
#SBATCH --job-name=melsi_corr_fix
#SBATCH --output=table6_corr_fix_%A_%a.out
#SBATCH --error=table6_corr_fix_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=2,3,4,6,7,8,10,11,12,14,15,16,18,19,20,22,23,24,26,27,29,30,32,33,35,36,38,39,41,42,44,45,47,48,50,51,52,53,55,56,57,59,60,61,63,64,65,67,68,69,71,72,73,75,76,77,79,80,81,83,84,85,87,88,89,91,92,93,95,96,97,99,100,101,102,103,104,106,107,109,110,112,113,114,116,117,118,120,121,122,124,125,126,128,129,130,132,133,134,136,137,138,140,141,142,144,145,146,148,149,150

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
