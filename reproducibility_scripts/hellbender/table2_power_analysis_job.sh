#!/bin/bash
#SBATCH --job-name=melsi_power
#SBATCH --output=table2_power_%A_%a.out
#SBATCH --error=table2_power_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-450

# Load R module (Hellbender uses lowercase 'r')
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Calculate total number of simulations
# 3 effect sizes × 3 sample sizes × 50 simulations = 450 total
TOTAL_SIMS=450

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}

# Run the R script in parallel mode
cd $SLURM_SUBMIT_DIR
Rscript table2_power_analysis.R $SIM_INDEX $TOTAL_SIMS
