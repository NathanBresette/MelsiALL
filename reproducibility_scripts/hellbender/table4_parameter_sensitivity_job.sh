#!/bin/bash
#SBATCH --job-name=melsi_sensitivity
#SBATCH --output=table4_sensitivity_%A_%a.out
#SBATCH --error=table4_sensitivity_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-275

# Load R module (Hellbender uses lowercase 'r')
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Calculate total number of simulations
# 25 replications × 11 parameter combinations (6 B values + 5 m_frac values) = 275 total
# Each replication generates ONE dataset and tests all parameters on it (matches original paper)
TOTAL_SIMS=275

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}

# Run the R script in parallel mode
cd $SLURM_SUBMIT_DIR
Rscript table4_parameter_sensitivity.R $SIM_INDEX $TOTAL_SIMS
