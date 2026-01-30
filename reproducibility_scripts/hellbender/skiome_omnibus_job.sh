#!/bin/bash
#SBATCH --job-name=skiome_omnibus
#SBATCH --output=skiome_omnibus_%j.out
#SBATCH --error=skiome_omnibus_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G

# Load R module (Hellbender uses lowercase 'r')
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

# Set number of threads for R (if packages support it)
export OMP_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6

# Run the omnibus-only validation script
cd $SLURM_SUBMIT_DIR
Rscript skiome_omnibus_only.R
