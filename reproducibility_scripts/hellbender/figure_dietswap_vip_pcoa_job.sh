#!/bin/bash
#SBATCH --job-name=dietswap_fig
#SBATCH --output=dietswap_fig_%j.out
#SBATCH --error=dietswap_fig_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Load R module (Hellbender uses lowercase 'r')
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null || echo "Warning: Could not load R module"

# Set R library path to include user-installed packages
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Run the figure generation script
cd $SLURM_SUBMIT_DIR
Rscript figure_dietswap_vip_pcoa.R
