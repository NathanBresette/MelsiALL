#!/bin/bash
# Find missing Table 6 simulations and create a job script to resubmit them

# Find missing files in High correlation range (151-200)
missing_sims=""
for i in {151..200}; do
    if [ ! -f "table6_sim_${i}.csv" ]; then
        missing_sims="${missing_sims}${i} "
    fi
done

echo "Missing simulations: $missing_sims"
echo "Count: $(echo $missing_sims | wc -w)"

# Create a resubmit job script for missing simulations
if [ -n "$missing_sims" ]; then
    cat > resubmit_table6_missing.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=melsi_corr_missing
#SBATCH --output=table6_corr_missing_%A_%a.out
#SBATCH --error=table6_corr_missing_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=MISSING_ARRAY

# Load R module - exit if it fails
module load r/4.4.0 2>/dev/null || module load r/4.3.0 2>/dev/null || module load r 2>/dev/null
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Module load failed."
    exit 1
fi

# Set R library path
export R_LIBS_USER=~/R:${R_LIBS_USER:-}

# Get the simulation index from SLURM array task ID
SIM_INDEX=${SLURM_ARRAY_TASK_ID}
TOTAL_SIMS=200

# Run the R script
cd $SLURM_SUBMIT_DIR
Rscript table6_feature_correlation.R $SIM_INDEX $TOTAL_SIMS
EOF
    
    # Replace MISSING_ARRAY with actual array indices
    missing_array=$(echo $missing_sims | tr ' ' ',')
    sed -i "s/MISSING_ARRAY/$missing_array/" resubmit_table6_missing.sh
    
    echo "Created resubmit_table6_missing.sh"
    echo "Array indices: $missing_array"
fi
