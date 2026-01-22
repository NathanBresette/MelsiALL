#!/bin/bash
# ==============================================================================
# Hellbender Setup and Check Script
# ==============================================================================
# This script checks your environment and helps set up for running simulations
# ==============================================================================

echo "================================================================================"
echo "MeLSI Hellbender Setup and Environment Check"
echo "================================================================================"
echo ""

# Check if we're on Hellbender
if [[ $(hostname) != *"hellbender"* ]]; then
    echo "⚠️  WARNING: This script should be run on Hellbender!"
    echo "   Current hostname: $(hostname)"
    echo ""
fi

# Check R availability
echo "1. Checking R availability..."
if command -v R &> /dev/null; then
    R_VERSION=$(R --version | head -n 1)
    echo "   ✅ R found: $R_VERSION"
else
    echo "   ❌ R not found. Attempting to load R module..."
    # Try common R module names (Hellbender uses lowercase 'r')
    if module load r/4.4.0 2>/dev/null; then
        echo "   ✅ Loaded r/4.4.0"
    elif module load r/4.3.0 2>/dev/null; then
        echo "   ✅ Loaded r/4.3.0"
    elif module load r 2>/dev/null; then
        echo "   ✅ Loaded r"
    elif module load R/4.3.0 2>/dev/null; then
        echo "   ✅ Loaded R/4.3.0"
    elif module load R 2>/dev/null; then
        echo "   ✅ Loaded R"
    else
        echo "   ⚠️  Could not load R module automatically."
        echo "   Checking available R modules..."
        R_MODULES=$(module avail 2>&1 | grep -iE '^[[:space:]]*[rR]/' | head -5)
        if [ -n "$R_MODULES" ]; then
            echo "   Available R modules:"
            echo "$R_MODULES" | sed 's/^/     /'
            FIRST_R=$(echo "$R_MODULES" | head -1 | awk '{print $1}')
            echo "   Attempting to load: $FIRST_R"
            module load $FIRST_R 2>/dev/null
        else
            echo "   ❌ No R modules found. R may need to be installed."
            echo "   Check Hellbender documentation or contact HPC support."
        fi
        echo "   Attempting to continue anyway..."
    fi
    
    # Verify R is now available
    if command -v R &> /dev/null; then
        R_VERSION=$(R --version | head -n 1)
        echo "   ✅ R is now available: $R_VERSION"
    else
        echo "   ⚠️  R still not available. Jobs may fail."
        echo "   You may need to install R or request it from HPC support."
    fi
fi
echo ""

# Check required R packages
echo "2. Checking required R packages..."
Rscript -e "
required <- c('vegan', 'microbiome', 'GUniFrac', 'ape', 'MeLSI')
missing <- c()
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
    cat('   ❌ Missing:', pkg, '\n')
  } else {
    cat('   ✅ Found:', pkg, '\n')
  }
}
if (length(missing) > 0) {
  cat('\n⚠️  Missing packages:', paste(missing, collapse=', '), '\n')
  cat('   Install with: install.packages(c(', paste0(\"'\", missing, \"'\", collapse=', '), '))\n')
  quit(status=1)
}
" || {
    echo "   ⚠️  Some packages may be missing. Check output above."
    echo ""
}
echo ""

# Check SLURM availability
echo "3. Checking SLURM availability..."
if command -v sbatch &> /dev/null; then
    echo "   ✅ SLURM found: $(sbatch --version 2>&1 | head -n 1)"
else
    echo "   ❌ SLURM not found. This is required for job submission."
    exit 1
fi
echo ""

# Check if scripts exist
echo "4. Checking required scripts..."
SCRIPTS=("table1_type1_error.R" "table2_power_analysis.R" 
         "table1_type1_error_job.sh" "table2_power_analysis_job.sh"
         "combine_table1_results.R" "combine_table2_results.R")
for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        echo "   ✅ Found: $script"
    else
        echo "   ❌ Missing: $script"
        exit 1
    fi
done
echo ""

# Check if scripts are executable
echo "5. Making scripts executable..."
chmod +x table1_type1_error_job.sh table2_power_analysis_job.sh
chmod +x combine_table1_results.R combine_table2_results.R
echo "   ✅ Scripts are executable"
echo ""

# Check current directory
echo "6. Current working directory:"
echo "   $(pwd)"
echo ""

# Summary
echo "================================================================================"
echo "Setup Check Complete!"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "  1. Submit Type I error jobs:    sbatch table1_type1_error_job.sh"
echo "  2. Submit Power analysis jobs:  sbatch table2_power_analysis_job.sh"
echo "  3. Monitor jobs:                squeue -u \$USER"
echo "  4. After completion, combine:    Rscript combine_table1_results.R"
echo "                                    Rscript combine_table2_results.R"
echo ""
