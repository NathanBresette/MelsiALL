#!/bin/bash
# ==============================================================================
# Quick Start Script - Run on Hellbender after transferring files
# ==============================================================================
# This script does everything needed to start the simulations
# ==============================================================================

echo "================================================================================"
echo "MeLSI Simulations Quick Start"
echo "================================================================================"
echo ""

# Run setup check
echo "Running setup check..."
bash setup_and_check.sh

SETUP_EXIT=$?

# If R is the issue, try to find it
if [ $SETUP_EXIT -ne 0 ]; then
    echo ""
    echo "Setup check had issues. Attempting to find R..."
    bash find_and_setup_r.sh
    
    # Re-run setup check after R setup attempt
    echo ""
    echo "Re-running setup check..."
    bash setup_and_check.sh
    SETUP_EXIT=$?
fi

if [ $SETUP_EXIT -ne 0 ]; then
    echo ""
    echo "❌ Setup check failed. Please fix issues before proceeding."
    echo "   Common issues:"
    echo "   - R not available: Run 'bash find_and_setup_r.sh' for help"
    echo "   - Missing packages: Install with R and install.packages()"
    echo ""
    exit 1
fi

echo ""
echo "================================================================================"
echo "Submitting Jobs"
echo "================================================================================"
echo ""

# Submit Type I error job
echo "1. Submitting Type I Error Analysis job..."
JOB1=$(sbatch table1_type1_error_job.sh | awk '{print $4}')
if [ -n "$JOB1" ]; then
    echo "   ✅ Type I Error job submitted: $JOB1"
else
    echo "   ❌ Failed to submit Type I Error job"
    exit 1
fi

# Submit Power analysis job
echo "2. Submitting Power Analysis job..."
JOB2=$(sbatch table2_power_analysis_job.sh | awk '{print $4}')
if [ -n "$JOB2" ]; then
    echo "   ✅ Power Analysis job submitted: $JOB2"
else
    echo "   ❌ Failed to submit Power Analysis job"
    exit 1
fi

echo ""
echo "================================================================================"
echo "Jobs Submitted Successfully!"
echo "================================================================================"
echo ""
echo "Job IDs:"
echo "  Type I Error: $JOB1"
echo "  Power Analysis: $JOB2"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo "  bash monitor_jobs.sh"
echo ""
echo "After completion, combine results with:"
echo "  Rscript combine_table1_results.R"
echo "  Rscript combine_table2_results.R"
echo ""
