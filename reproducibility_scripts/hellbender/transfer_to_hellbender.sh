#!/bin/bash
# ==============================================================================
# Transfer Script - Run from your LOCAL machine
# ==============================================================================
# This script transfers all necessary files to Hellbender
# ==============================================================================

HELLBENDER_USER="nbhtd"
HELLBENDER_HOST="hellbender.rnet.missouri.edu"
REMOTE_DIR="~/melsi_simulations/hellbender"

echo "================================================================================"
echo "Transferring MeLSI Simulation Files to Hellbender"
echo "================================================================================"
echo ""
echo "Target: ${HELLBENDER_USER}@${HELLBENDER_HOST}:${REMOTE_DIR}"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Transferring files from: $SCRIPT_DIR"
echo ""

# Create remote directory
echo "1. Creating remote directory..."
ssh ${HELLBENDER_USER}@${HELLBENDER_HOST} "mkdir -p ${REMOTE_DIR}"
echo "   ✅ Remote directory created"
echo ""

# Transfer all files in hellbender directory
echo "2. Transferring files..."
rsync -avz --progress \
    table1_type1_error.R \
    table2_power_analysis.R \
    table1_type1_error_job.sh \
    table2_power_analysis_job.sh \
    combine_table1_results.R \
    combine_table2_results.R \
    setup_and_check.sh \
    HELLBENDER_SETUP.md \
    ${HELLBENDER_USER}@${HELLBENDER_HOST}:${REMOTE_DIR}/

if [ $? -eq 0 ]; then
    echo "   ✅ Files transferred successfully"
else
    echo "   ❌ Transfer failed"
    exit 1
fi
echo ""

echo "================================================================================"
echo "Transfer Complete!"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "  1. SSH into Hellbender:"
echo "     ssh ${HELLBENDER_USER}@${HELLBENDER_HOST}"
echo ""
echo "  2. Navigate to directory:"
echo "     cd ${REMOTE_DIR}"
echo ""
echo "  3. Run setup check:"
echo "     bash setup_and_check.sh"
echo ""
echo "  4. Submit jobs:"
echo "     sbatch table1_type1_error_job.sh"
echo "     sbatch table2_power_analysis_job.sh"
echo ""
