#!/bin/bash
# ==============================================================================
# Monitor SKIOME Validation Job in Real-Time
# ==============================================================================
# Run this script to see live output from the Hellbender job
# ==============================================================================

HELLBENDER_USER="nbhtd"
HELLBENDER_HOST="hellbender.rnet.missouri.edu"
REMOTE_DIR="~/melsi_simulations/hellbender/"

echo "================================================================================"
echo "Monitoring SKIOME Validation Job"
echo "================================================================================"
echo ""
echo "This will show live output from the job running on Hellbender."
echo "Press Ctrl+C to stop monitoring."
echo ""
echo "Connecting to: ${HELLBENDER_USER}@${HELLBENDER_HOST}"
echo ""

# Monitor both output and error logs
ssh ${HELLBENDER_USER}@${HELLBENDER_HOST} "cd ${REMOTE_DIR} && tail -f skiome_run2.out skiome_run2.err 2>/dev/null || tail -f skiome_validation_direct.out skiome_validation_direct.err 2>/dev/null"
