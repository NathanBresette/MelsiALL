#!/bin/bash
# ==============================================================================
# Kill Current SKIOME Job and Resubmit via SLURM
# ==============================================================================
# This script kills any running SKIOME processes and resubmits via SLURM
# ==============================================================================

HELLBENDER_USER="nbhtd"
HELLBENDER_HOST="hellbender.rnet.missouri.edu"
REMOTE_DIR="~/melsi_simulations/hellbender"

echo "================================================================================"
echo "Killing Current SKIOME Job and Resubmitting via SLURM"
echo "================================================================================"
echo ""

# Commands to run on Hellbender
COMMANDS="cd ${REMOTE_DIR} && \
echo 'Step 1: Finding and killing running SKIOME SLURM jobs...' && \
JOB_IDS=\$(squeue -u \$USER -n skiome_melsi --format='%i' --noheader) && \
if [ -n \"\$JOB_IDS\" ]; then \
  echo \"Found job IDs: \$JOB_IDS\" && \
  for JOB_ID in \$JOB_IDS; do \
    echo \"Cancelling job \$JOB_ID...\" && \
    scancel \$JOB_ID && \
    sleep 1; \
  done && \
  echo 'Jobs cancelled. Waiting for cleanup...' && \
  sleep 3; \
else \
  echo 'No running SKIOME jobs found.'; \
fi && \
echo '' && \
echo 'Step 2: Checking for any remaining processes...' && \
ps aux | grep 'Rscript.*skiome' | grep -v grep || echo 'No processes found' && \
echo '' && \
echo 'Step 3: Submitting SKIOME job via SLURM with updated settings (15h, 6 CPUs)...' && \
JOB_ID=\$(sbatch skiome_validation_job.sh | awk '{print \$4}') && \
echo 'Job submitted with ID:' \$JOB_ID && \
echo '' && \
echo 'Step 4: Checking job status...' && \
sleep 2 && \
squeue -j \$JOB_ID && \
echo '' && \
echo 'Step 5: Job output files will be:' && \
echo '  - skiome_validation_\${JOB_ID}.out' && \
echo '  - skiome_validation_\${JOB_ID}.err'"

echo "Connecting to Hellbender and executing commands..."
echo "You will be prompted for your password."
echo ""

ssh "${HELLBENDER_USER}@${HELLBENDER_HOST}" "$COMMANDS"

echo ""
echo "================================================================================"
echo "Done! Monitor the job with:"
echo "  ssh ${HELLBENDER_USER}@${HELLBENDER_HOST}"
echo "  cd ${REMOTE_DIR}"
echo "  squeue -u \$USER"
echo "  tail -f skiome_validation_*.err"
echo "================================================================================"
