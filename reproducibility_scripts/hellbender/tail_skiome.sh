#!/bin/bash
# Tail SKIOME job output and error files
# Usage: bash tail_skiome.sh [job_id]
# If job_id not provided, will use latest SKIOME job

JOB_ID=${1:-12356610}

echo "Tailing SKIOME job $JOB_ID..."
echo "Press Ctrl+C to stop"
echo ""
echo "=== Output (stdout) ==="
echo ""

ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && tail -f skiome_validation_${JOB_ID}.out skiome_validation_${JOB_ID}.err"
