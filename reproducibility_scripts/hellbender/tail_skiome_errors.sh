#!/bin/bash
# Tail only SKIOME job errors (stderr)
# Usage: bash tail_skiome_errors.sh [job_id]

JOB_ID=${1:-12356610}

echo "Tailing SKIOME job $JOB_ID errors..."
echo "Press Ctrl+C to stop"
echo ""

ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && tail -f skiome_validation_${JOB_ID}.err"
