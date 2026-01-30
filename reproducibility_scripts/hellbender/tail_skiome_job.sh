#!/bin/bash
# Tail SKIOME omnibus job output in real-time

JOB_ID=${1:-12402757}  # Use provided job ID or default to current job

echo "Tailing output for job $JOB_ID..."
echo "Press Ctrl+C to stop"
echo ""

ssh nbhtd@hellbender.rnet.missouri.edu "tail -f ~/melsi_simulations/hellbender/skiome_omnibus_${JOB_ID}.out"
