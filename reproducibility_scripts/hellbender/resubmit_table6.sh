#!/bin/bash
# Find and resubmit failed Table 6 tasks

cd ~/melsi_simulations/hellbender

failed_tasks=""
for i in {1..150}; do
    if [ ! -f "table6_sim_${i}.csv" ]; then
        failed_tasks="${failed_tasks}${i},"
    fi
done

failed_tasks=${failed_tasks%,}

if [ -n "$failed_tasks" ]; then
    echo "Failed tasks: $failed_tasks"
    echo "Resubmitting..."
    sbatch --array=$failed_tasks table6_feature_correlation_job.sh
else
    echo "All 150 tasks completed successfully!"
fi
