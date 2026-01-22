#!/bin/bash
# ==============================================================================
# Job Monitoring Script - Run on Hellbender
# ==============================================================================
# This script monitors the progress of your simulation jobs
# ==============================================================================

echo "================================================================================"
echo "MeLSI Simulation Job Monitor"
echo "================================================================================"
echo ""

# Check job status
echo "1. Current Job Status:"
echo "   $(squeue -u $USER -o '%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R' 2>/dev/null || echo '   No jobs found or squeue not available')"
echo ""

# Count completed simulations
echo "2. Simulation Progress:"
TABLE1_COUNT=$(ls -1 table1_sim_*.csv 2>/dev/null | wc -l | tr -d ' ')
TABLE2_COUNT=$(ls -1 table2_sim_*.csv 2>/dev/null | wc -l | tr -d ' ')

echo "   Type I Error simulations: ${TABLE1_COUNT}/300 completed"
echo "   Power Analysis simulations: ${TABLE2_COUNT}/270 completed"
echo ""

# Check for errors
echo "3. Checking for errors..."
ERROR_COUNT=$(ls -1 table1_type1_*_*.err table2_power_*_*.err 2>/dev/null | xargs grep -l "Error\|error\|ERROR" 2>/dev/null | wc -l | tr -d ' ')
if [ "$ERROR_COUNT" -gt 0 ]; then
    echo "   ⚠️  Found errors in $ERROR_COUNT error files"
    echo "   Check recent errors:"
    ls -t table1_type1_*_*.err table2_power_*_*.err 2>/dev/null | head -3 | while read errfile; do
        echo "     - $errfile (last 3 lines):"
        tail -3 "$errfile" | sed 's/^/       /'
    done
else
    echo "   ✅ No obvious errors found"
fi
echo ""

# Check if all jobs are complete
if [ "$TABLE1_COUNT" -eq 300 ] && [ "$TABLE2_COUNT" -eq 270 ]; then
    echo "================================================================================"
    echo "✅ All Simulations Complete!"
    echo "================================================================================"
    echo ""
    echo "Next steps:"
    echo "  1. Combine Type I Error results:"
    echo "     Rscript combine_table1_results.R"
    echo ""
    echo "  2. Combine Power Analysis results:"
    echo "     Rscript combine_table2_results.R"
    echo ""
    echo "  3. Download results to local machine:"
    echo "     scp table*_summary.csv nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/"
    echo ""
else
    echo "================================================================================"
    echo "Jobs Still Running..."
    echo "================================================================================"
    echo ""
    echo "Run this script again to check progress:"
    echo "  bash monitor_jobs.sh"
    echo ""
fi
