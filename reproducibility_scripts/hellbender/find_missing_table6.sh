#!/bin/bash
# Find missing Table 6 simulation files

cd ~/melsi_simulations/hellbender

echo "=== FINDING MISSING HIGH CORRELATION FILES (151-200) ==="
missing=""
for i in {151..200}; do
    if [ ! -f "table6_sim_${i}.csv" ]; then
        missing="${missing}${i} "
        echo "Missing: table6_sim_${i}.csv"
    fi
done

echo ""
echo "Missing count: $(echo $missing | wc -w)"
echo "Missing indices: $missing"

# Create array string for SLURM
if [ -n "$missing" ]; then
    missing_array=$(echo $missing | tr ' ' ',')
    echo ""
    echo "SLURM array string: $missing_array"
fi
