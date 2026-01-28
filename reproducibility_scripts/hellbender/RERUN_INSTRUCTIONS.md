# Instructions for Rerunning Tables on Hellbender

## Summary of Changes

1. **Table 1 (Type I Error)**: Updated combine script to include all 5 traditional methods (Jaccard, Weighted UniFrac, Unweighted UniFrac). The individual simulation files already contain all methods, so we just need to rerun the combine script.

2. **Table 3 (Scalability)**: Updated script to include all 5 traditional methods instead of just 2. This requires a full rerun of all simulations.

3. **Table 6 (Correlation)**: Already exists in data files. Added to manuscript with ranks.

## Steps to Rerun

### Table 1: Rerun Combine Script Only

The individual simulation files already have all 5 methods. Just rerun the combine script:

```bash
cd /path/to/reproducibility_scripts/hellbender
sbatch rerun_table1_combine.sh
```

This will regenerate `table1_type1_error_summary.csv` with all 5 traditional methods.

### Table 3: Full Rerun Required

Since we updated the script to include all 5 methods, we need to rerun all simulations:

```bash
cd /path/to/reproducibility_scripts/hellbender
sbatch table3_scalability_job.sh
```

After all array jobs complete, combine the results:

```bash
Rscript combine_table3_results.R
```

This will generate `table3_scalability_summary.csv` with ranks (1/6 to 6/6) instead of just comparing to best traditional.

## Updated Files

- `table1_type1_error.R` - Already runs all 5 methods (no changes needed)
- `combine_table1_results.R` - Updated to include all 5 methods
- `table3_scalability.R` - Updated to include all 5 methods
- `combine_table3_results.R` - Updated to calculate ranks (1/6 to 6/6)
- `manuscript/MeLSI_Research_Paper_mSystems.md` - Updated Table 1, Table 3, and added Table 6
