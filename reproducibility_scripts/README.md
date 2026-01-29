# MeLSI Results Reproducibility Scripts

This folder contains scripts to reproduce all results from the MeLSI paper.

## Directory Structure

- **`hellbender/`** - All simulation scripts and job submission scripts for running on HPC
  - See `hellbender/HELLBENDER_GUIDE.md` for detailed instructions
  - Contains R scripts for all tables (Table 1-6) and SKIOME validation
  - Contains SLURM job scripts for parallel execution
- **`figure_atlas1006_vip_pcoa.R`** - Script to generate Atlas1006 VIP and PCoA figures

## Running Simulations

All table simulations are designed to run on HPC (Hellbender cluster). See `hellbender/HELLBENDER_GUIDE.md` for complete instructions.

### Quick Overview

1. Transfer files to Hellbender
2. Submit SLURM jobs for each table
3. Combine results after jobs complete
4. Download summary CSV files

## Scripts Overview

| Script | Table | Description |
|--------|-------|-------------|
| `hellbender/table1_type1_error.R` | Table 1 | Type I Error Control on Null Data |
| `hellbender/table2_power_analysis.R` | Table 2 | Statistical Power Analysis |
| `hellbender/table3_scalability.R` | Table 3 | Scalability Across Sample Size and Dimensionality |
| `hellbender/table4_parameter_sensitivity.R` | Table 4 | Parameter Sensitivity Analysis |
| `hellbender/table5_prefiltering.R` | Table 5 | Benefit of Conservative Pre-filtering |
| `hellbender/table6_feature_correlation.R` | Table 6 | Feature Correlation Robustness |
| `hellbender/skiome_validation.R` | - | SKIOME Real Data Validation |

## Notes

- All simulations use 200 permutations for MeLSI (conservative compared to 999 for traditional methods)
- Results are saved as CSV files for easy analysis
- Summary files are generated after combining individual simulation results
- All scripts use `set.seed(42)` for reproducibility

