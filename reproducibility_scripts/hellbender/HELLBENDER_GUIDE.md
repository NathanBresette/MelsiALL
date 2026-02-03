# Hellbender HPC Guide for MeLSI Simulations

This guide explains how to run MeLSI validation simulations on Mizzou's Hellbender HPC cluster.

## Prerequisites

- SSH access to Hellbender: `ssh <your_username>@hellbender.rnet.missouri.edu`
- R module available (typically `r/4.4.0`)
- MeLSI package installed in your R library

## Quick Start

### 1. Transfer Files to Hellbender

From your local machine:

```bash
cd <local_path>/reproducibility_scripts/hellbender
scp -r *.R *.sh <your_username>@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
```

### 2. SSH into Hellbender

```bash
ssh <your_username>@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
```

### 3. Load R Module and Set Environment

```bash
module load r/4.4.0
export R_LIBS_USER=~/R
```

### 4. Submit Jobs

Each table has its own job script. Submit them individually:

```bash
# Table 1: Type I Error Analysis
sbatch table1_type1_error_job.sh

# Table 2: Power Analysis
sbatch table2_power_analysis_job.sh

# Table 3: Scalability Analysis
sbatch table3_scalability_job.sh

# Table 4: Parameter Sensitivity
sbatch table4_parameter_sensitivity_job.sh

# Table 5: Pre-filtering Analysis
sbatch table5_prefiltering_job.sh

# Table 6: Feature Correlation Analysis
sbatch table6_feature_correlation_job.sh

# SKIOME Validation
sbatch skiome_validation_job.sh
```

### 5. Monitor Jobs

```bash
# Check job status
squeue -u $USER

# View job output (replace JOBID with actual job ID)
tail -f table1_type1_error_JOBID.out
tail -f table1_type1_error_JOBID.err
```

### 6. Combine Results

After jobs complete, combine parallel results:

```bash
module load r/4.4.0
export R_LIBS_USER=~/R

# Combine results for each table
Rscript combine_table1_results.R
Rscript combine_table2_results.R
Rscript combine_table3_results.R
Rscript combine_table4_results.R
Rscript combine_table5_results.R
Rscript combine_table6_results.R
```

### 7. Download Results

From your local machine:

```bash
scp <your_username>@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_summary.csv ./
```

## R Scripts Overview

### Simulation Scripts
- `table1_type1_error.R` - Type I error control validation
- `table2_power_analysis.R` - Statistical power analysis
- `table3_scalability.R` - Scalability across sample sizes and dimensions
- `table4_parameter_sensitivity.R` - Hyperparameter sensitivity
- `table5_prefiltering.R` - Pre-filtering benefit analysis
- `table6_feature_correlation.R` - Feature correlation robustness
- `skiome_validation.R` - SKIOME real data validation

### Combine Scripts
- `combine_table1_results.R` - Combines Table 1 simulation results
- `combine_table2_results.R` - Combines Table 2 simulation results
- `combine_table3_results.R` - Combines Table 3 simulation results
- `combine_table4_results.R` - Combines Table 4 simulation results
- `combine_table5_results.R` - Combines Table 5 simulation results
- `combine_table6_results.R` - Combines Table 6 simulation results

## Monitoring Jobs

### Check Job Status
```bash
squeue -u $USER
```

### View Recent Output
```bash
tail -50 table1_type1_error_JOBID.out
tail -50 table1_type1_error_JOBID.err
```

### Follow Output in Real-Time
```bash
tail -f table1_type1_error_JOBID.out table1_type1_error_JOBID.err
```

## Troubleshooting

### R Module Not Found
```bash
module avail r
module load r/4.4.0
```

### R Packages Not Found
```bash
module load r/4.4.0
export R_LIBS_USER=~/R
Rscript -e "library(MeLSI)"
```

### Job Fails Immediately
- Check error file: `cat table1_type1_error_JOBID.err`
- Verify R module loads: `module load r/4.4.0 && which Rscript`

## Notes

- All simulations use 200 permutations for MeLSI (conservative compared to 999 for traditional methods)
- Results are saved as CSV files for easy analysis
- Summary files are generated after combining individual simulation results
