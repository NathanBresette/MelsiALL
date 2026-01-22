# Hellbender HPC Setup Guide for MeLSI Simulations

## Overview
This guide explains how to run the Type I error and power analysis simulations on Mizzou's Hellbender HPC cluster using parallel job arrays.

## Prerequisites
- SSH access to Hellbender: `ssh nbhtd@hellbender.rnet.missouri.edu`
- R and required packages installed (or available via modules)
- MeLSI package installed

## Quick Start (Easiest Method)

### From Your Local Machine:

```bash
# Navigate to the hellbender directory
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender

# Transfer all files to Hellbender
bash transfer_to_hellbender.sh
```

### On Hellbender:

```bash
# SSH in
ssh nbhtd@hellbender.rnet.missouri.edu

# Navigate to directory
cd ~/melsi_simulations/hellbender

# Run quick start (checks setup and submits jobs)
bash quick_start.sh
```

That's it! The jobs are now running. Monitor with `bash monitor_jobs.sh`

---

## Detailed Manual Setup

### Step 1: Transfer Files to Hellbender

From your local machine:

```bash
# Navigate to the hellbender directory
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender

# Use the transfer script (recommended)
bash transfer_to_hellbender.sh

# OR manually transfer
scp -r * nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
```

## Step 2: SSH into Hellbender and Set Up

```bash
ssh nbhtd@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
```

## Step 3: Check Available R Modules

```bash
# Check available R versions
module avail R

# Load R (adjust version as needed)
module load R/4.3.0  # or whatever version is available

# Verify R is loaded
which R
R --version
```

## Step 4: Install Required R Packages (if needed)

```bash
# Start R
R

# Install packages (if not already installed)
install.packages(c("vegan", "microbiome", "GUniFrac", "ape"))

# Install MeLSI package (adjust path as needed)
# install.packages("~/melsi_simulations/github", repos=NULL, type="source")

# Exit R
q()
```

## Step 5: Run Setup Check

```bash
# Check environment and dependencies
bash setup_and_check.sh
```

This will verify:
- R is available
- Required packages are installed
- All scripts are present
- SLURM is available

## Step 6: Submit Jobs

### Option A: Use Quick Start (Recommended)
```bash
bash quick_start.sh
```

### Option B: Submit Manually

### Submit Type I Error Analysis (300 simulations)
```bash
sbatch table1_type1_error_job.sh
```

### Submit Power Analysis (270 simulations)
```bash
sbatch table2_power_analysis_job.sh
```

## Step 7: Monitor Jobs

### Easy Monitoring:
```bash
# Use the monitoring script (shows progress and errors)
bash monitor_jobs.sh
```

### Manual Monitoring:
```bash
# Check job status
squeue -u nbhtd

# Check specific job details
squeue -j <JOB_ID>

# Cancel a job if needed
scancel <JOB_ID>
```

## Step 8: Check Job Output

```bash
# List output files
ls -lh table1_type1_*.out table1_type1_*.err
ls -lh table2_power_*.out table2_power_*.err

# Check if simulations are completing
ls -1 table1_sim_*.csv | wc -l  # Should approach 300
ls -1 table2_sim_*.csv | wc -l  # Should approach 270
```

## Step 9: Combine Results (After All Jobs Complete)

### Combine Type I Error Results
```bash
module load R
Rscript combine_table1_results.R
```

### Combine Power Analysis Results
```bash
Rscript combine_table2_results.R
```

## Step 10: Download Results

From your local machine:

```bash
# Download combined results
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/table1_type1_error_results.csv ./
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/table1_type1_error_summary.csv ./
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/table2_power_analysis_results.csv ./
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/table2_power_analysis_summary.csv ./
```

Or download all results at once:
```bash
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_summary.csv ./
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_results.csv ./
```

## Troubleshooting

### Jobs Not Starting
- Check queue status: `squeue`
- Check job requirements: `scontrol show job <JOB_ID>`
- Verify R module is loaded in the script

### Missing R Packages
- Install packages in your home directory or request system-wide installation
- Check R library paths: `Rscript -e ".libPaths()"`

### Simulation Files Not Generated
- Check error files: `cat table1_type1_*_*.err | head -20`
- Verify R script paths are correct
- Check that MeLSI package is accessible

### Need More Resources
- Adjust `--time`, `--mem`, or `--cpus-per-task` in the job scripts
- Contact HPC support if you need more resources

## Expected Runtime

- **Type I Error**: 300 simulations × ~2-5 min each = ~10-25 hours (but parallelized to ~15-40 min with 50+ cores)
- **Power Analysis**: 270 simulations × ~2-5 min each = ~9-22 hours (but parallelized to ~15-35 min with 50+ cores)

## Documentation

- Hellbender Documentation: https://docs.itrss.umsystem.edu/pub/hpc/hellbender
- Hellbender DOI: https://doi.org/10.32469/10355/97710
