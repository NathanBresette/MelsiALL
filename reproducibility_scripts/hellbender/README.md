# Hellbender Simulation Setup - Quick Start

## 🚀 Get Running in 3 Steps

### Step 1: Transfer Files (from your local machine)

```bash
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender
scp -r * nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
```

You'll be prompted for your password. Enter it when asked.

### Step 2: SSH into Hellbender and Run Setup

```bash
ssh nbhtd@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
bash quick_start.sh
```

This will:
- ✅ Check your environment (R, packages, SLURM)
- ✅ Submit both job arrays automatically
- ✅ Show you the job IDs

### Step 3: Monitor Progress

```bash
# Check job status
squeue -u $USER

# Or use the monitoring script
bash monitor_jobs.sh
```

## 📊 After Jobs Complete

```bash
# Combine results
Rscript combine_table1_results.R
Rscript combine_table2_results.R

# Download to local machine (from your local terminal)
scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_summary.csv ./
```

## 📁 Files in This Directory

- **R Scripts**: `table1_type1_error.R`, `table2_power_analysis.R` - Main simulation scripts
- **Job Scripts**: `table1_type1_error_job.sh`, `table2_power_analysis_job.sh` - SLURM job submissions
- **Combine Scripts**: `combine_table1_results.R`, `combine_table2_results.R` - Combine parallel results
- **Helper Scripts**:
  - `quick_start.sh` - Automated setup and job submission
  - `setup_and_check.sh` - Environment verification
  - `monitor_jobs.sh` - Progress monitoring
  - `transfer_to_hellbender.sh` - File transfer (run from local machine)

## ⚙️ What Gets Run

- **Type I Error Analysis**: 300 simulations (2 dataset types × 3 sample sizes × 50 sims each)
- **Power Analysis**: 270 simulations (3 effect sizes × 3 sample sizes × 30 sims each)

**Expected Runtime**: ~15-40 minutes (parallelized) vs 20+ hours sequential

## 📖 Full Documentation

See `HELLBENDER_SETUP.md` for detailed instructions and troubleshooting.
