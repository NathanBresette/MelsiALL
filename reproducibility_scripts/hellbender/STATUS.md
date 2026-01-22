# Hellbender Simulation Status

## ✅ Setup Complete

All required packages have been installed:
- ✅ R 4.4.0 loaded
- ✅ vegan, GUniFrac, ape installed
- ✅ microbiome, phyloseq installed  
- ✅ MeLSI package installed
- ✅ All dependencies resolved

## 📊 Jobs Running ✅

**Type I Error Analysis**: Job array submitted (300 simulations) - **Job ID: 12315381**
**Power Analysis**: Job array submitted (270 simulations) - **Job ID: 12315382**

**Status**: Simulations are running successfully! Result files are being generated.

## 🔍 Monitor Progress

From your local machine:
```bash
ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && bash monitor_jobs.sh"
```

Or check directly:
```bash
ssh nbhtd@hellbender.rnet.missouri.edu "cd ~/melsi_simulations/hellbender && squeue -u \$USER"
```

## 📥 After Completion

Once all jobs complete (check with `monitor_jobs.sh`):

1. **Combine results on Hellbender:**
   ```bash
   ssh nbhtd@hellbender.rnet.missouri.edu
   cd ~/melsi_simulations/hellbender
   module load r/4.4.0
   export R_LIBS_USER=~/R
   Rscript combine_table1_results.R
   Rscript combine_table2_results.R
   ```

2. **Download results to local machine:**
   ```bash
   cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender
   scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_summary.csv ./
   scp nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/*_results.csv ./
   ```

## ⏱️ Expected Runtime

- **Sequential**: ~20+ hours
- **Parallel (Hellbender)**: ~15-40 minutes (depending on cluster load)

**Current Status**: Jobs are running. Each simulation takes 1-3 minutes due to computational intensity (200 permutations × 30 bootstrap samples per simulation).

## 📝 Notes

- All packages installed in `~/R` directory
- Job scripts set `R_LIBS_USER=~/R` to find packages
- MeLSI DESCRIPTION updated to require R >= 4.4.0 (was 4.5.0)
