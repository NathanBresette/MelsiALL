# SKIOME Validation Job Status

## ✅ Files Successfully Transferred

All files have been transferred to Hellbender:
- `skiome_validation.R` - Main validation script
- `skiome_validation_job.sh` - SLURM job script  
- `skiome_data_loaded.RData` - Data with real group labels (511 samples, 3 groups)

## ⚠️ Job Submission Issue

**Status:** Files transferred, but job submission blocked by partition access

**Error:** `User's group not permitted to use this partition`

This is the same partition access issue we've encountered before. The files are ready on Hellbender, but SLURM is defaulting to the "requeue" partition which your group doesn't have access to.

## Solutions

### Option 1: Contact Hellbender Support
Ask them to grant your group access to a partition, or find out which partition your group can use.

### Option 2: Try Manual Submission Later
When partition access is available, SSH and submit:
```bash
ssh nbhtd@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
sbatch skiome_validation_job.sh
```

### Option 3: Run Directly (if partition access isn't available soon)
You could run the script directly without SLURM (though it may be slower):
```bash
ssh nbhtd@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
module load r/4.4.0
export R_LIBS_USER=~/R
Rscript skiome_validation.R
```

## What Will Run

Once the job runs, it will:
1. Load SKIOME data (511 samples, 3 groups: Atopic_Dermatitis, Healthy, Psoriasis)
2. Run MeLSI multi-group analysis (omnibus + pairwise comparisons)
3. Compare to traditional methods (Euclidean, Bray-Curtis, Jaccard PERMANOVA)
4. Save results to `skiome_validation_results.csv`

**Expected runtime:** ~40-70 minutes

## Files Location on Hellbender

All files are in: `~/melsi_simulations/hellbender/`

You can verify with:
```bash
ssh nbhtd@hellbender.rnet.missouri.edu
cd ~/melsi_simulations/hellbender
ls -lh skiome*
```
