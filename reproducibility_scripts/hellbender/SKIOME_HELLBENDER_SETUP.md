# SKIOME Validation on Hellbender - Setup Instructions

## Files Created
- `skiome_validation.R` - Main validation script (multi-group MeLSI analysis)
- `skiome_validation_job.sh` - SLURM job script
- `skiome_data_loaded.RData` - Pre-loaded data with real group labels
- `SraRunTable.csv` - SRA metadata (for reference)

## Data Summary
- **Samples:** 511
- **Taxa:** 1,856
- **Groups:** 3 (Atopic_Dermatitis: 109, Healthy: 250, Psoriasis: 152)
- **Analysis Type:** Multi-group (omnibus + pairwise comparisons)

## Transfer Files to Hellbender

### Option 1: Manual Transfer (Recommended)
```bash
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender

# Transfer files one by one (you'll be prompted for password)
scp skiome_validation.R skiome_validation_job.sh nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
scp skiome_data_loaded.RData nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
scp SraRunTable.csv nbhtd@hellbender.rnet.missouri.edu:~/melsi_simulations/hellbender/
```

### Option 2: Use Transfer Script
```bash
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender
bash transfer_skiome_to_hellbender.sh
# (You'll need to enter password when prompted)
```

## Run on Hellbender

### Step 1: SSH into Hellbender
```bash
ssh nbhtd@hellbender.rnet.missouri.edu
```

### Step 2: Navigate to Directory
```bash
cd ~/melsi_simulations/hellbender
```

### Step 3: Verify Files
```bash
ls -lh skiome*
ls -lh skiome_data_loaded.RData SraRunTable.csv
```

### Step 4: Submit Job
```bash
sbatch skiome_validation_job.sh
```

### Step 5: Monitor Job
```bash
# Check job status
squeue -u $USER

# Check output (once job starts)
tail -f skiome_validation_*.out

# Check for errors
tail -f skiome_validation_*.err
```

## Expected Output

The script will:
1. Load the SKIOME data (511 samples, 3 groups)
2. Run MeLSI multi-group analysis (omnibus + pairwise)
3. Run traditional methods (Euclidean, Bray-Curtis, Jaccard PERMANOVA)
4. Save results to `skiome_validation_results.csv`

## Results File

`skiome_validation_results.csv` will contain:
- Method (MeLSI, Euclidean, Bray-Curtis, Jaccard)
- F-statistic
- P-value
- Time (seconds)
- Significant (TRUE/FALSE)

## Multi-Group Analysis

MeLSI will automatically:
- Run **omnibus test** (all 3 groups together)
- Run **pairwise comparisons** (Atopic_Dermatitis vs Healthy, Atopic_Dermatitis vs Psoriasis, Healthy vs Psoriasis)
- Apply **multiple testing correction** (Benjamini-Hochberg FDR)

## Estimated Runtime

- **MeLSI:** ~30-60 minutes (200 permutations, 3 groups)
- **Traditional methods:** ~5-10 minutes total
- **Total:** ~40-70 minutes

## Troubleshooting

If job fails:
1. Check error log: `cat skiome_validation_*.err`
2. Verify R packages are installed: `module load r/4.4.0 && Rscript -e "library(MeLSI)"`
3. Check data file exists: `ls -lh skiome_data_loaded.RData`
4. Verify group labels: `Rscript -e "load('skiome_data_loaded.RData'); print(table(metadata\$Group))"`
