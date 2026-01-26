# How to Get Real Group Labels for SKIOME Data

## ⚠️ IMPORTANT: We Need REAL Labels, Not Inferred Ones

The automated download isn't working, so you need to manually download the SRA metadata file which contains the **real** disease/condition labels for each sample.

## Step-by-Step Instructions

### Step 1: Open SRA Run Selector
1. Open this URL in your browser:
   **https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545**

2. The page should automatically show the "Run Selector" interface with a table of all 511 runs

### Step 2: Download Metadata
1. Look for the **"Download"** button (usually in the top right area of the Run Selector table)
2. Click the **"Download"** button
3. A dropdown menu should appear
4. Select **"Metadata"** from the dropdown
5. This will download a CSV file (usually named something like `SraRunInfo.csv` or `metadata.csv`)

### Step 3: Save the File
1. Save the downloaded file as **`sra_metadata.csv`** in this directory:
   `/Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender/`

### Step 4: Parse the Metadata
Once the file is saved, run:
```bash
cd /Users/nathanbresette/Desktop/MeLSI/reproducibility_scripts/hellbender
Rscript parse_metadata.R
```

This will:
- Read the metadata file
- Extract group labels (Atopic Dermatitis vs Psoriasis)
- Match them to your sample IDs (SRR numbers)
- Update the `skiome_data_loaded.RData` file with real group labels

## What the Metadata File Should Contain

The metadata file should have columns like:
- `Run` (SRR numbers, e.g., SRR9678960)
- `SampleName` or similar
- `BioSample`
- Disease/condition information (may be in `SampleName`, `sample_title`, or a dedicated column)

## Alternative: Check Publication Supplementary Data

If the SRA metadata doesn't have clear group labels, check the publication:
- **DOI:** 10.1038/s41467-019-12253-y
- **Link:** https://www.nature.com/articles/s41467-019-12253-y
- Look for "Supplementary Data" or "Data Availability" section
- May contain a sample metadata table with group assignments

## After Getting Metadata

Once you have the real group labels:
1. The `parse_metadata.R` script will automatically match them to your samples
2. The updated data will be saved to `skiome_data_loaded.RData`
3. You can then proceed with MeLSI validation using the real group labels
