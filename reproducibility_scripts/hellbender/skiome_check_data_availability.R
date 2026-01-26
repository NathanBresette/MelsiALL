#!/usr/bin/env Rscript
# ==============================================================================
# Check SKIOME Data Availability
# ==============================================================================
# Checks SRA, MGnify, and publication for pre-processed data availability
# for PRJNA554499 (SRP214545)
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("Checking Data Availability for PRJNA554499 (SKIOME)\n")
cat("==============================================================================\n\n")

# Dataset info
bioproject <- "PRJNA554499"
sra_study <- "SRP214545"
mgnify_id <- "MGYS00005299"
doi <- "10.1038/s41467-019-12253-y"
publication <- "Nature Communications, 2019"

cat("Dataset Information:\n")
cat("  BioProject:", bioproject, "\n")
cat("  SRA Study:", sra_study, "\n")
cat("  MGnify ID:", mgnify_id, "\n")
cat("  Publication:", publication, "\n")
cat("  DOI:", doi, "\n")
cat("  Expected samples: 531\n")
cat("  Disease: Atopic Dermatitis and Psoriasis\n\n")

cat("Data Sources to Check:\n")
cat("  1. SRA supplementary files\n")
cat("  2. MGnify processed data (MGYS00005299)\n")
cat("  3. Publication supplementary materials\n")
cat("  4. Raw FASTQ files in SRA\n\n")

cat("==============================================================================\n")
cat("RECOMMENDED APPROACH:\n")
cat("==============================================================================\n")
cat("1. Check MGnify first: https://www.ebi.ac.uk/metagenomics/projects/MGYS00005299\n")
cat("   - Look for 'Download' or 'Processed data' section\n")
cat("   - May have OTU/ASV tables available\n\n")

cat("2. Check publication supplementary:\n")
cat("   - Visit: https://www.nature.com/articles/s41467-019-12253-y\n")
cat("   - Look for 'Supplementary Data' or 'Data Availability' section\n")
cat("   - May have processed OTU tables or links to data repositories\n\n")

cat("3. If no pre-processed data, download from SRA:\n")
cat("   - SRA Run Selector: https://www.ncbi.nlm.nih.gov/sra?term=PRJNA554499\n")
cat("   - Download raw FASTQ files\n")
cat("   - Process with DADA2/QIIME2 on Hellbender\n\n")

cat("4. Alternative: Use a subset of samples (300-400) for computational efficiency\n\n")

cat("==============================================================================\n")
cat("NEXT STEPS:\n")
cat("==============================================================================\n")
cat("Run this script to check availability, then proceed with download/processing\n")
cat("based on what's available.\n\n")
