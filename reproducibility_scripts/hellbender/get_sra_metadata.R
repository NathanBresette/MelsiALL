#!/usr/bin/env Rscript
# ==============================================================================
# Get SRA Metadata for Group Labels
# ==============================================================================
# Downloads SRA metadata and extracts group labels (Atopic Dermatitis vs Psoriasis)
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("Getting SRA Metadata for SRP214545\n")
cat("==============================================================================\n\n")

# Try to download metadata using SRA Run Selector format
# URL: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545

cat("Attempting to download SRA metadata...\n")

# Method 1: Try direct download from SRA Run Selector
metadata_url <- "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545&go=go"

# Method 2: Try efetch API
efetch_url <- "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP214545"

cat("Trying efetch API...\n")
tryCatch({
  download.file(efetch_url, "sra_runinfo.csv", quiet = TRUE)
  if (file.exists("sra_runinfo.csv") && file.info("sra_runinfo.csv")$size > 0) {
    cat("✅ Downloaded SRA run info\n")
    runinfo <- read.csv("sra_runinfo.csv", stringsAsFactors = FALSE)
    cat("  Rows:", nrow(runinfo), "\n")
    cat("  Columns:", paste(colnames(runinfo), collapse = ", "), "\n")
    
    # Check if we have sample information
    if ("SampleName" %in% colnames(runinfo)) {
      cat("\nSample names found in runinfo\n")
      print(head(runinfo[, c("Run", "SampleName")], 10))
    }
  }
}, error = function(e) {
  cat("❌ Could not download via efetch:", e$message, "\n")
})

# Method 3: Query individual SRR samples for BioSample attributes
cat("\nQuerying individual samples for BioSample attributes...\n")

# Load the sample IDs we have
load("skiome_data_loaded.RData")
sample_ids <- rownames(counts)

# Query a few samples to see structure
test_samples <- head(sample_ids, 5)
cat("Testing with samples:", paste(test_samples, collapse = ", "), "\n")

# For each sample, try to get BioSample info
# We'll use the SRA web interface or API
cat("\nTo get group labels, we need to:\n")
cat("1. Download SRA Run Selector metadata manually from:\n")
cat("   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545\n")
cat("   (Click 'Run Selector' → 'Download' → 'Metadata')\n\n")
cat("2. Or check publication supplementary data:\n")
cat("   https://www.nature.com/articles/s41467-019-12253-y\n\n")
cat("3. Or query BioSample records for disease/condition attributes\n\n")

# Create a script to parse metadata once we have it
cat("Creating metadata parser script...\n")
parser_script <- '
# Once you have the SRA metadata file, run this to extract group labels
# Expected format: CSV with columns including Run, SampleName, or disease/condition info

if (!file.exists("sra_metadata.csv")) {
  stop("Please download SRA metadata first from SRA Run Selector")
}

metadata <- read.csv("sra_metadata.csv", stringsAsFactors = FALSE)
print(colnames(metadata))

# Look for disease/condition columns
disease_cols <- grep("disease|condition|diagnosis|group|phenotype", 
                     colnames(metadata), ignore.case = TRUE, value = TRUE)
cat("Disease-related columns:", paste(disease_cols, collapse = ", "), "\n")

# Match to our sample IDs
load("skiome_data_loaded.RData")
sample_ids <- rownames(counts)

# Create group labels
# This will need to be customized based on actual metadata structure
'
writeLines(parser_script, "parse_metadata.R")
cat("✅ Created parse_metadata.R\n")

cat("\n==============================================================================\n")
cat("NEXT STEPS:\n")
cat("==============================================================================\n")
cat("1. Download SRA metadata manually:\n")
cat("   - Go to: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545\n")
cat("   - Click 'Run Selector' button\n")
cat("   - Click 'Download' → 'Metadata'\n")
cat("   - Save as 'sra_metadata.csv' in this directory\n\n")
cat("2. Or check publication for sample metadata table\n\n")
cat("3. Then run: Rscript parse_metadata.R\n\n")
