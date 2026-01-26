#!/usr/bin/env Rscript
# ==============================================================================
# Load SKIOME Data (PRJNA554499) from MGnify
# ==============================================================================
# Loads the taxonomy abundance table from MGnify and prepares it for MeLSI
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("Loading SKIOME Data: PRJNA554499 (SRP214545)\n")
cat("Source: MGnify processed data\n")
cat("==============================================================================\n\n")

# Data file location
data_file <- "../extdata/SRP214545_taxonomy_abundances_SSU_v5.0.tsv"

if (!file.exists(data_file)) {
  # Try alternative locations
  alt_locations <- c(
    "extdata/SRP214545_taxonomy_abundances_SSU_v5.0.tsv",
    "../../extdata/SRP214545_taxonomy_abundances_SSU_v5.0.tsv",
    "SRP214545_taxonomy_abundances_SSU_v5.0.tsv"
  )
  
  for (alt in alt_locations) {
    if (file.exists(alt)) {
      data_file <- alt
      break
    }
  }
}

if (!file.exists(data_file)) {
  cat("❌ Error: Could not find data file\n")
  cat("Expected location:", data_file, "\n")
  cat("\nPlease download from MGnify:\n")
  cat("https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00005299/pipelines/5.0/file/SRP214545_taxonomy_abundances_SSU_v5.0.tsv\n")
  quit(status = 1)
}

cat("Loading abundance table from:", data_file, "\n")

# Read the abundance table
# First column is taxonomy (#SampleID), rest are samples (SRR numbers)
# The header starts with #SampleID, then has SRR numbers
counts_raw <- readLines(data_file)

# Remove the # from the first line to make it a proper header
if (startsWith(counts_raw[1], "#")) {
  counts_raw[1] <- sub("^#", "", counts_raw[1])
}

# Write to temp file and read properly
temp_file <- tempfile()
writeLines(counts_raw, temp_file)

counts <- read.delim(temp_file, sep = "\t", header = TRUE, 
                     check.names = FALSE, stringsAsFactors = FALSE,
                     fill = TRUE, quote = "", row.names = 1)

unlink(temp_file)

cat("Raw data dimensions:", nrow(counts), "taxa x", ncol(counts), "samples\n")

# Transpose to samples x taxa (standard format)
counts <- t(counts)

cat("Transposed dimensions:", nrow(counts), "samples x", ncol(counts), "taxa\n")

# Remove taxa with zero counts across all samples
taxa_sums <- colSums(counts)
non_zero_taxa <- taxa_sums > 0
counts <- counts[, non_zero_taxa, drop = FALSE]

cat("After removing zero taxa:", nrow(counts), "samples x", ncol(counts), "taxa\n\n")

# Now we need metadata - check SRA for sample information
# For now, create a placeholder metadata structure
# We'll need to get actual group labels from SRA or publication

cat("Creating metadata structure...\n")
cat("NOTE: We need to get group labels (Atopic Dermatitis vs Psoriasis) from SRA metadata\n")
cat("For now, creating placeholder metadata...\n\n")

# Create basic metadata
metadata <- data.frame(
  SampleID = rownames(counts),
  stringsAsFactors = FALSE
)

# Try to infer groups from sample names or load from file if available
# For PRJNA554499, we need to identify which samples are AD vs Psoriasis
# This will require downloading SRA sample metadata

cat("Sample IDs (first 10):\n")
print(head(rownames(counts), 10))

cat("\nTo complete metadata, we need to:\n")
cat("1. Download SRA sample metadata for SRP214545\n")
cat("2. Match sample IDs (SRR numbers) to group labels\n")
cat("3. Identify Atopic Dermatitis vs Psoriasis samples\n\n")

# Save what we have so far
save(counts, metadata, file = "skiome_data_loaded.RData")
cat("✅ Data loaded and saved to: skiome_data_loaded.RData\n")
cat("  Samples:", nrow(counts), "\n")
cat("  Taxa:", ncol(counts), "\n")
cat("  Metadata columns:", ncol(metadata), "\n\n")

cat("Next step: Download SRA metadata to get group labels\n")
