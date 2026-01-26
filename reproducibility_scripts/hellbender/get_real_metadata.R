#!/usr/bin/env Rscript
# ==============================================================================
# Get Real Group Labels from SRA/BioSample
# ==============================================================================
# This script provides instructions and attempts to download real metadata
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("Getting REAL Group Labels for SRP214545\n")
cat("==============================================================================\n\n")

cat("To get the REAL group labels, you need to download the SRA metadata file.\n")
cat("This file contains the actual disease/condition information for each sample.\n\n")

cat("METHOD 1: Manual Download from SRA Run Selector (RECOMMENDED)\n")
cat("----------------------------------------------------------------------\n")
cat("1. Open this URL in your browser:\n")
cat("   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP214545\n\n")
cat("2. The page should show a 'Run Selector' interface\n")
cat("3. Look for a 'Download' button (usually top right)\n")
cat("4. Click 'Download' → Select 'Metadata' from dropdown\n")
cat("5. Save the file as 'sra_metadata.csv' in this directory\n\n")

cat("METHOD 2: Try Direct Download URL\n")
cat("----------------------------------------------------------------------\n")
cat("Attempting to construct download URL...\n")

# Try common SRA metadata download patterns
base_url <- "https://www.ncbi.nlm.nih.gov/Traces/study/"
study_acc <- "SRP214545"

# Pattern 1: Direct metadata download
url1 <- paste0(base_url, "?acc=", study_acc, "&go=go")
cat("URL 1:", url1, "\n")

# Pattern 2: Run info download
url2 <- paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", study_acc)
cat("URL 2:", url2, "\n")

cat("\nTrying URL 2...\n")
tryCatch({
  download.file(url2, "sra_runinfo_temp.csv", quiet = TRUE, method = "libcurl")
  if (file.exists("sra_runinfo_temp.csv") && file.info("sra_runinfo_temp.csv")$size > 100) {
    cat("✅ Successfully downloaded run info!\n")
    runinfo <- read.csv("sra_runinfo_temp.csv", stringsAsFactors = FALSE)
    cat("  Rows:", nrow(runinfo), "\n")
    cat("  Columns:", paste(colnames(runinfo), collapse = ", "), "\n")
    
    # Check for disease/group columns
    disease_cols <- grep("disease|condition|diagnosis|group|phenotype|sample_name|sample_title", 
                         colnames(runinfo), ignore.case = TRUE, value = TRUE)
    if (length(disease_cols) > 0) {
      cat("\n✅ Found disease/group columns:", paste(disease_cols, collapse = ", "), "\n")
      cat("\nFirst few rows:\n")
      print(head(runinfo[, c("Run", disease_cols[1])], 10))
      
      # Rename to final file
      file.rename("sra_runinfo_temp.csv", "sra_metadata.csv")
      cat("\n✅ Saved as sra_metadata.csv\n")
      cat("Now run: Rscript parse_metadata.R\n")
    } else {
      cat("\n⚠️  Downloaded but no disease columns found.\n")
      cat("Columns available:", paste(colnames(runinfo), collapse = ", "), "\n")
      cat("First row:\n")
      print(runinfo[1, ])
    }
  } else {
    cat("❌ Download failed or file too small\n")
  }
}, error = function(e) {
  cat("❌ Error:", e$message, "\n")
  cat("\nPlease use METHOD 1 (manual download) instead.\n")
})

cat("\n==============================================================================\n")
cat("Once you have sra_metadata.csv, run: Rscript parse_metadata.R\n")
cat("==============================================================================\n\n")
