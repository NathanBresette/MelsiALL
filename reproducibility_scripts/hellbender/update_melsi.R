#!/usr/bin/env Rscript
# Update MeLSI package with new melsi_robust.R

# Find MeLSI package location
melsi_path <- find.package("MeLSI")
cat("MeLSI package found at:", melsi_path, "\n")

# Copy updated file
source_file <- "melsi_robust.R"
target_file <- file.path(melsi_path, "R", "melsi_robust.R")

if (file.exists(source_file)) {
  file.copy(source_file, target_file, overwrite = TRUE)
  cat("Successfully updated:", target_file, "\n")
  
  # Verify update
  if (file.exists(target_file)) {
    cat("Verification: File exists at target location\n")
  } else {
    stop("ERROR: File was not copied successfully")
  }
} else {
  stop("ERROR: Source file 'melsi_robust.R' not found in current directory")
}

cat("\nMeLSI package updated successfully!\n")
