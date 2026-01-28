#!/usr/bin/env Rscript
# Check Table 6 simulation file validity and identify missing/invalid files

# Check High correlation files (151-200)
cat("=== CHECKING HIGH CORRELATION FILES (151-200) ===\n")

invalid_files <- character()
missing_files <- character()
valid_count <- 0

for (i in 151:200) {
  file <- paste0("table6_sim_", i, ".csv")
  
  if (!file.exists(file)) {
    missing_files <- c(missing_files, as.character(i))
    cat("Missing:", file, "\n")
  } else {
    # Check if file is valid
    tryCatch({
      data <- read.csv(file, stringsAsFactors = FALSE)
      
      # Check if it has data rows and correct correlation level
      if (nrow(data) == 0) {
        invalid_files <- c(invalid_files, as.character(i))
        cat("Invalid (empty):", file, "\n")
      } else if (!"Correlation_Level" %in% colnames(data)) {
        invalid_files <- c(invalid_files, as.character(i))
        cat("Invalid (no Correlation_Level column):", file, "\n")
      } else if (!any(data$Correlation_Level == "High")) {
        invalid_files <- c(invalid_files, as.character(i))
        cat("Invalid (wrong correlation level):", file, "\n")
        cat("  Found levels:", unique(data$Correlation_Level), "\n")
      } else {
        valid_count <- valid_count + 1
      }
    }, error = function(e) {
      invalid_files <<- c(invalid_files, as.character(i))
      cat("Invalid (error reading):", file, "-", e$message, "\n")
    })
  }
}

cat("\n=== SUMMARY ===\n")
cat("Valid High correlation files:", valid_count, "/ 50\n")
cat("Missing files:", length(missing_files), "\n")
cat("Invalid files:", length(invalid_files), "\n")

if (length(missing_files) > 0) {
  cat("\nMissing file indices:", paste(missing_files, collapse = ", "), "\n")
}

if (length(invalid_files) > 0) {
  cat("\nInvalid file indices:", paste(invalid_files, collapse = ", "), "\n")
  cat("\nSLURM array for resubmission:", paste(invalid_files, collapse = ","), "\n")
}
