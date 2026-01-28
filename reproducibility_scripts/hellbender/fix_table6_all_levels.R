#!/usr/bin/env Rscript
# Identify missing simulations for each correlation level and delete invalid files

cat("=== IDENTIFYING MISSING SIMULATIONS FOR EACH LEVEL ===\n\n")

correlation_levels <- c(0.0, 0.3, 0.6, 0.8)
correlation_names <- c("None", "Low", "Moderate", "High")
n_simulations_per_condition <- 50

# Expected file ranges for each correlation level
expected_ranges <- list(
  "None" = 1:50,
  "Low" = 51:100,
  "Moderate" = 101:150,
  "High" = 151:200
)

for (i in seq_along(correlation_names)) {
  level <- correlation_names[i]
  corr_value <- correlation_levels[i]
  expected_files <- expected_ranges[[level]]
  
  cat(sprintf("\n=== %s correlation (r=%.1f) ===\n", level, corr_value))
  
  valid_files <- c()
  invalid_files <- c()
  missing_files <- c()
  
  for (file_idx in expected_files) {
    file <- paste0("table6_sim_", file_idx, ".csv")
    
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file_idx)
    } else {
      data <- tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
      
      if (is.null(data) || nrow(data) == 0 || !"Correlation_Level" %in% colnames(data)) {
        invalid_files <- c(invalid_files, file_idx)
      } else {
        file_level <- unique(data$Correlation_Level)
        if (any(file_level == level)) {
          valid_files <- c(valid_files, file_idx)
        } else {
          invalid_files <- c(invalid_files, file_idx)
          cat(sprintf("  File %d has wrong level: %s (should be %s)\n", file_idx, paste(file_level, collapse=", "), level))
        }
      }
    }
  }
  
  cat(sprintf("Valid: %d / 50\n", length(valid_files)))
  cat(sprintf("Invalid: %d\n", length(invalid_files)))
  cat(sprintf("Missing: %d\n", length(missing_files)))
  
  if (length(invalid_files) > 0) {
    cat(sprintf("Invalid file indices: %s\n", paste(invalid_files, collapse = ", ")))
  }
  
  if (length(missing_files) > 0) {
    cat(sprintf("Missing file indices: %s\n", paste(missing_files, collapse = ", ")))
  }
  
  # Files that need to be resubmitted (invalid + missing)
  need_resubmit <- sort(c(invalid_files, missing_files))
  if (length(need_resubmit) > 0) {
    cat(sprintf("Need to resubmit: %s\n", paste(need_resubmit, collapse = ", ")))
  }
}
