#!/usr/bin/env Rscript
# Check all Table 6 correlation levels

cat("=== CHECKING ALL CORRELATION LEVELS ===\n\n")

for (level in c("None", "Low", "Moderate", "High")) {
  files <- list.files(pattern = "table6_sim_.*\\.csv")
  valid <- 0
  
  for (f in files) {
    data <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(data) && nrow(data) > 0 && "Correlation_Level" %in% colnames(data) && any(data$Correlation_Level == level)) {
      valid <- valid + 1
    }
  }
  
  cat(sprintf("%s: %d / 50\n", level, valid))
}

cat("\n=== SUMMARY ===\n")
cat("All levels should have 50/50 simulations\n")
