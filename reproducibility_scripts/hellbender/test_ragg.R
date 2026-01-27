#!/usr/bin/env Rscript
# Test ragg package for PNG generation on headless systems

cat("=== Testing ragg package ===\n")

# Install ragg if needed
if (!requireNamespace("ragg", quietly = TRUE)) {
  cat("Installing ragg package...\n")
  install.packages("ragg", repos = "https://cloud.r-project.org", quiet = FALSE)
}

library(ragg)
cat("ragg package loaded successfully\n")

# Test creating a simple plot
cat("\nTesting agg_png() device...\n")
tryCatch({
  agg_png("test_ragg.png", width = 4, height = 4, units = "in", res = 300, background = "white")
  plot(1:10, main = "Test Plot")
  dev.off()
  
  if (file.exists("test_ragg.png")) {
    file_size <- file.info("test_ragg.png")$size
    cat("SUCCESS: PNG file created!\n")
    cat("File size:", file_size, "bytes\n")
    file.remove("test_ragg.png")
    cat("Test file cleaned up\n")
  } else {
    cat("FAILED: File was not created\n")
    stop("File creation failed")
  }
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  stop("Test failed: ", e$message)
})

cat("\n=== ragg test completed successfully ===\n")
