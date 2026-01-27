#!/usr/bin/env Rscript
# Install and test Cairo package

if (!requireNamespace("Cairo", quietly = TRUE)) {
  cat("Installing Cairo package...\n")
  install.packages("Cairo", repos = "https://cloud.r-project.org", quiet = FALSE)
  cat("Installation complete.\n")
} else {
  cat("Cairo package already installed.\n")
}

# Test CairoPNG
library(Cairo)
cat("Testing CairoPNG...\n")
CairoPNG("test_cairo_pkg.png", width = 100, height = 100)
plot(1:10)
dev.off()

if (file.exists("test_cairo_pkg.png")) {
  cat("SUCCESS: CairoPNG works! File created.\n")
  file.remove("test_cairo_pkg.png")
} else {
  cat("FAILED: No file created\n")
}
