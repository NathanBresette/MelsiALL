#!/usr/bin/env Rscript
# Test ggsave with default device

library(ggplot2)

cat("Testing ggsave with default device...\n")

# Create a simple plot
p <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x, y)) + geom_point()

# Try to save
test_file <- "test_ggsave.png"
cat("Attempting to save to:", test_file, "\n")

tryCatch({
  ggsave(test_file, p, width = 4, height = 4, dpi = 300, device = NULL, bg = "white")
  
  if (file.exists(test_file)) {
    file_size <- file.info(test_file)$size
    cat("SUCCESS: File created! Size:", file_size, "bytes\n")
    file.remove(test_file)
    cat("Test file cleaned up\n")
  } else {
    cat("FAILED: File was not created\n")
  }
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})

cat("Test complete\n")
