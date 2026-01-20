#!/usr/bin/env Rscript
# Quick test to verify VIP directionality plotting works

suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(dplyr)
})

# Source the MeLSI function directly from source to ensure we have the latest version
source("github/R/melsi_robust.R")

set.seed(42)

cat("Testing VIP directionality plotting...\n\n")

# Load a small subset of Atlas1006 for quick test
data(atlas1006)
counts <- t(abundances(atlas1006))
metadata <- meta(atlas1006)
sex <- metadata$sex

# Remove samples with missing labels
valid_idx <- !is.na(sex)
counts <- counts[valid_idx, , drop = FALSE]
sex <- sex[valid_idx]

# Use only first 50 samples and 100 features for speed
counts <- counts[1:50, 1:100]
sex <- sex[1:50]

cat("Test data: ", nrow(counts), "samples, ", ncol(counts), "features\n")
cat("Groups:", paste(unique(sex), collapse = ", "), "\n\n")

# CLR transformation
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)

# Run MeLSI with fewer permutations for quick test
cat("Running MeLSI (n_perms = 10 for quick test)...\n")
melsi_res <- melsi(
  counts_clr,
  sex,
  n_perms = 10,
  B = 10,
  m_frac = 0.8,
  show_progress = TRUE,
  plot_vip = FALSE
)

cat("\n=== Checking Results ===\n")
cat("F-statistic:", round(melsi_res$F_observed, 3), "\n")
cat("P-value:", round(melsi_res$p_value, 3), "\n")

# Check directionality
has_directionality <- !is.null(melsi_res$directionality) && length(melsi_res$directionality) > 0
cat("Directionality available:", has_directionality, "\n")

if (has_directionality) {
  cat("Directionality length:", length(melsi_res$directionality), "\n")
  cat("Unique directions:", paste(unique(melsi_res$directionality), collapse = ", "), "\n")
  cat("First 5 directionality values:\n")
  print(head(melsi_res$directionality, 5))
} else {
  cat("ERROR: Directionality not found in results!\n")
  cat("Available result names:", paste(names(melsi_res), collapse = ", "), "\n")
}

# Test VIP plot creation
cat("\n=== Testing VIP Plot Creation ===\n")
vip_df <- data.frame(
  Feature = names(melsi_res$feature_weights),
  Weight = as.numeric(melsi_res$feature_weights)
) %>%
  arrange(desc(Weight)) %>%
  slice_head(n = 10) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

cat("Top 10 features selected\n")

if (has_directionality) {
  # Try to match directionality
  vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$Feature)]
  cat("Directionality matched for", sum(!is.na(vip_df$Directionality)), "out of", nrow(vip_df), "features\n")
  
  if (any(is.na(vip_df$Directionality))) {
    cat("WARNING: Some features have NA directionality\n")
    cat("Features with NA:", paste(vip_df$Feature[is.na(vip_df$Directionality)], collapse = ", "), "\n")
  }
  
  vip_df$Directionality <- ifelse(is.na(vip_df$Directionality), "Unknown", vip_df$Directionality)
  unique_dirs <- unique(vip_df$Directionality)
  cat("Unique directions in plot:", paste(unique_dirs, collapse = ", "), "\n")
  
  # Try to create plot
  tryCatch({
    if (length(unique_dirs) >= 2) {
      color_map <- setNames(c("#E63946", "#457B9D"), unique_dirs[1:2])
      if (any(vip_df$Directionality == "Unknown")) {
        color_map["Unknown"] <- "#95A5A6"
      }
    } else {
      color_map <- setNames(c("#1f77b4", "#95A5A6"), unique_dirs)
    }
    
    vip_plot <- ggplot(vip_df, aes(x = Weight, y = Feature, fill = Directionality)) +
      geom_col() +
      scale_fill_manual(values = color_map, name = "Higher in", drop = FALSE) +
      labs(title = "Test VIP Plot", x = "Feature Weight", y = NULL) +
      theme_bw()
    
    cat("SUCCESS: VIP plot with directionality created successfully!\n")
    cat("Saving test plot to test_vip_directionality.png\n")
    ggsave("test_vip_directionality.png", vip_plot, width = 6, height = 4, dpi = 150)
    cat("Test passed!\n")
  }, error = function(e) {
    cat("ERROR creating plot:", e$message, "\n")
  })
} else {
  cat("SKIP: Cannot test plotting without directionality\n")
}

cat("\n=== Test Complete ===\n")

