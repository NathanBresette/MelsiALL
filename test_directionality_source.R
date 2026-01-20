#!/usr/bin/env Rscript
# Test directionality by sourcing the file directly

suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(dplyr)
  library(vegan)
})

# Source the MeLSI functions directly
source("github/R/melsi_robust.R")

set.seed(42)

cat("Testing directionality with sourced functions...\n\n")

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
  cat("✓ SUCCESS: Directionality is available!\n")
  cat("Directionality length:", length(melsi_res$directionality), "\n")
  cat("Unique directions:", paste(unique(melsi_res$directionality), collapse = ", "), "\n")
  cat("First 5 directionality values:\n")
  print(head(melsi_res$directionality, 5))
  
  # Create and save VIP plot with directionality
  cat("\n=== Creating VIP Plot with Directionality ===\n")
  vip_df <- data.frame(
    Feature = names(melsi_res$feature_weights),
    Weight = as.numeric(melsi_res$feature_weights)
  ) %>%
    arrange(desc(Weight)) %>%
    slice_head(n = 10) %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)))
  
  # Add directionality
  vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$Feature)]
  vip_df$Directionality <- ifelse(is.na(vip_df$Directionality), "Unknown", vip_df$Directionality)
  
  unique_dirs <- unique(vip_df$Directionality)
  cat("Unique directions in plot:", paste(unique_dirs, collapse = ", "), "\n")
  
  # Create color map
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
    geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 3) +
    labs(
      title = "Test VIP Plot with Directionality",
      x = "Feature Weight",
      y = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      legend.position = "right"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # Save the plot
  ggsave("test_vip_directionality.png", vip_plot, width = 8, height = 5, dpi = 300)
  cat("✓ VIP plot saved to: test_vip_directionality.png\n")
  cat("\n✓ Test PASSED: Directionality is working correctly!\n")
} else {
  cat("✗ ERROR: Directionality not found in results!\n")
  cat("Available result names:", paste(names(melsi_res), collapse = ", "), "\n")
  cat("\n✗ Test FAILED\n")
}

cat("\n=== Test Complete ===\n")

