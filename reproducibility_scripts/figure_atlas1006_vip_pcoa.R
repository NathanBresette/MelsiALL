#!/usr/bin/env Rscript
# ==============================================================================
# Atlas1006 VIP and PCoA Figure Generation
# ==============================================================================
# Creates a VIP bar plot and PCoA visualization for the Atlas1006 dataset using
# the MeLSI learned metric. Outputs are saved to the figures/ directory.
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(vegan)
  library(scales)
  library(gridExtra)
})

set.seed(42)

# Determine consistent project-relative paths ---------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
project_dir <- dirname(script_dir)
output_dir <- file.path(project_dir, "figures")

# Source MeLSI functions directly to ensure we have the latest version with directionality
# This ensures directionality is always available
melsi_source_path <- file.path(project_dir, "github/R/melsi_robust.R")
if (file.exists(melsi_source_path)) {
  source(melsi_source_path)
  cat("Loaded MeLSI functions from source file\n")
} else {
  # Fallback to installed package
  library(MeLSI)
  cat("Using installed MeLSI package (source file not found)\n")
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n")
cat("==============================================================================\n")
cat("Generating VIP and PCoA figures for Atlas1006 dataset\n")
cat("==============================================================================\n\n")

# Load Atlas1006 dataset -------------------------------------------------------
data(atlas1006)

counts <- t(abundances(atlas1006))
metadata <- meta(atlas1006)
sex <- metadata$sex

# Remove samples with missing labels
valid_idx <- !is.na(sex)
counts <- counts[valid_idx, , drop = FALSE]
sex <- sex[valid_idx]

cat("Samples after filtering:", length(sex), "\n")
cat("Number of taxa:", ncol(counts), "\n\n")

# CLR transformation -----------------------------------------------------------
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)

# Run MeLSI --------------------------------------------------------------------
cat("Running MeLSI (n_perms = 200)...\n")
cat("Progress will be shown for each permutation.\n\n")
flush.console()

melsi_res <- melsi(
  counts_clr,
  sex,
  n_perms = 200,
  B = 30,
  m_frac = 0.8,
  show_progress = TRUE,
  plot_vip = FALSE
)

cat("\nMeLSI Results:\n")
cat("  F-statistic:", round(melsi_res$F_observed, 3), "\n")
cat("  P-value    :", round(melsi_res$p_value, 3), "\n")

# Check if directionality is available
has_directionality <- !is.null(melsi_res$directionality) && length(melsi_res$directionality) > 0

if (has_directionality) {
  cat("  Directionality: Available\n")
  cat("  Directionality values:", length(melsi_res$directionality), "features\n")
  # Show unique directionality values
  unique_dirs <- unique(melsi_res$directionality)
  cat("  Unique directions:", paste(unique_dirs, collapse = ", "), "\n")
} else {
  cat("  Directionality: Not available (this should not happen for 2-group analysis)\n")
  cat("  Available result names:", paste(names(melsi_res), collapse = ", "), "\n")
}
cat("\n")

# Conservative pre-filtering to align with MeLSI feature weights ---------------
counts_filtered <- MeLSI:::apply_conservative_prefiltering(counts_clr, sex, filter_frac = 0.7)

# Create consistent color scheme for both plots
groups <- unique(sex)
# Use ggplot2 default colors to match PCoA plot
group_colors <- scales::hue_pal()(length(groups))
names(group_colors) <- groups

# VIP Plot Data Preparation -------------------------------------------------
vip_df <- data.frame(
  Feature = names(melsi_res$feature_weights),
  Weight = as.numeric(melsi_res$feature_weights)
) %>%
  arrange(desc(Weight)) %>%
  slice_head(n = 15) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

# Create VIP plot WITHOUT directionality (left panel)
vip_plot_no_dir <- ggplot(vip_df, aes(x = Weight, y = Feature)) +
  geom_col(fill = "#1f77b4") +
  geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 3) +
  labs(
    title = "Top Atlas1006 Taxa Weighted by MeLSI",
    x = "Feature Weight",
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

# Create VIP plot WITH directionality (right panel)
if (has_directionality) {
  # Match directionality to features
  vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$Feature)]
  
  # Handle any NAs
  na_idx <- is.na(vip_df$Directionality)
  if (any(na_idx)) {
    feature_names_in_directionality <- names(melsi_res$directionality)
    for (i in which(na_idx)) {
      feature_name <- as.character(vip_df$Feature[i])
      if (feature_name %in% feature_names_in_directionality) {
        vip_df$Directionality[i] <- melsi_res$directionality[feature_name]
      }
    }
  }
  
  vip_df$Directionality <- ifelse(is.na(vip_df$Directionality), "Unknown", vip_df$Directionality)
  unique_directions <- unique(vip_df$Directionality[vip_df$Directionality != "Unknown"])
  
  if (length(unique_directions) >= 2) {
    # Map directionality to group colors
    color_map <- character(length(unique_directions))
    names(color_map) <- unique_directions
    
    for (dir in unique_directions) {
      group_name <- gsub("Higher in ", "", dir)
      if (group_name %in% names(group_colors)) {
        color_map[dir] <- group_colors[group_name]
      } else {
        dir_idx <- which(unique_directions == dir)
        if (dir_idx <= length(group_colors)) {
          color_map[dir] <- group_colors[dir_idx]
        } else {
          color_map[dir] <- "#95A5A6"
        }
      }
    }
    
    if (any(vip_df$Directionality == "Unknown")) {
      color_map["Unknown"] <- "#95A5A6"
  }
  
    vip_plot_with_dir <- ggplot(vip_df, aes(x = Weight, y = Feature, fill = Directionality)) +
    geom_col() +
      scale_fill_manual(values = color_map, name = "Directionality", drop = FALSE) +
    geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 3) +
    labs(
      title = "Top Atlas1006 Taxa Weighted by MeLSI",
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
} else {
    # Fallback if only one direction
    vip_plot_with_dir <- vip_plot_no_dir
  }
} else {
  # If no directionality, use same plot for both
  vip_plot_with_dir <- vip_plot_no_dir
}

# Save individual plots
vip_path_no_dir <- file.path(output_dir, "atlas1006_vip_no_directionality.png")
vip_path_with_dir <- file.path(output_dir, "atlas1006_vip.png")
vip_path_combined <- file.path(output_dir, "atlas1006_vip_combined.png")

ggsave(vip_path_no_dir, vip_plot_no_dir, width = 6, height = 4, dpi = 300)
cat("VIP plot (no directionality) saved to:", vip_path_no_dir, "\n")

ggsave(vip_path_with_dir, vip_plot_with_dir, width = 6, height = 4, dpi = 300)
cat("VIP plot (with directionality) saved to:", vip_path_with_dir, "\n")

# Create combined side-by-side plot
combined_vip <- grid.arrange(
  vip_plot_no_dir,
  vip_plot_with_dir,
  ncol = 2,
  widths = c(1, 1.15)  # Slightly wider right panel to accommodate legend
)

ggsave(vip_path_combined, combined_vip, width = 12, height = 4, dpi = 300)
cat("Combined VIP plot (side-by-side) saved to:", vip_path_combined, "\n")

# PCoA Plot using MeLSI learned metric -----------------------------------------
cat("\nComputing PCoA using MeLSI learned distance...\n")
dist_melsi <- MeLSI:::calculate_mahalanobis_dist_robust(counts_filtered, melsi_res$metric_matrix)
pcoa_res <- ape::pcoa(dist_melsi)

explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Group <- sex

# Use same color scheme as VIP plot for consistency
pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.4) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "Atlas1006 PCoA with MeLSI Distance",
    x = paste0("PCoA1 (", explained[1], "%)"),
    y = paste0("PCoA2 (", explained[2], "%)"),
    color = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

pcoa_path <- file.path(output_dir, "atlas1006_pcoa.png")
ggsave(pcoa_path, pcoa_plot, width = 6, height = 4, dpi = 300)
cat("PCoA plot saved to:", pcoa_path, "\n")

cat("\nFigure generation complete.\n")
cat("==============================================================================\n\n")

