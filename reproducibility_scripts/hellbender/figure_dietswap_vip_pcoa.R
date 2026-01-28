#!/usr/bin/env Rscript
# ==============================================================================
# DietSwap VIP and PCoA Figure Generation
# ==============================================================================
# Creates a VIP bar plot and PCoA visualization for the DietSwap dataset using
# the MeLSI learned metric. Outputs are saved to the figures/ directory.
# Matches the format and style of Atlas1006 figures.
# ==============================================================================

# Load required packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(vegan)
  library(scales)
  library(phyloseq)
})

# Install and load gridExtra if needed
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  cat("Installing gridExtra package...\n")
  install.packages("gridExtra", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(gridExtra)
cat("gridExtra package loaded\n")

# Use R's built-in PostScript device (no packages needed, works on headless systems without Cairo)
cat("Using R's built-in PostScript device for output\n")

set.seed(42)

# Use PostScript output (works on headless systems without Cairo/X11)
cat("Using PostScript output format (headless-compatible)\n")

# Determine output directory - save to current directory or figures/ if it exists
output_dir <- "figures"
if (!dir.exists(output_dir)) {
  # Try to create figures directory, or use current directory
  tryCatch({
    dir.create(output_dir, recursive = TRUE)
  }, error = function(e) {
    output_dir <<- "."
  })
}
if (!dir.exists(output_dir)) {
  output_dir <- "."
}

# Try to source MeLSI functions, otherwise use installed package
melsi_source_path <- "../../github/R/melsi_robust.R"  # Relative to hellbender/
if (file.exists(melsi_source_path)) {
  source(melsi_source_path)
  cat("Loaded MeLSI functions from source file\n")
} else {
  # Try absolute path from project root
  melsi_source_path2 <- "~/melsi_simulations/github/R/melsi_robust.R"
  if (file.exists(melsi_source_path2)) {
    source(melsi_source_path2)
    cat("Loaded MeLSI functions from absolute path\n")
  } else {
    # Fallback to installed package
    library(MeLSI)
    cat("Using installed MeLSI package\n")
    # Define helper functions if not available
    if (!exists("apply_conservative_prefiltering")) {
      apply_conservative_prefiltering <- function(X, y, filter_frac = 0.7) {
        group1_idx <- which(y == unique(y)[1])
        group2_idx <- which(y == unique(y)[2])
        importance_scores <- numeric(ncol(X))
        for (j in 1:ncol(X)) {
          mu1 <- mean(X[group1_idx, j])
          mu2 <- mean(X[group2_idx, j])
          sigma1_sq <- var(X[group1_idx, j])
          sigma2_sq <- var(X[group2_idx, j])
          denominator <- sqrt(sigma1_sq + sigma2_sq + 1e-6)
          importance_scores[j] <- abs(mu1 - mu2) / denominator
        }
        n_keep <- max(1, floor(filter_frac * ncol(X)))
        top_features <- order(importance_scores, decreasing = TRUE)[1:n_keep]
        return(X[, top_features, drop = FALSE])
      }
    }
    if (!exists("calculate_mahalanobis_dist_robust")) {
      calculate_mahalanobis_dist_robust <- function(X, M) {
        # Compute Mahalanobis distance: sqrt((x_i - x_j)^T M (x_i - x_j))
        n <- nrow(X)
        dist_matrix <- matrix(0, nrow = n, ncol = n)
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            diff <- X[i, ] - X[j, ]
            dist_sq <- as.numeric(t(diff) %*% M %*% diff)
            dist_matrix[i, j] <- sqrt(max(0, dist_sq))
            dist_matrix[j, i] <- dist_matrix[i, j]
          }
        }
        return(as.dist(dist_matrix))
      }
    }
  }
}

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("\n")
cat("==============================================================================\n")
cat("Generating VIP and PCoA figures for DietSwap dataset\n")
cat("==============================================================================\n\n")

# Load DietSwap dataset -------------------------------------------------------
cat("Loading DietSwap dataset...\n")
data(dietswap)

# Subset to baseline samples (timepoint.within.group == 1) and two groups (DI, HE)
dietswap_subset <- subset_samples(dietswap, timepoint.within.group == 1)
dietswap_subset <- subset_samples(dietswap_subset, group %in% c("DI", "HE"))
dietswap_subset <- prune_taxa(taxa_sums(dietswap_subset) > 0, dietswap_subset)

# Extract counts and metadata
counts <- as(otu_table(dietswap_subset), "matrix")
if (taxa_are_rows(dietswap_subset)) {
  counts <- t(counts)
}
metadata <- data.frame(sample_data(dietswap_subset))
metadata$group <- droplevels(metadata$group)
counts <- counts[rownames(metadata), ]

# Get group labels
groups <- metadata$group

# Remove samples with missing labels
valid_idx <- !is.na(groups)
counts <- counts[valid_idx, , drop = FALSE]
groups <- groups[valid_idx]

cat("Samples after filtering:", length(groups), "\n")
cat("Number of taxa:", ncol(counts), "\n")
cat("Groups:", paste(unique(groups), collapse = ", "), "\n")
cat("Group distribution:\n")
print(table(groups))
cat("\n")

# CLR transformation -----------------------------------------------------------
cat("Applying CLR transformation...\n")
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)
colnames(counts_clr) <- colnames(counts)

# Run MeLSI --------------------------------------------------------------------
cat("Running MeLSI (n_perms = 200)...\n")
cat("Progress will be shown for each permutation.\n\n")
flush.console()

melsi_res <- melsi(
  counts_clr,
  groups,
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
} else {
  cat("  Directionality: Not available (this should not happen for 2-group analysis)\n")
}
cat("\n")

# Conservative pre-filtering to align with MeLSI feature weights ---------------
# Try MeLSI::: first, then direct call if sourced
if (exists("apply_conservative_prefiltering", mode = "function")) {
  counts_filtered <- apply_conservative_prefiltering(counts_clr, groups, filter_frac = 0.7)
} else {
  counts_filtered <- MeLSI:::apply_conservative_prefiltering(counts_clr, groups, filter_frac = 0.7)
}

# Create consistent color scheme for both plots
group_names <- unique(groups)
# Use ggplot2 default colors to match PCoA plot
group_colors <- scales::hue_pal()(length(group_names))
names(group_colors) <- group_names

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
    title = "Top DietSwap Taxa Weighted by MeLSI",
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
        title = "Top DietSwap Taxa Weighted by MeLSI",
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

# Save individual plots as PostScript (works on headless systems without Cairo/X11)
vip_path_no_dir <- file.path(output_dir, "dietswap_vip_no_directionality.ps")
vip_path_with_dir <- file.path(output_dir, "dietswap_vip.ps")
vip_path_combined <- file.path(output_dir, "dietswap_vip_combined.ps")

# Save plots using R's built-in PostScript device (works on headless systems)
tryCatch({
  postscript(vip_path_no_dir, width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special")
  print(vip_plot_no_dir)
  dev.off()
  if (file.exists(vip_path_no_dir)) {
    cat("VIP plot (no directionality) saved to:", vip_path_no_dir, "\n")
  } else {
    stop("File was not created")
  }
}, error = function(e) {
  cat("Error saving VIP plot (no directionality):", e$message, "\n")
  stop("Failed to save VIP plot: ", e$message)
})

tryCatch({
  postscript(vip_path_with_dir, width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special")
  print(vip_plot_with_dir)
  dev.off()
  if (file.exists(vip_path_with_dir)) {
    cat("VIP plot (with directionality) saved to:", vip_path_with_dir, "\n")
  } else {
    stop("File was not created")
  }
}, error = function(e) {
  cat("Error saving VIP plot (with directionality):", e$message, "\n")
  stop("Failed to save VIP plot: ", e$message)
})

# Create combined side-by-side plot using arrangeGrob (works better with PostScript)
# Load grid package for grid.draw()
library(grid)
tryCatch({
  combined_vip <- arrangeGrob(
    vip_plot_no_dir,
    vip_plot_with_dir,
    ncol = 2,
    widths = c(1, 1.15)  # Slightly wider right panel to accommodate legend
  )
  
  postscript(vip_path_combined, width = 12, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special")
  grid.draw(combined_vip)
  dev.off()
  if (file.exists(vip_path_combined)) {
    cat("Combined VIP plot (side-by-side) saved to:", vip_path_combined, "\n")
  } else {
    stop("File was not created")
  }
}, error = function(e) {
  cat("Error saving combined VIP plot:", e$message, "\n")
  stop("Failed to save combined VIP plot: ", e$message)
})

# PCoA Plot using MeLSI learned metric -----------------------------------------
cat("\nComputing PCoA using MeLSI learned distance...\n")
# Try direct call first, then MeLSI:::
if (exists("calculate_mahalanobis_dist_robust", mode = "function")) {
  dist_melsi <- calculate_mahalanobis_dist_robust(counts_filtered, melsi_res$metric_matrix)
} else {
  dist_melsi <- MeLSI:::calculate_mahalanobis_dist_robust(counts_filtered, melsi_res$metric_matrix)
}
pcoa_res <- ape::pcoa(dist_melsi)

explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Group <- groups

# Use same color scheme as VIP plot for consistency
# Note: PostScript doesn't support transparency well, so use solid points
pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 5.5, alpha = 1.0, stroke = 1.0, shape = 19) +  # Large, solid, visible points (shape 19 = filled circle)
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "DietSwap PCoA with MeLSI Distance",
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

pcoa_path <- file.path(output_dir, "dietswap_pcoa.ps")
tryCatch({
  postscript(pcoa_path, width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special")
  print(pcoa_plot)
  dev.off()
  if (file.exists(pcoa_path)) {
    cat("PCoA plot saved to:", pcoa_path, "\n")
  } else {
    stop("File was not created")
  }
}, error = function(e) {
  cat("Error saving PCoA plot:", e$message, "\n")
  stop("Failed to save PCoA plot: ", e$message)
})

cat("\nFigure generation complete.\n")
cat("==============================================================================\n\n")
