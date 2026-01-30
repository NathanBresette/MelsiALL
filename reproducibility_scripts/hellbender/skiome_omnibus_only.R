#!/usr/bin/env Rscript
# ==============================================================================
# SKIOME Validation: Omnibus-Only MeLSI Analysis
# ==============================================================================
# Validates MeLSI on SKIOME skin microbiome data (PRJNA554499)
# Three groups: Atopic_Dermatitis, Healthy, Psoriasis
# Only runs omnibus test (no pairwise comparisons)
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("SKIOME Omnibus-Only Validation: PRJNA554499\n")
cat("Groups: Atopic_Dermatitis, Healthy, Psoriasis\n")
cat("==============================================================================\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(vegan)
  library(MeLSI)
  library(GUniFrac)
  library(ape)
  library(ggplot2)
})

# Use PostScript device for headless operation (required on HPC)
cat("Using PostScript device for figure generation (headless-compatible)\n")

# Load data
cat("Loading SKIOME data...\n")
load("skiome_data_loaded.RData")

cat("Data loaded:\n")
cat("  Samples:", nrow(counts), "\n")
cat("  Taxa:", ncol(counts), "\n")
cat("  Groups:", paste(unique(metadata$Group), collapse = ", "), "\n")
cat("\nGroup distribution:\n")
print(table(metadata$Group))

# Check if we have taxa names
original_taxa_names <- colnames(counts)
if (is.null(original_taxa_names) || all(grepl("^Taxa_", original_taxa_names))) {
  cat("\n⚠️  Warning: No meaningful taxa names found, using generic names\n")
  has_taxa_names <- FALSE
} else {
  cat("\n✓ Taxa names preserved from original data\n")
  has_taxa_names <- TRUE
  cat("  First 5 taxa names:", paste(head(original_taxa_names, 5), collapse = ", "), "\n")
}

# Check group sizes
group_counts <- table(metadata$Group)
if (any(group_counts < 5)) {
  cat("\n⚠️  Warning: Some groups have < 5 samples\n")
  print(group_counts)
}

# CLR transformation
cat("\nApplying CLR transformation...\n")
X_clr <- log(counts + 1)
X_clr <- X_clr - rowMeans(X_clr)

# Preserve taxa names if available
if (has_taxa_names && length(original_taxa_names) == ncol(X_clr)) {
  colnames(X_clr) <- original_taxa_names
  cat("✓ Preserved original taxa names in CLR-transformed data\n")
} else {
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  cat("⚠️  Using generic taxa names\n")
}

# Get group labels
groups <- metadata$Group

# Remove any samples with missing groups
valid_idx <- !is.na(groups)
X_clr <- X_clr[valid_idx, , drop = FALSE]
groups <- groups[valid_idx]

cat("\nAfter filtering:\n")
cat("  Samples:", nrow(X_clr), "\n")
cat("  Groups:", paste(unique(groups), collapse = ", "), "\n")
cat("\nFinal group distribution:\n")
print(table(groups))

# Run MeLSI (omnibus only - no pairwise)
cat("\n")
cat("==============================================================================\n")
cat("Running MeLSI (Omnibus Analysis Only)\n")
cat("==============================================================================\n")
cat("Running omnibus test only (no pairwise comparisons).\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()
# Use analysis_type = "omnibus" to skip pairwise comparisons
melsi_result <- melsi(
  X_clr,
  groups,
  analysis_type = "omnibus",  # Only omnibus, no pairwise
  n_perms = 200,
  B = 30,
  m_frac = 0.8,
  show_progress = TRUE,
  plot_vip = FALSE
)
melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\nMeLSI Results:\n")
# For omnibus-only, results should be in omnibus slot or directly
if ("omnibus" %in% names(melsi_result)) {
  melsi_omnibus <- melsi_result$omnibus
  cat("  Omnibus F-statistic:", round(melsi_omnibus$F_observed, 3), "\n")
  cat("  Omnibus P-value    :", round(melsi_omnibus$p_value, 4), "\n")
} else {
  # Fallback to standard format
  melsi_omnibus <- melsi_result
  cat("  F-statistic:", round(melsi_omnibus$F_observed, 3), "\n")
  cat("  P-value    :", round(melsi_omnibus$p_value, 4), "\n")
}
cat("  Time       :", round(melsi_time, 1), "seconds\n")

# Run traditional methods for comparison
cat("\n")
cat("==============================================================================\n")
cat("Running Traditional Methods for Comparison\n")
cat("==============================================================================\n\n")

# Euclidean PERMANOVA
cat("Euclidean PERMANOVA...\n")
start_time <- Sys.time()
dist_euc <- dist(X_clr)
perm_euc <- adonis2(dist_euc ~ groups, permutations = 999)
euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  F-statistic:", round(perm_euc$F[1], 3), "\n")
cat("  P-value    :", round(perm_euc$`Pr(>F)`[1], 4), "\n")
cat("  Time       :", round(euc_time, 1), "seconds\n\n")

# Bray-Curtis PERMANOVA
cat("Bray-Curtis PERMANOVA...\n")
start_time <- Sys.time()
dist_bray <- vegdist(counts[valid_idx, ], method = "bray")
perm_bray <- adonis2(dist_bray ~ groups, permutations = 999)
bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  F-statistic:", round(perm_bray$F[1], 3), "\n")
cat("  P-value    :", round(perm_bray$`Pr(>F)`[1], 4), "\n")
cat("  Time       :", round(bray_time, 1), "seconds\n\n")

# Jaccard PERMANOVA
cat("Jaccard PERMANOVA...\n")
start_time <- Sys.time()
dist_jaccard <- vegdist(counts[valid_idx, ], method = "jaccard")
perm_jaccard <- adonis2(dist_jaccard ~ groups, permutations = 999)
jaccard_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  F-statistic:", round(perm_jaccard$F[1], 3), "\n")
cat("  P-value    :", round(perm_jaccard$`Pr(>F)`[1], 4), "\n")
cat("  Time       :", round(jaccard_time, 1), "seconds\n\n")

# Compile results
melsi_F <- melsi_omnibus$F_observed
melsi_p <- melsi_omnibus$p_value

results <- data.frame(
  Method = c("MeLSI", "Euclidean", "Bray-Curtis", "Jaccard"),
  F_statistic = c(
    melsi_F,
    perm_euc$F[1],
    perm_bray$F[1],
    perm_jaccard$F[1]
  ),
  P_value = c(
    melsi_p,
    perm_euc$`Pr(>F)`[1],
    perm_bray$`Pr(>F)`[1],
    perm_jaccard$`Pr(>F)`[1]
  ),
  Time_sec = c(melsi_time, euc_time, bray_time, jaccard_time),
  Significant = c(
    melsi_p < 0.05,
    perm_euc$`Pr(>F)`[1] < 0.05,
    perm_bray$`Pr(>F)`[1] < 0.05,
    perm_jaccard$`Pr(>F)`[1] < 0.05
  ),
  stringsAsFactors = FALSE
)

# Save results
write.csv(results, "skiome_omnibus_results.csv", row.names = FALSE)
cat("Results saved to: skiome_omnibus_results.csv\n\n")

cat("==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n")
print(results)
cat("\n")

# Feature importance (if available)
if (!is.null(melsi_omnibus$feature_weights) && length(melsi_omnibus$feature_weights) > 0) {
  cat("Top 10 Feature Weights:\n")
  top_features <- head(sort(melsi_omnibus$feature_weights, decreasing = TRUE), 10)
  print(top_features)
  cat("\n")
}

# Generate and save plots
cat("==============================================================================\n")
cat("Generating Plots\n")
cat("==============================================================================\n\n")

# Generate VIP plot with 3-group directionality
cat("Generating VIP plot (3-group directionality enabled)...\n")
tryCatch({
  if (!is.null(melsi_omnibus$feature_weights) && length(melsi_omnibus$feature_weights) > 0) {
    # For 3+ groups, directionality shows which group has highest abundance for each feature
    # Note: plot_vip currently only colors for 2 groups, but directionality info is still useful
    # We'll enable it - the function will show directionality labels even if colors are limited
    vip_plot <- plot_vip(melsi_omnibus, top_n = 15, 
                         title = "SKIOME: Variable Importance (MeLSI)",
                         directionality = TRUE)  # Enable directionality - shows which group has highest abundance
    
    postscript("skiome_omnibus_vip_plot.ps", width = 10, height = 8, 
               horizontal = FALSE, onefile = FALSE, paper = "special")
    print(vip_plot)
    dev.off()
    
    if (file.exists("skiome_omnibus_vip_plot.ps")) {
      cat("  ✓ VIP plot saved to: skiome_omnibus_vip_plot.ps\n")
    } else {
      cat("  ✗ VIP plot file was not created\n")
    }
  } else {
    cat("  ⚠️  No feature weights available for VIP plot\n")
  }
}, error = function(e) {
  cat("  ✗ Error generating VIP plot:", e$message, "\n")
  traceback()
})

# Generate PCoA plots using multiple methods (since job takes long)
cat("Generating PCoA plots (trying multiple methods)...\n")
tryCatch({
  if (!is.null(melsi_omnibus$distance_matrix)) {
    dist_matrix <- melsi_omnibus$distance_matrix
    
    # Method 1: Use MeLSI's plot_pcoa function (cmdscale) - fix alpha for PostScript
    cat("  Method 1: cmdscale (MeLSI default)...\n")
    pcoa_plot1 <- plot_pcoa(melsi_omnibus, X_clr, groups, 
                            title = "SKIOME: PCoA using MeLSI Distance (cmdscale)")
    
    # Fix alpha transparency for PostScript: replace point layer
    pcoa_plot1_fixed <- pcoa_plot1
    # Remove first layer (points with alpha) and add solid points
    pcoa_plot1_fixed$layers <- pcoa_plot1_fixed$layers[-1]  # Remove point layer
    pcoa_plot1_fixed <- pcoa_plot1_fixed + 
      ggplot2::geom_point(size = 3, alpha = 1.0)  # Add solid points
    
    postscript("skiome_omnibus_pcoa_cmdscale.ps", width = 10, height = 8, 
               horizontal = FALSE, onefile = FALSE, paper = "special")
    print(pcoa_plot1_fixed)
    dev.off()
    cat("    ✓ Saved: skiome_omnibus_pcoa_cmdscale.ps\n")
    
    # Method 2: Try ape::pcoa
    cat("  Method 2: ape::pcoa...\n")
    tryCatch({
      pcoa_result2 <- ape::pcoa(dist_matrix)
      var_explained2 <- pcoa_result2$values$Relative_eig[1:2] * 100
      
      plot_data2 <- data.frame(
        PC1 = pcoa_result2$vectors[, 1],
        PC2 = pcoa_result2$vectors[, 2],
        Group = as.factor(groups)
      )
      
      pcoa_plot2 <- ggplot2::ggplot(plot_data2, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
        ggplot2::geom_point(size = 3, alpha = 1.0) +
        ggplot2::stat_ellipse(level = 0.95, linetype = 2, alpha = 1.0) +
        ggplot2::labs(
          title = "SKIOME: PCoA using MeLSI Distance (ape::pcoa)",
          x = sprintf("PCoA1 (%.1f%%)", var_explained2[1]),
          y = sprintf("PCoA2 (%.1f%%)", var_explained2[2])
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right"
        )
      
      postscript("skiome_omnibus_pcoa_ape.ps", width = 10, height = 8, 
                 horizontal = FALSE, onefile = FALSE, paper = "special")
      print(pcoa_plot2)
      dev.off()
      cat("    ✓ Saved: skiome_omnibus_pcoa_ape.ps\n")
    }, error = function(e) {
      cat("    ✗ ape::pcoa failed:", e$message, "\n")
    })
    
    # Copy main one with standard name
    file.copy("skiome_omnibus_pcoa_cmdscale.ps", "skiome_omnibus_pcoa_plot.ps", overwrite = TRUE)
    cat("  ✓ Main plot saved as: skiome_omnibus_pcoa_plot.ps\n")
    
    if (file.exists("skiome_omnibus_pcoa_plot.ps")) {
      cat("  ✓ PCoA plot saved to: skiome_omnibus_pcoa_plot.ps\n")
    } else {
      cat("  ✗ PCoA plot file was not created\n")
    }
  } else {
    cat("  ⚠️  No distance matrix available for PCoA plot\n")
  }
}, error = function(e) {
  cat("  ✗ Error generating PCoA plot:", e$message, "\n")
  traceback()
})

cat("\n")
cat("==============================================================================\n")
cat("Omnibus-Only Validation Complete!\n")
cat("==============================================================================\n\n")
