#!/usr/bin/env Rscript
# ==============================================================================
# SKIOME Validation: Multi-Group MeLSI Analysis
# ==============================================================================
# Validates MeLSI on SKIOME skin microbiome data (PRJNA554499)
# Three groups: Atopic_Dermatitis, Healthy, Psoriasis
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("SKIOME Multi-Group Validation: PRJNA554499\n")
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
# PostScript works without Cairo/X11 dependencies
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
colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))

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

# Run MeLSI (multi-group)
cat("\n")
cat("==============================================================================\n")
cat("Running MeLSI (Multi-Group Analysis)\n")
cat("==============================================================================\n")
cat("This will run both omnibus (all groups) and pairwise comparisons.\n")
cat("This may take several minutes...\n\n")

start_time <- Sys.time()
# Use analysis_type = "both" for 3+ groups (omnibus + pairwise)
# Or "auto" will automatically choose "both" for 3+ groups
melsi_result <- melsi(
  X_clr,
  groups,
  analysis_type = "auto",  # Will automatically use "both" for 3 groups
  n_perms = 200,
  B = 30,
  m_frac = 0.8,
  show_progress = TRUE,
  plot_vip = FALSE
)
melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat("\nMeLSI Results:\n")
# For multi-group, results may have omnibus and pairwise components
if ("omnibus" %in% names(melsi_result)) {
  cat("  Omnibus F-statistic:", round(melsi_result$omnibus$F_observed, 3), "\n")
  cat("  Omnibus P-value    :", round(melsi_result$omnibus$p_value, 4), "\n")
  if ("pairwise" %in% names(melsi_result)) {
    cat("\n  Pairwise Comparisons:\n")
    if ("summary" %in% names(melsi_result$pairwise)) {
      print(melsi_result$pairwise$summary)
    }
  }
} else {
  # Fallback to standard format
  cat("  F-statistic:", round(melsi_result$F_observed, 3), "\n")
  cat("  P-value    :", round(melsi_result$p_value, 4), "\n")
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

# Weighted UniFrac (if we have a tree - for now skip or use mock tree)
cat("Weighted UniFrac PERMANOVA...\n")
cat("  (Skipping - requires phylogenetic tree)\n\n")

# Unweighted UniFrac
cat("Unweighted UniFrac PERMANOVA...\n")
cat("  (Skipping - requires phylogenetic tree)\n\n")

# Compile results - handle multi-group MeLSI results
if ("omnibus" %in% names(melsi_result)) {
  melsi_F <- melsi_result$omnibus$F_observed
  melsi_p <- melsi_result$omnibus$p_value
} else {
  melsi_F <- melsi_result$F_observed
  melsi_p <- melsi_result$p_value
}

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
write.csv(results, "skiome_validation_results.csv", row.names = FALSE)
cat("Results saved to: skiome_validation_results.csv\n\n")

cat("==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n")
print(results)
cat("\n")

# Feature importance (if available)
if (!is.null(melsi_result$feature_weights) && length(melsi_result$feature_weights) > 0) {
  cat("Top 10 Feature Weights:\n")
  top_features <- head(sort(melsi_result$feature_weights, decreasing = TRUE), 10)
  print(top_features)
  cat("\n")
}

# Generate and save plots
cat("==============================================================================\n")
cat("Generating Plots\n")
cat("==============================================================================\n\n")

# Extract omnibus results for plotting (multi-group analysis)
if ("omnibus" %in% names(melsi_result)) {
  melsi_omnibus <- melsi_result$omnibus
} else {
  melsi_omnibus <- melsi_result
}

# Generate VIP plot using PostScript (works on headless systems)
cat("Generating VIP plot...\n")
tryCatch({
  if (!is.null(melsi_omnibus$feature_weights) && length(melsi_omnibus$feature_weights) > 0) {
    vip_plot <- plot_vip(melsi_omnibus, top_n = 15, 
                         title = "SKIOME: Variable Importance (MeLSI)")
    postscript("skiome_vip_plot.ps", width = 10, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special")
    print(vip_plot)
    dev.off()
    if (file.exists("skiome_vip_plot.ps")) {
      cat("  ✓ VIP plot saved to: skiome_vip_plot.ps\n")
    } else {
      cat("  ✗ VIP plot file was not created\n")
    }
  } else {
    cat("  ⚠️  No feature weights available for VIP plot\n")
  }
}, error = function(e) {
  cat("  ✗ Error generating VIP plot:", e$message, "\n")
})

# Generate PCoA plot using PostScript (works on headless systems)
# Create plot with solid points for PostScript compatibility
cat("Generating PCoA plot...\n")
tryCatch({
  if (!is.null(melsi_omnibus$distance_matrix)) {
    # Get distance matrix and compute PCoA
    dist_matrix <- melsi_omnibus$distance_matrix
    pcoa_result <- stats::cmdscale(dist_matrix, k = 2, eig = TRUE)
    
    # Calculate variance explained
    var_explained <- pcoa_result$eig / sum(abs(pcoa_result$eig)) * 100
    
    # Create plot data
    plot_data <- data.frame(
      PC1 = pcoa_result$points[, 1],
      PC2 = pcoa_result$points[, 2],
      Group = as.factor(groups)
    )
    
    # Create ggplot with solid points (PostScript doesn't support transparency well)
    pcoa_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3.5) +  # Solid points, no alpha for PostScript
      stat_ellipse(level = 0.95, linetype = 2, linewidth = 1.2) +
      labs(
        title = "SKIOME: PCoA using MeLSI Distance",
        x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
        y = sprintf("PCoA2 (%.1f%%)", var_explained[2]),
        color = "Group"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
      )
    
    postscript("skiome_pcoa_plot.ps", width = 10, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special")
    print(pcoa_plot)
    dev.off()
    if (file.exists("skiome_pcoa_plot.ps")) {
      cat("  ✓ PCoA plot saved to: skiome_pcoa_plot.ps\n")
    } else {
      cat("  ✗ PCoA plot file was not created\n")
    }
  } else {
    cat("  ⚠️  No distance matrix available for PCoA plot\n")
  }
}, error = function(e) {
  cat("  ✗ Error generating PCoA plot:", e$message, "\n")
})

cat("\n")

cat("==============================================================================\n")
cat("Validation Complete!\n")
cat("==============================================================================\n\n")
