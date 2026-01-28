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

# Generate PCoA plot using MeLSI's plot_pcoa function (like DietSwap)
cat("Generating PCoA plot...\n")
tryCatch({
  if (!is.null(melsi_omnibus$distance_matrix)) {
    # Use MeLSI's plot_pcoa function (consistent with DietSwap)
    pcoa_plot <- plot_pcoa(melsi_omnibus, X_clr, groups, 
                           title = "SKIOME: PCoA using MeLSI Distance")
    
    postscript("skiome_pcoa_plot.ps", width = 10, height = 8, 
               horizontal = FALSE, onefile = FALSE, paper = "special")
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

# Generate pairwise plots (VIP and PCoA for each comparison)
if ("pairwise" %in% names(melsi_result) && "pairwise_results" %in% names(melsi_result$pairwise)) {
  cat("\n")
  cat("==============================================================================\n")
  cat("Generating Pairwise Comparison Plots\n")
  cat("==============================================================================\n\n")
  
  pairwise_results <- melsi_result$pairwise$pairwise_results
  
  for (comparison_name in names(pairwise_results)) {
    pair_result <- pairwise_results[[comparison_name]]
    
    if (is.null(pair_result)) {
      cat("Skipping", comparison_name, "- no results available\n")
      next
    }
    
    # Extract group names from comparison name (e.g., "Atopic_Dermatitis_vs_Psoriasis")
    group_names <- strsplit(comparison_name, "_vs_")[[1]]
    if (length(group_names) != 2) {
      cat("Skipping", comparison_name, "- invalid format\n")
      next
    }
    
    group1 <- group_names[1]
    group2 <- group_names[2]
    
    # Subset data for this pairwise comparison
    pair_indices <- groups %in% c(group1, group2)
    X_pair <- X_clr[pair_indices, , drop = FALSE]
    groups_pair <- groups[pair_indices]
    
    cat("Generating plots for:", comparison_name, "\n")
    
    # Generate VIP plot for this pairwise comparison
    tryCatch({
      if (!is.null(pair_result$feature_weights) && length(pair_result$feature_weights) > 0) {
        vip_title <- paste0("SKIOME: VIP (", group1, " vs ", group2, ")")
        vip_plot <- plot_vip(pair_result, top_n = 15, title = vip_title)
        
        vip_filename <- paste0("skiome_vip_", comparison_name, ".ps")
        postscript(vip_filename, width = 10, height = 8, 
                   horizontal = FALSE, onefile = FALSE, paper = "special")
        print(vip_plot)
        dev.off()
        
        if (file.exists(vip_filename)) {
          cat("  ✓ VIP plot saved to:", vip_filename, "\n")
        } else {
          cat("  ✗ VIP plot file was not created\n")
        }
      } else {
        cat("  ⚠️  No feature weights available for VIP plot\n")
      }
    }, error = function(e) {
      cat("  ✗ Error generating VIP plot:", e$message, "\n")
    })
    
    # Generate PCoA plot for this pairwise comparison
    tryCatch({
      if (!is.null(pair_result$distance_matrix)) {
        pcoa_title <- paste0("SKIOME: PCoA (", group1, " vs ", group2, ")")
        pcoa_plot <- plot_pcoa(pair_result, X_pair, groups_pair, title = pcoa_title)
        
        pcoa_filename <- paste0("skiome_pcoa_", comparison_name, ".ps")
        postscript(pcoa_filename, width = 10, height = 8, 
                   horizontal = FALSE, onefile = FALSE, paper = "special")
        print(pcoa_plot)
        dev.off()
        
        if (file.exists(pcoa_filename)) {
          cat("  ✓ PCoA plot saved to:", pcoa_filename, "\n")
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
  }
}

cat("\n")

cat("==============================================================================\n")
cat("Validation Complete!\n")
cat("==============================================================================\n\n")
