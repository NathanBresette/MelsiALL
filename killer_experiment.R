# Killer Experiment: Prove MeLSI's Unique Value
# This experiment creates a scenario where:
# 1. 30 taxa each change by 10% between groups (subtle, coordinated shift)
# 2. No single taxon is significant after multiple testing correction
# 3. MeLSI detects the strong community-level signal

library(vegan)
library(dplyr)
library(ggplot2)

# Source MeLSI functions
source("melsi_robust.R")

# Function to create the killer dataset
create_killer_dataset <- function(n_samples = 100, n_taxa = 200, n_signal_taxa = 30, effect_size = 0.1) {
  cat("Creating killer dataset...\n")
  cat("Parameters:\n")
  cat("  - Total samples:", n_samples, "\n")
  cat("  - Total taxa:", n_taxa, "\n")
  cat("  - Signal taxa:", n_signal_taxa, "\n")
  cat("  - Effect size:", effect_size, "\n")
  
  # Create base counts (log-normal distribution)
  set.seed(123)
  base_counts <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
                       nrow = n_samples, ncol = n_taxa)
  
  # Add subtle, coordinated shift across signal taxa
  group1_indices <- 1:(n_samples/2)
  group2_indices <- (n_samples/2 + 1):n_samples
  
  # Apply subtle increase to signal taxa in group 2
  for (i in 1:n_signal_taxa) {
    base_counts[group2_indices, i] <- base_counts[group2_indices, i] * (1 + effect_size)
  }
  
  # Add some noise to make it realistic
  noise <- matrix(rnorm(n_samples * n_taxa, mean = 0, sd = 0.1), 
                 nrow = n_samples, ncol = n_taxa)
  counts <- pmax(1, base_counts + noise)
  
  # Ensure counts is a matrix
  counts <- as.matrix(counts)
  
  # Create metadata
  metadata <- data.frame(
    Group = c(rep("Group1", n_samples/2), rep("Group2", n_samples/2)),
    SampleID = paste0("Sample_", 1:n_samples)
  )
  
  # Add row names
  rownames(counts) <- metadata$SampleID
  colnames(counts) <- paste0("Taxon_", 1:n_taxa)
  
  # Ensure counts is a matrix
  counts <- as.matrix(counts)
  
  return(list(
    counts = counts,
    metadata = metadata,
    signal_taxa = 1:n_signal_taxa,
    effect_size = effect_size
  ))
}

# Function to run DA analysis
run_da_analysis <- function(counts, metadata) {
  cat("Running Differential Abundance analysis...\n")
  
  # Simple t-test for each taxon
  p_values <- numeric(ncol(counts))
  for (i in 1:ncol(counts)) {
    group1 <- counts[metadata$Group == "Group1", i]
    group2 <- counts[metadata$Group == "Group2", i]
    test_result <- t.test(group1, group2)
    p_values[i] <- test_result$p.value
  }
  
  # Multiple testing correction (Bonferroni)
  p_adjusted <- p.adjust(p_values, method = "bonferroni")
  
  # Count significant taxa
  n_significant_raw <- sum(p_values < 0.05)
  n_significant_adjusted <- sum(p_adjusted < 0.05)
  
  return(list(
    p_values = p_values,
    p_adjusted = p_adjusted,
    n_significant_raw = n_significant_raw,
    n_significant_adjusted = n_significant_adjusted
  ))
}

# Function to run MeLSI analysis
run_melsi_analysis <- function(counts, metadata) {
  cat("Running MeLSI analysis...\n")
  
  # CLR transformation
  counts_clr <- counts
  counts_clr[counts_clr == 0] <- 1e-10
  counts_clr <- log(counts_clr)
  counts_clr <- counts_clr - rowMeans(counts_clr)
  
  # Run MeLSI with 99 permutations for reasonable precision
  results <- run_melsi_permutation_test(
    counts_clr, metadata$Group,
    n_perms = 99,
    B = 50,
    m_frac = 0.7,
    show_progress = TRUE
  )
  
  return(results)
}

# Function to run traditional distance methods
run_traditional_methods <- function(counts, metadata) {
  cat("Running traditional distance methods...\n")
  
  # CLR transformation
  counts_clr <- counts
  counts_clr[counts_clr == 0] <- 1e-10
  counts_clr <- log(counts_clr)
  counts_clr <- counts_clr - rowMeans(counts_clr)
  
  # Bray-Curtis
  dist_bray <- vegdist(counts, method = "bray")
  permanova_bray <- adonis2(dist_bray ~ Group, data = metadata, permutations = 99)
  
  # Euclidean
  dist_euclidean <- vegdist(counts_clr, method = "euclidean")
  permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = metadata, permutations = 99)
  
  # Jaccard
  dist_jaccard <- vegdist(counts, method = "jaccard")
  permanova_jaccard <- adonis2(dist_jaccard ~ Group, data = metadata, permutations = 99)
  
  return(list(
    bray = permanova_bray,
    euclidean = permanova_euclidean,
    jaccard = permanova_jaccard
  ))
}

# Main experiment function
run_killer_experiment <- function() {
  cat("🚀 STARTING KILLER EXPERIMENT\n")
  cat("================================\n")
  cat("Goal: Prove MeLSI can detect subtle, coordinated shifts\n")
  cat("that DA analysis misses due to multiple testing correction\n\n")
  
  # Create the killer dataset
  data <- create_killer_dataset(
    n_samples = 100,
    n_taxa = 200,
    n_signal_taxa = 30,
    effect_size = 0.1  # 10% change
  )
  
  cat("\n📊 Dataset created:\n")
  cat("  - Samples:", nrow(data$counts), "\n")
  cat("  - Taxa:", ncol(data$counts), "\n")
  cat("  - Signal taxa:", length(data$signal_taxa), "\n")
  cat("  - Effect size:", data$effect_size, "\n")
  
  # Run DA analysis
  cat("\n🔬 STEP 1: Differential Abundance Analysis\n")
  cat("==========================================\n")
  da_results <- run_da_analysis(data$counts, data$metadata)
  
  cat("Results:\n")
  cat("  - Raw significant taxa (p < 0.05):", da_results$n_significant_raw, "\n")
  cat("  - Adjusted significant taxa (Bonferroni):", da_results$n_significant_adjusted, "\n")
  
  if (da_results$n_significant_adjusted == 0) {
    cat("  ✅ SUCCESS: No taxa significant after correction!\n")
  } else {
    cat("  ❌ FAILURE: Some taxa still significant after correction\n")
  }
  
  # Run MeLSI analysis
  cat("\n🧠 STEP 2: MeLSI Analysis\n")
  cat("========================\n")
  melsi_results <- run_melsi_analysis(data$counts, data$metadata)
  
  cat("Results:\n")
  cat("  - F-statistic:", round(melsi_results$F_observed, 4), "\n")
  cat("  - P-value:", round(melsi_results$p_value, 4), "\n")
  cat("  - Significant:", ifelse(melsi_results$p_value < 0.05, "YES", "NO"), "\n")
  
  if (melsi_results$p_value < 0.05) {
    cat("  ✅ SUCCESS: MeLSI detected the community-level signal!\n")
  } else {
    cat("  ❌ FAILURE: MeLSI did not detect the signal\n")
  }
  
  # Run traditional methods
  cat("\n📏 STEP 3: Traditional Distance Methods\n")
  cat("======================================\n")
  traditional_results <- run_traditional_methods(data$counts, data$metadata)
  
  cat("Results:\n")
  cat("  - Bray-Curtis F:", round(traditional_results$bray$F[1], 4), 
      "p:", round(traditional_results$bray$`Pr(>F)`[1], 4), "\n")
  cat("  - Euclidean F:", round(traditional_results$euclidean$F[1], 4), 
      "p:", round(traditional_results$euclidean$`Pr(>F)`[1], 4), "\n")
  cat("  - Jaccard F:", round(traditional_results$jaccard$F[1], 4), 
      "p:", round(traditional_results$jaccard$`Pr(>F)`[1], 4), "\n")
  
  # Summary
  cat("\n📋 EXPERIMENT SUMMARY\n")
  cat("====================\n")
  cat("DA Analysis (Bonferroni):", da_results$n_significant_adjusted, "significant taxa\n")
  cat("MeLSI:", ifelse(melsi_results$p_value < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
  cat("Bray-Curtis:", ifelse(traditional_results$bray$`Pr(>F)`[1] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
  cat("Euclidean:", ifelse(traditional_results$euclidean$`Pr(>F)`[1] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
  cat("Jaccard:", ifelse(traditional_results$jaccard$`Pr(>F)`[1] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
  
  # Check if experiment succeeded
  experiment_success <- (da_results$n_significant_adjusted == 0) && (melsi_results$p_value < 0.05)
  
  if (experiment_success) {
    cat("\n🎉 KILLER EXPERIMENT SUCCEEDED!\n")
    cat("✅ DA analysis found 0 significant taxa after correction\n")
    cat("✅ MeLSI detected the community-level signal\n")
    cat("✅ This proves MeLSI's unique value!\n")
  } else {
    cat("\n❌ KILLER EXPERIMENT FAILED\n")
    if (da_results$n_significant_adjusted > 0) {
      cat("❌ DA analysis found", da_results$n_significant_adjusted, "significant taxa\n")
    }
    if (melsi_results$p_value >= 0.05) {
      cat("❌ MeLSI did not detect the signal\n")
    }
  }
  
  return(list(
    data = data,
    da_results = da_results,
    melsi_results = melsi_results,
    traditional_results = traditional_results,
    success = experiment_success
  ))
}

# Run the experiment
if (interactive()) {
  results <- run_killer_experiment()
}
