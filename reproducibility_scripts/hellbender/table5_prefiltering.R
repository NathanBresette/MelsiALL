#!/usr/bin/env Rscript
# ==============================================================================
# Table 5: Benefit of Conservative Pre-filtering
# ==============================================================================
# Rigorous version with repeated simulations for proper power estimation
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)

# ==============================================================================
# Configuration: Support parallel execution
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  sim_index <- as.integer(args[1])
  total_sims <- as.integer(args[2])
  parallel_mode <- TRUE
} else {
  sim_index <- NULL
  total_sims <- NULL
  parallel_mode <- FALSE
}

# Simulation parameters
n_simulations_per_condition <- 50  # Rigorous: 50 simulations for proper power estimation

# Test cases: effect sizes with different dimensionalities
test_cases <- list(
  list(n = 100, p = 500, effect = "small"),
  list(n = 100, p = 200, effect = "medium"),
  list(n = 100, p = 100, effect = "large")
)

# ==============================================================================
# Helper Function: Generate Data with Varying Sparsity
# ==============================================================================
generate_sparse_data <- function(n_samples = 100, n_taxa = 500, effect_size = "medium") {
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                  nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  # Add realistic sparsity (many zero-inflated features)
  for (i in 1:n_samples) {
    zero_taxa <- sample(1:n_taxa, n_taxa * 0.7)  # 70% sparse
    counts[i, zero_taxa] <- 0
  }
  
  # Add signal
  group1_size <- n_samples %/% 2
  
  if (effect_size == "small") {
    signal_taxa <- 5
    fold_change <- 1.5
  } else if (effect_size == "medium") {
    signal_taxa <- 10
    fold_change <- 2.0
  } else if (effect_size == "large") {
    signal_taxa <- 20
    fold_change <- 3.0
  }
  
  signal_indices <- sample(1:n_taxa, signal_taxa)
  for (i in 1:group1_size) {
    counts[i, signal_indices] <- counts[i, signal_indices] * fold_change
  }
  
  group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))
  
  return(list(counts = counts, group = group))
}

# ==============================================================================
# Run Pre-filtering Tests - REPEATED SIMULATIONS
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 5: Benefit of Conservative Pre-filtering (Rigorous)\n")
cat("Running multiple simulations per effect size for proper power estimation\n")
cat("==============================================================================\n\n")

results <- data.frame()

# Determine which simulation to run (for parallel execution)
if (parallel_mode) {
  # Calculate which condition this simulation belongs to
  conditions <- expand.grid(
    test_case = 1:length(test_cases),
    sim = 1:n_simulations_per_condition
  )
  total_conditions <- nrow(conditions)
  
  if (sim_index < 1 || sim_index > total_conditions) {
    stop("sim_index out of range")
  }
  
  current_condition <- conditions[sim_index, ]
  tc_idx <- as.integer(current_condition$test_case)
  sim <- as.integer(current_condition$sim)
  tc <- test_cases[[tc_idx]]
  
  cat(sprintf("Running simulation %d/%d: %s effect, n=%d, p=%d\n", 
              sim_index, total_conditions, tc$effect, tc$n, tc$p))
  
  # Set unique seed for this simulation
  set.seed(42 + sim_index)
  
  # Generate data
  data <- generate_sparse_data(tc$n, tc$p, tc$effect)
  
  # CLR transformation
  X_clr <- log(data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  # Count features before filtering
  n_features_before <- ncol(X_clr)
  
  # Apply conservative variance-based pre-filtering (retain top 70% by importance score)
  # Importance score: I_j = |μ_1j - μ_2j| / sqrt(σ_1j² + σ_2j²)
  group1_idx <- which(data$group == unique(data$group)[1])
  group2_idx <- which(data$group == unique(data$group)[2])
  
  # Calculate importance scores for each feature
  importance_scores <- numeric(ncol(X_clr))
  for (j in 1:ncol(X_clr)) {
    mu1 <- mean(X_clr[group1_idx, j])
    mu2 <- mean(X_clr[group2_idx, j])
    sigma1_sq <- var(X_clr[group1_idx, j])
    sigma2_sq <- var(X_clr[group2_idx, j])
    denom <- sqrt(sigma1_sq + sigma2_sq)
    # Handle division by zero (constant features)
    if (denom < 1e-10) {
      importance_scores[j] <- 0
    } else {
      importance_scores[j] <- abs(mu1 - mu2) / denom
    }
  }
  
  # Retain top 70% of features by importance score
  n_retain <- round(0.7 * ncol(X_clr))
  top_features <- order(importance_scores, decreasing = TRUE)[1:n_retain]
  X_clr_filtered <- X_clr[, top_features]
  n_features_after <- length(top_features)
  
  # Test WITH pre-filtering
  start_time <- Sys.time()
  melsi_with <- melsi(X_clr_filtered, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                     show_progress = FALSE, plot_vip = FALSE)
  time_with <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Test WITHOUT pre-filtering
  start_time <- Sys.time()
  melsi_without <- melsi(X_clr, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                        show_progress = FALSE, plot_vip = FALSE)
  time_without <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Calculate improvements
  F_improvement <- melsi_with$F_observed - melsi_without$F_observed
  F_improvement_pct <- 100 * F_improvement / melsi_without$F_observed
  time_reduction_pct <- 100 * (time_without - time_with) / time_without
  
  # Store single result
  result <- data.frame(
    Effect_Size = tc$effect,
    n_samples = tc$n,
    n_taxa = tc$p,
    Simulation = sim,
    Features_Before = n_features_before,
    Features_After = n_features_after,
    With_Prefilter_F = melsi_with$F_observed,
    With_Prefilter_p = melsi_with$p_value,
    With_Prefilter_significant = melsi_with$p_value < 0.05,
    Without_Prefilter_F = melsi_without$F_observed,
    Without_Prefilter_p = melsi_without$p_value,
    Without_Prefilter_significant = melsi_without$p_value < 0.05,
    F_Improvement = F_improvement,
    F_Improvement_Pct = F_improvement_pct,
    Time_With = time_with,
    Time_Without = time_without,
    Time_Reduction_Pct = time_reduction_pct
  )
  
  # Save single result to file (for parallel collection)
  output_file <- sprintf("table5_sim_%d.csv", sim_index)
  write.csv(result, output_file, row.names = FALSE)
  cat("Result saved to:", output_file, "\n")
  
  # Exit early in parallel mode - don't calculate summary statistics
  quit(status = 0)
  
} else {
  # Sequential mode: run all simulations
  for (i in 1:length(test_cases)) {
    tc <- test_cases[[i]]
    cat(sprintf("\nTest Case %d: %s effect, n=%d, p=%d (%d simulations)...\n", 
                i, tc$effect, tc$n, tc$p, n_simulations_per_condition))
    
    for (sim in 1:n_simulations_per_condition) {
      if (sim %% 10 == 0) cat("  Simulation", sim, "of", n_simulations_per_condition, "\n")
      
      set.seed(42 + sim + (i - 1) * 100)
      
      # Generate data
      data <- generate_sparse_data(tc$n, tc$p, tc$effect)
      
      # CLR transformation
      X_clr <- log(data$counts + 1)
      X_clr <- X_clr - rowMeans(X_clr)
      colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
      
      # Count features before filtering
      n_features_before <- ncol(X_clr)
      
      # Apply conservative variance-based pre-filtering (retain top 70% by importance score)
      # Importance score: I_j = |μ_1j - μ_2j| / sqrt(σ_1j² + σ_2j²)
      group1_idx <- which(data$group == unique(data$group)[1])
      group2_idx <- which(data$group == unique(data$group)[2])
      
      # Calculate importance scores for each feature
      importance_scores <- numeric(ncol(X_clr))
      for (j in 1:ncol(X_clr)) {
        mu1 <- mean(X_clr[group1_idx, j])
        mu2 <- mean(X_clr[group2_idx, j])
        sigma1_sq <- var(X_clr[group1_idx, j])
        sigma2_sq <- var(X_clr[group2_idx, j])
        denom <- sqrt(sigma1_sq + sigma2_sq)
        # Handle division by zero (constant features)
        if (denom < 1e-10) {
          importance_scores[j] <- 0
        } else {
          importance_scores[j] <- abs(mu1 - mu2) / denom
        }
      }
      
      # Retain top 70% of features by importance score
      n_retain <- round(0.7 * ncol(X_clr))
      top_features <- order(importance_scores, decreasing = TRUE)[1:n_retain]
      X_clr_filtered <- X_clr[, top_features]
      n_features_after <- length(top_features)
      
      # Test WITH pre-filtering
      start_time <- Sys.time()
      melsi_with <- melsi(X_clr_filtered, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                         show_progress = FALSE, plot_vip = FALSE)
      time_with <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Test WITHOUT pre-filtering
      start_time <- Sys.time()
      melsi_without <- melsi(X_clr, data$group, n_perms = 200, B = 30, m_frac = 0.8,
                            show_progress = FALSE, plot_vip = FALSE)
      time_without <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Calculate improvements
      F_improvement <- melsi_with$F_observed - melsi_without$F_observed
      F_improvement_pct <- 100 * F_improvement / melsi_without$F_observed
      time_reduction_pct <- 100 * (time_without - time_with) / time_without
      
      results <- rbind(results, data.frame(
        Effect_Size = tc$effect,
        n_samples = tc$n,
        n_taxa = tc$p,
        Simulation = sim,
        Features_Before = n_features_before,
        Features_After = n_features_after,
        With_Prefilter_F = melsi_with$F_observed,
        With_Prefilter_p = melsi_with$p_value,
        With_Prefilter_significant = melsi_with$p_value < 0.05,
        Without_Prefilter_F = melsi_without$F_observed,
        Without_Prefilter_p = melsi_without$p_value,
        Without_Prefilter_significant = melsi_without$p_value < 0.05,
        F_Improvement = F_improvement,
        F_Improvement_Pct = F_improvement_pct,
        Time_With = time_with,
        Time_Without = time_without,
        Time_Reduction_Pct = time_reduction_pct
      ))
    }
  }
}

# ==============================================================================
# Summary Statistics (only in sequential mode)
# ==============================================================================
if (!parallel_mode && nrow(results) > 0) {
  cat("\n")
  cat("==============================================================================\n")
  cat("SUMMARY STATISTICS\n")
  cat("==============================================================================\n\n")
  
  # Summary by effect size
  for (effect in unique(results$Effect_Size)) {
    subset_results <- results[results$Effect_Size == effect, ]
    
    if (nrow(subset_results) == 0) next
    
    cat(sprintf("\n%s Effect:\n", toupper(substring(effect, 1, 1))))
    cat(sprintf("  With Pre-filtering:\n"))
    cat(sprintf("    Mean F = %.3f (SD = %.3f)\n",
                mean(subset_results$With_Prefilter_F), sd(subset_results$With_Prefilter_F)))
    cat(sprintf("    Power = %.1f%% (%d/%d)\n",
                mean(subset_results$With_Prefilter_significant) * 100,
                sum(subset_results$With_Prefilter_significant), nrow(subset_results)))
    cat(sprintf("  Without Pre-filtering:\n"))
    cat(sprintf("    Mean F = %.3f (SD = %.3f)\n",
                mean(subset_results$Without_Prefilter_F), sd(subset_results$Without_Prefilter_F)))
    cat(sprintf("    Power = %.1f%% (%d/%d)\n",
                mean(subset_results$Without_Prefilter_significant) * 100,
                sum(subset_results$Without_Prefilter_significant), nrow(subset_results)))
    cat(sprintf("  Improvement:\n"))
    cat(sprintf("    Mean F improvement = %.3f (%.1f%%)\n",
                mean(subset_results$F_Improvement), mean(subset_results$F_Improvement_Pct)))
    cat(sprintf("    Mean time reduction = %.1f%%\n",
                mean(subset_results$Time_Reduction_Pct)))
  }
  
  # Save detailed results
  write.csv(results, "table5_prefiltering_results.csv", row.names = FALSE)
  cat("\nDetailed results saved to: table5_prefiltering_results.csv\n")
  cat("==============================================================================\n")
}
