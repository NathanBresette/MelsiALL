#!/usr/bin/env Rscript
# ==============================================================================
# Table 4: Parameter Sensitivity Analysis
# ==============================================================================
# Rigorous version with repeated simulations for variance estimation
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
n_replications <- 25  # Rigorous: 25 replications for variance estimation

# Parameter values to test
B_values <- c(1, 10, 20, 30, 50, 100)
m_frac_values <- c(0.5, 0.7, 0.8, 0.9, 1.0)

# Total simulations: n_replications × (n_B_values + n_m_frac_values) = 25 × (6 + 5) = 275
total_simulations <- n_replications * (length(B_values) + length(m_frac_values))

# ==============================================================================
# Helper Function: Generate Test Dataset
# ==============================================================================
generate_test_dataset <- function(seed) {
  set.seed(seed)  # Set seed for reproducibility
  n_samples <- 100
  n_taxa <- 200
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                  nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  # Add signal (medium effect)
  group1_size <- n_samples %/% 2
  signal_taxa <- sample(1:n_taxa, 10)
  for (i in 1:group1_size) {
    counts[i, signal_taxa] <- counts[i, signal_taxa] * 2
  }
  
  group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))
  
  # CLR transformation
  X_clr <- log(counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  return(list(X_clr = X_clr, group = group))
}

# ==============================================================================
# Run Parameter Sensitivity Tests - REPEATED SIMULATIONS
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 4: Parameter Sensitivity Analysis (Rigorous)\n")
cat("Running multiple simulations per parameter value for variance estimation\n")
cat("==============================================================================\n\n")

results <- data.frame()

# Determine which simulation to run (for parallel execution)
if (parallel_mode) {
  # Structure: Each replication generates ONE dataset, then tests all parameters on it
  # Replication 1: dataset seed = 42 + 1, test all 11 parameter combinations
  # Replication 2: dataset seed = 42 + 2, test all 11 parameter combinations
  # ...
  # Replication 25: dataset seed = 42 + 25, test all 11 parameter combinations
  # Total: 25 replications × 11 parameter combinations = 275 simulations
  
  # Calculate which replication and parameter this simulation belongs to
  n_params_per_replication <- length(B_values) + length(m_frac_values)  # 11
  
  replication <- ((sim_index - 1) %/% n_params_per_replication) + 1
  param_index_within_replication <- ((sim_index - 1) %% n_params_per_replication) + 1
  
  if (replication < 1 || replication > n_replications) {
    stop("replication out of range")
  }
  
  # Determine which parameter to test
  if (param_index_within_replication <= length(B_values)) {
    param_type <- "B"
    param_val <- B_values[param_index_within_replication]
  } else {
    param_type <- "m_frac"
    param_val <- m_frac_values[param_index_within_replication - length(B_values)]
  }
  
  cat(sprintf("Running simulation %d/%d: Replication %d, %s = %s\n", 
              sim_index, total_sims, replication, param_type, param_val))
  
  # Generate dataset for this replication (same dataset for all parameters in this replication)
  dataset_seed <- 42 + replication
  data <- generate_test_dataset(dataset_seed)
  
  # Run MeLSI with appropriate parameters
  if (param_type == "B") {
    start_time <- Sys.time()
    melsi_result <- melsi(data$X_clr, data$group, n_perms = 200, B = param_val, 
                         m_frac = 0.8, show_progress = FALSE, plot_vip = FALSE)
    time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    result <- data.frame(
      Replication = replication,
      Parameter_Type = "B",
      Parameter_Value = param_val,
      F_statistic = melsi_result$F_observed,
      p_value = melsi_result$p_value,
      time_sec = time_taken
    )
  } else {
    start_time <- Sys.time()
    melsi_result <- melsi(data$X_clr, data$group, n_perms = 200, B = 30, 
                         m_frac = param_val, show_progress = FALSE, plot_vip = FALSE)
    time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    result <- data.frame(
      Replication = replication,
      Parameter_Type = "m_frac",
      Parameter_Value = param_val,
      F_statistic = melsi_result$F_observed,
      p_value = melsi_result$p_value,
      time_sec = time_taken
    )
  }
  
  # Save single result to file (for parallel collection)
  output_file <- sprintf("table4_sim_%d.csv", sim_index)
  write.csv(result, output_file, row.names = FALSE)
  cat("Result saved to:", output_file, "\n")
  
  # Exit early in parallel mode - don't calculate summary statistics
  quit(status = 0)
  
} else {
  # Sequential mode: run all simulations
  # Structure matches original: Generate one dataset per replication, test all parameters on it
  for (rep in 1:n_replications) {
    cat(sprintf("\nReplication %d/%d: Generating dataset and testing all parameters...\n", 
                rep, n_replications))
    
    # Generate ONE dataset for this replication (matches original approach)
    dataset_seed <- 42 + rep
    data <- generate_test_dataset(dataset_seed)
    
    # Test Ensemble Size (B) on this dataset
    cat("  Testing B values...\n")
    for (B in B_values) {
      start_time <- Sys.time()
      melsi_result <- melsi(data$X_clr, data$group, n_perms = 200, B = B, m_frac = 0.8,
                           show_progress = FALSE, plot_vip = FALSE)
      time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      results <- rbind(results, data.frame(
        Replication = rep,
        Parameter_Type = "B",
        Parameter_Value = B,
        F_statistic = melsi_result$F_observed,
        p_value = melsi_result$p_value,
        time_sec = time_taken
      ))
    }
    
    # Test Feature Fraction (m_frac) on this same dataset
    cat("  Testing m_frac values...\n")
    for (m_frac in m_frac_values) {
      start_time <- Sys.time()
      melsi_result <- melsi(data$X_clr, data$group, n_perms = 200, B = 30, m_frac = m_frac,
                           show_progress = FALSE, plot_vip = FALSE)
      time_taken <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      results <- rbind(results, data.frame(
        Replication = rep,
        Parameter_Type = "m_frac",
        Parameter_Value = m_frac,
        F_statistic = melsi_result$F_observed,
        p_value = melsi_result$p_value,
        time_sec = time_taken
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
  
  # Summary by parameter type (matches original structure)
  cat("\nPart A: Ensemble Size (B) Sensitivity\n")
  cat("------------------------------------------------------------------------------\n")
  for (B in B_values) {
    B_subset <- results[results$Parameter_Type == "B" & results$Parameter_Value == B, ]
    if (nrow(B_subset) > 0) {
      cat(sprintf("B = %d: Mean F = %.3f (SD = %.3f), Mean Time = %.1fs (SD = %.1fs)\n",
                  B, mean(B_subset$F_statistic), sd(B_subset$F_statistic),
                  mean(B_subset$time_sec), sd(B_subset$time_sec)))
    }
  }
  
  cat("\nPart B: Feature Fraction (m_frac) Sensitivity\n")
  cat("------------------------------------------------------------------------------\n")
  for (m_frac in m_frac_values) {
    m_subset <- results[results$Parameter_Type == "m_frac" & results$Parameter_Value == m_frac, ]
    if (nrow(m_subset) > 0) {
      cat(sprintf("m_frac = %.1f: Mean F = %.3f (SD = %.3f), Mean Time = %.1fs (SD = %.1fs)\n",
                  m_frac, mean(m_subset$F_statistic), sd(m_subset$F_statistic),
                  mean(m_subset$time_sec), sd(m_subset$time_sec)))
    }
  }
  
  # Save detailed results
  write.csv(results, "table4_parameter_sensitivity_results.csv", row.names = FALSE)
  cat("\nDetailed results saved to: table4_parameter_sensitivity_results.csv\n")
  cat("==============================================================================\n")
}
