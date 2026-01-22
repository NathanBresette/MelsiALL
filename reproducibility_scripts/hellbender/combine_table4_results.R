#!/usr/bin/env Rscript
# ==============================================================================
# Combine Table 4 Results from Parallel Simulations
# ==============================================================================
# This script combines all individual simulation result files into a single
# comprehensive results file for analysis.
# ==============================================================================

# Get all simulation result files
result_files <- list.files(pattern = "table4_sim_.*\\.csv", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No simulation result files found. Make sure the parallel jobs completed successfully.")
}

cat("Found", length(result_files), "simulation result files\n")
cat("Combining results...\n")

# Read and combine all results
all_results <- do.call(rbind, lapply(result_files, function(f) {
  tryCatch({
    read.csv(f, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("Warning: Could not read", f, "\n")
    return(NULL)
  })
}))

# Remove any NULL entries
all_results <- all_results[!sapply(all_results, is.null)]

if (nrow(all_results) == 0) {
  stop("No valid results found in simulation files.")
}

cat("Combined", nrow(all_results), "simulation results\n")

# Save combined results
write.csv(all_results, "table4_parameter_sensitivity_results.csv", row.names = FALSE)
cat("Combined results saved to: table4_parameter_sensitivity_results.csv\n")

# Generate summary statistics
cat("\nGenerating summary statistics...\n")

# Summary by parameter type and value
# Structure: Each replication generates one dataset, tests all parameters on it
summary_table <- data.frame()
for (param_type in unique(all_results$Parameter_Type)) {
  for (param_val in unique(all_results$Parameter_Value)) {
    subset_results <- all_results[all_results$Parameter_Type == param_type & 
                                  all_results$Parameter_Value == param_val, ]
    
    if (nrow(subset_results) == 0) next
    
    summary_table <- rbind(summary_table, data.frame(
      Parameter_Type = param_type,
      Parameter_Value = param_val,
      n_replications = length(unique(subset_results$Replication)),
      Mean_F = round(mean(subset_results$F_statistic), 3),
      SD_F = round(sd(subset_results$F_statistic), 3),
      Mean_p = round(mean(subset_results$p_value), 3),
      SD_p = round(sd(subset_results$p_value), 3),
      Mean_Time = round(mean(subset_results$time_sec), 1),
      SD_Time = round(sd(subset_results$time_sec), 1)
    ))
  }
}

write.csv(summary_table, "table4_parameter_sensitivity_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table4_parameter_sensitivity_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)
