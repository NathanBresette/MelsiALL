#!/usr/bin/env Rscript
# ==============================================================================
# Combine Table 3 Results from Parallel Simulations
# ==============================================================================
# This script combines all individual simulation result files into a single
# comprehensive results file for analysis.
# ==============================================================================

# Get all simulation result files
result_files <- list.files(pattern = "table3_sim_.*\\.csv", full.names = TRUE)

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
write.csv(all_results, "table3_scalability_results.csv", row.names = FALSE)
cat("Combined results saved to: table3_scalability_results.csv\n")

# Generate summary statistics
cat("\nGenerating summary statistics...\n")

# Summary by condition type and value
summary_table <- data.frame()
for (cond_type in unique(all_results$Condition_Type)) {
  for (cond_val in unique(all_results$Condition_Value)) {
    subset_results <- all_results[all_results$Condition_Type == cond_type & 
                                  all_results$Condition_Value == cond_val, ]
    
    if (nrow(subset_results) == 0) next
    
    summary_table <- rbind(summary_table, data.frame(
      Condition_Type = cond_type,
      Condition_Value = cond_val,
      n_samples = subset_results$n_samples[1],
      n_taxa = subset_results$n_taxa[1],
      n_simulations = nrow(subset_results),
      MeLSI_Mean_F = round(mean(subset_results$MeLSI_F), 3),
      MeLSI_SD_F = round(sd(subset_results$MeLSI_F), 3),
      MeLSI_Mean_Time = round(mean(subset_results$MeLSI_time_sec), 1),
      MeLSI_SD_Time = round(sd(subset_results$MeLSI_time_sec), 1),
      Best_Traditional_Mean_F = round(mean(subset_results$Best_Traditional_F), 3),
      Best_Traditional_SD_F = round(sd(subset_results$Best_Traditional_F), 3),
      Best_Traditional_Mean_Time = round(mean(subset_results$Best_Traditional_time_sec), 1),
      Best_Traditional_SD_Time = round(sd(subset_results$Best_Traditional_time_sec), 1)
    ))
  }
}

write.csv(summary_table, "table3_scalability_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table3_scalability_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)
