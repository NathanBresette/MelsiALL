#!/usr/bin/env Rscript
# ==============================================================================
# Combine Table 5 Results from Parallel Simulations
# ==============================================================================
# This script combines all individual simulation result files into a single
# comprehensive results file for analysis.
# ==============================================================================

# Get all simulation result files
result_files <- list.files(pattern = "table5_sim_.*\\.csv", full.names = TRUE)

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
write.csv(all_results, "table5_prefiltering_results.csv", row.names = FALSE)
cat("Combined results saved to: table5_prefiltering_results.csv\n")

# Generate summary statistics
cat("\nGenerating summary statistics...\n")

# Summary by effect size
summary_table <- data.frame()
for (effect in unique(all_results$Effect_Size)) {
  subset_results <- all_results[all_results$Effect_Size == effect, ]
  
  if (nrow(subset_results) == 0) next
  
  summary_table <- rbind(summary_table, data.frame(
    Effect_Size = effect,
    n_samples = subset_results$n_samples[1],
    n_taxa = subset_results$n_taxa[1],
    n_simulations = nrow(subset_results),
    With_Prefilter_Mean_F = round(mean(subset_results$With_Prefilter_F), 3),
    With_Prefilter_SD_F = round(sd(subset_results$With_Prefilter_F), 3),
    With_Prefilter_Power = round(mean(subset_results$With_Prefilter_significant) * 100, 1),
    Without_Prefilter_Mean_F = round(mean(subset_results$Without_Prefilter_F), 3),
    Without_Prefilter_SD_F = round(sd(subset_results$Without_Prefilter_F), 3),
    Without_Prefilter_Power = round(mean(subset_results$Without_Prefilter_significant) * 100, 1),
    Mean_F_Improvement = round(mean(subset_results$F_Improvement), 3),
    Mean_F_Improvement_Pct = round(mean(subset_results$F_Improvement_Pct), 1),
    Mean_Time_Reduction_Pct = round(mean(subset_results$Time_Reduction_Pct), 1)
  ))
}

write.csv(summary_table, "table5_prefiltering_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table5_prefiltering_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)
