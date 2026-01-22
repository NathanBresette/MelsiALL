#!/usr/bin/env Rscript
# ==============================================================================
# Combine Table 1 Results from Parallel Simulations
# ==============================================================================
# This script combines all individual simulation result files into a single
# comprehensive results file for analysis.
# ==============================================================================

# Get all simulation result files
result_files <- list.files(pattern = "table1_sim_.*\\.csv", full.names = TRUE)

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
write.csv(all_results, "table1_type1_error_results.csv", row.names = FALSE)
cat("Combined results saved to: table1_type1_error_results.csv\n")

# Generate summary statistics
cat("\nGenerating summary statistics...\n")

# Summary by dataset type and sample size
summary_table <- data.frame()
for (dataset_type in unique(all_results$Dataset_Type)) {
  for (n_size in unique(all_results$Sample_Size)) {
    subset_results <- all_results[all_results$Dataset_Type == dataset_type & 
                                  all_results$Sample_Size == n_size, ]
    
    if (nrow(subset_results) == 0) next
    
    summary_table <- rbind(summary_table, data.frame(
      Dataset_Type = dataset_type,
      Sample_Size = n_size,
      n_simulations = nrow(subset_results),
      MeLSI_TypeI_Rate = round(mean(subset_results$MeLSI_significant) * 100, 2),
      Euclidean_TypeI_Rate = round(mean(subset_results$Euclidean_significant) * 100, 2),
      BrayCurtis_TypeI_Rate = round(mean(subset_results$BrayCurtis_significant) * 100, 2)
    ))
  }
}

write.csv(summary_table, "table1_type1_error_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table1_type1_error_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)
