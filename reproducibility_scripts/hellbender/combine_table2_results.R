#!/usr/bin/env Rscript
# ==============================================================================
# Combine Table 2 Results from Parallel Simulations
# ==============================================================================
# This script combines all individual simulation result files into a single
# comprehensive results file for analysis.
# ==============================================================================

# Get all simulation result files
result_files <- list.files(pattern = "table2_sim_.*\\.csv", full.names = TRUE)

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
write.csv(all_results, "table2_power_analysis_results.csv", row.names = FALSE)
cat("Combined results saved to: table2_power_analysis_results.csv\n")

# Generate summary statistics
cat("\nGenerating summary statistics...\n")

# Summary by effect size and sample size
summary_table <- data.frame()
for (effect in unique(all_results$Effect_Size)) {
  for (n_size in unique(all_results$Sample_Size)) {
    subset_results <- all_results[all_results$Effect_Size == effect & 
                                  all_results$Sample_Size == n_size, ]
    
    if (nrow(subset_results) == 0) next
    
    # Find best traditional method by mean F-statistic
    trad_means <- c(
      Euclidean = mean(subset_results$Euclidean_F),
      BrayCurtis = mean(subset_results$BrayCurtis_F),
      Jaccard = mean(subset_results$Jaccard_F),
      WeightedUniFrac = mean(subset_results$WeightedUniFrac_F, na.rm = TRUE),
      UnweightedUniFrac = mean(subset_results$UnweightedUniFrac_F, na.rm = TRUE)
    )
    best_trad_method <- names(which.max(trad_means))
    
    if (best_trad_method == "Euclidean") {
      best_trad_power <- mean(subset_results$Euclidean_significant)
      best_trad_mean_F <- mean(subset_results$Euclidean_F)
    } else if (best_trad_method == "BrayCurtis") {
      best_trad_power <- mean(subset_results$BrayCurtis_significant)
      best_trad_mean_F <- mean(subset_results$BrayCurtis_F)
    } else if (best_trad_method == "Jaccard") {
      best_trad_power <- mean(subset_results$Jaccard_significant)
      best_trad_mean_F <- mean(subset_results$Jaccard_F)
    } else if (best_trad_method == "WeightedUniFrac") {
      best_trad_power <- mean(subset_results$WeightedUniFrac_significant)
      best_trad_mean_F <- mean(subset_results$WeightedUniFrac_F)
    } else {
      best_trad_power <- mean(subset_results$UnweightedUniFrac_significant)
      best_trad_mean_F <- mean(subset_results$UnweightedUniFrac_F)
    }
    
    summary_table <- rbind(summary_table, data.frame(
      Effect_Size = effect,
      Sample_Size = n_size,
      n_simulations = nrow(subset_results),
      MeLSI_Power = round(mean(subset_results$MeLSI_significant) * 100, 1),
      MeLSI_Mean_F = round(mean(subset_results$MeLSI_F), 3),
      Best_Traditional_Power = round(best_trad_power * 100, 1),
      Best_Traditional_Mean_F = round(best_trad_mean_F, 3)
    ))
  }
}

write.csv(summary_table, "table2_power_analysis_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table2_power_analysis_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)
