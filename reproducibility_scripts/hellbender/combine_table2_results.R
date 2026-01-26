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

# Summary by effect size and sample size - Compare to each traditional method individually
summary_table <- data.frame()
for (effect in unique(all_results$Effect_Size)) {
  for (n_size in unique(all_results$Sample_Size)) {
    subset_results <- all_results[all_results$Effect_Size == effect & 
                                  all_results$Sample_Size == n_size, ]
    
    if (nrow(subset_results) == 0) next
    
    # Calculate MeLSI statistics
    melsi_power <- mean(subset_results$MeLSI_significant)
    melsi_mean_F <- mean(subset_results$MeLSI_F)
    melsi_sd_F <- sd(subset_results$MeLSI_F)
    
    # Compare to each traditional method individually
    traditional_methods <- c("Euclidean", "BrayCurtis", "Jaccard", "WeightedUniFrac", "UnweightedUniFrac")
    
    for (method in traditional_methods) {
      # Check if this method has results (some may be NA for real datasets)
      method_F_col <- paste0(method, "_F")
      method_power_col <- paste0(method, "_significant")
      
      if (!method_F_col %in% colnames(subset_results)) next
      
      # Calculate statistics for this traditional method
      method_F_values <- subset_results[[method_F_col]]
      method_power_values <- subset_results[[method_power_col]]
      
      # Skip if all NA
      if (all(is.na(method_F_values))) next
      
      method_power <- mean(method_power_values, na.rm = TRUE)
      method_mean_F <- mean(method_F_values, na.rm = TRUE)
      method_sd_F <- sd(method_F_values, na.rm = TRUE)
      
      # Calculate power difference (MeLSI - Traditional)
      power_diff <- melsi_power - method_power
      
      summary_table <- rbind(summary_table, data.frame(
        Effect_Size = effect,
        Sample_Size = n_size,
        n_simulations = nrow(subset_results),
        Traditional_Method = method,
        MeLSI_Power = round(melsi_power * 100, 1),
        MeLSI_Mean_F = round(melsi_mean_F, 3),
        MeLSI_SD_F = round(melsi_sd_F, 3),
        Traditional_Power = round(method_power * 100, 1),
        Traditional_Mean_F = round(method_mean_F, 3),
        Traditional_SD_F = round(method_sd_F, 3),
        Power_Difference = round(power_diff * 100, 1),
        F_Difference = round(melsi_mean_F - method_mean_F, 3)
      ))
    }
  }
}

write.csv(summary_table, "table2_power_analysis_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: table2_power_analysis_summary.csv\n")

cat("\nSummary Table:\n")
print(summary_table)

# Generate recovery metrics summary (interpretability validation)
if (any(c("Precision_5", "Precision_10", "Recall_10", "AUC_ROC") %in% colnames(all_results))) {
  cat("\n==============================================================================\n")
  cat("RECOVERY METRICS SUMMARY (Interpretability Validation)\n")
  cat("==============================================================================\n\n")
  
  recovery_summary <- data.frame()
  for (effect in unique(all_results$Effect_Size)) {
    for (n_size in unique(all_results$Sample_Size)) {
      subset_results <- all_results[all_results$Effect_Size == effect & 
                                    all_results$Sample_Size == n_size, ]
      
      if (nrow(subset_results) == 0) next
      
      recovery_summary <- rbind(recovery_summary, data.frame(
        Effect_Size = effect,
        Sample_Size = n_size,
        n_simulations = nrow(subset_results),
        Precision_5 = round(mean(subset_results$Precision_5, na.rm = TRUE), 3),
        Precision_10 = round(mean(subset_results$Precision_10, na.rm = TRUE), 3),
        Precision_20 = round(mean(subset_results$Precision_20, na.rm = TRUE), 3),
        Recall_5 = round(mean(subset_results$Recall_5, na.rm = TRUE), 3),
        Recall_10 = round(mean(subset_results$Recall_10, na.rm = TRUE), 3),
        Recall_20 = round(mean(subset_results$Recall_20, na.rm = TRUE), 3),
        Mean_Rank = round(mean(subset_results$Mean_Rank, na.rm = TRUE), 1),
        AUC_ROC = round(mean(subset_results$AUC_ROC, na.rm = TRUE), 3)
      ))
    }
  }
  
  write.csv(recovery_summary, "table2_recovery_metrics_summary.csv", row.names = FALSE)
  cat("Recovery metrics summary saved to: table2_recovery_metrics_summary.csv\n\n")
  print(recovery_summary)
}
