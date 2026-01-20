# Null Hypothesis Test: Prove MeLSI has proper Type I error control
# This test creates datasets with NO real differences between groups
# and verifies that MeLSI (and other methods) correctly find no significance

library(vegan)
source("melsi_robust.R")

# Function to create null datasets (no real group differences)
create_null_dataset <- function(n_samples = 100, n_taxa = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cat("Creating null dataset (no real group differences)...\n")
  cat("Parameters:\n")
  cat("  - Samples:", n_samples, "\n")
  cat("  - Taxa:", n_taxa, "\n")
  cat("  - Group differences: NONE (null hypothesis)\n")
  
  # Create counts from the same distribution for both groups
  # (no systematic differences between groups)
  counts <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
                  nrow = n_samples, ncol = n_taxa)
  
  # Create metadata with random group assignments
  # (groups are assigned randomly, so no real biological difference)
  metadata <- data.frame(
    Group = sample(c("Group1", "Group2"), n_samples, replace = TRUE)
  )
  
  return(list(
    counts = counts,
    metadata = metadata,
    description = "Null dataset with no real group differences"
  ))
}

# Function to run comprehensive null hypothesis testing
run_null_hypothesis_tests <- function(n_datasets = 10, n_samples = 100, n_taxa = 200) {
  cat("🚀 NULL HYPOTHESIS TESTING\n")
  cat("===========================\n")
  cat("Testing Type I error control across", n_datasets, "null datasets\n")
  cat("Goal: Verify MeLSI doesn't have inflated false positive rates\n\n")
  
  results <- data.frame(
    dataset = character(),
    melsi_significant = logical(),
    melsi_p_value = numeric(),
    bray_significant = logical(),
    bray_p_value = numeric(),
    euclidean_significant = logical(),
    euclidean_p_value = numeric(),
    jaccard_significant = logical(),
    jaccard_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_datasets) {
    cat("\n--- Dataset", i, "of", n_datasets, "---\n")
    
    # Create null dataset
    data <- create_null_dataset(n_samples, n_taxa, seed = i)
    
    # Test MeLSI
    cat("Testing MeLSI...\n")
    counts_clr <- log(data$counts + 1)
    counts_clr <- counts_clr - rowMeans(counts_clr)
    
    melsi_results <- run_melsi_permutation_test(
      counts_clr, data$metadata$Group,
      n_perms = 99,  # Reasonable precision for null testing
      B = 50,
      m_frac = 0.7,
      show_progress = FALSE
    )
    
    # Test traditional methods
    cat("Testing traditional methods...\n")
    
    # Bray-Curtis
    dist_bray <- vegdist(data$counts, method = "bray")
    permanova_bray <- adonis2(dist_bray ~ Group, data = data$metadata, permutations = 99)
    
    # Euclidean
    dist_euclidean <- vegdist(counts_clr, method = "euclidean")
    permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = data$metadata, permutations = 99)
    
    # Jaccard
    dist_jaccard <- vegdist(data$counts, method = "jaccard")
    permanova_jaccard <- adonis2(dist_jaccard ~ Group, data = data$metadata, permutations = 99)
    
    # Store results
    results <- rbind(results, data.frame(
      dataset = paste0("Null_", i),
      melsi_significant = melsi_results$p_value < 0.05,
      melsi_p_value = melsi_results$p_value,
      bray_significant = permanova_bray$`Pr(>F)`[1] < 0.05,
      bray_p_value = permanova_bray$`Pr(>F)`[1],
      euclidean_significant = permanova_euclidean$`Pr(>F)`[1] < 0.05,
      euclidean_p_value = permanova_euclidean$`Pr(>F)`[1],
      jaccard_significant = permanova_jaccard$`Pr(>F)`[1] < 0.05,
      jaccard_p_value = permanova_jaccard$`Pr(>F)`[1]
    ))
    
    # Show current results
    cat("Results:\n")
    cat("  MeLSI: p =", round(melsi_results$p_value, 4), 
        ifelse(melsi_results$p_value < 0.05, "***", ""), "\n")
    cat("  Bray-Curtis: p =", round(permanova_bray$`Pr(>F)`[1], 4), 
        ifelse(permanova_bray$`Pr(>F)`[1] < 0.05, "***", ""), "\n")
    cat("  Euclidean: p =", round(permanova_euclidean$`Pr(>F)`[1], 4), 
        ifelse(permanova_euclidean$`Pr(>F)`[1] < 0.05, "***", ""), "\n")
    cat("  Jaccard: p =", round(permanova_jaccard$`Pr(>F)`[1], 4), 
        ifelse(permanova_jaccard$`Pr(>F)`[1] < 0.05, "***", ""), "\n")
  }
  
  # Summary statistics
  cat("\n📊 TYPE I ERROR ANALYSIS\n")
  cat("========================\n")
  
  melsi_false_positives <- sum(results$melsi_significant)
  bray_false_positives <- sum(results$bray_significant)
  euclidean_false_positives <- sum(results$euclidean_significant)
  jaccard_false_positives <- sum(results$jaccard_significant)
  
  melsi_type1_rate <- melsi_false_positives / n_datasets
  bray_type1_rate <- bray_false_positives / n_datasets
  euclidean_type1_rate <- euclidean_false_positives / n_datasets
  jaccard_type1_rate <- jaccard_false_positives / n_datasets
  
  cat("False Positive Rates (should be ~5%):\n")
  cat("  MeLSI:", sprintf("%.1f%% (%d/%d)", melsi_type1_rate*100, melsi_false_positives, n_datasets), "\n")
  cat("  Bray-Curtis:", sprintf("%.1f%% (%d/%d)", bray_type1_rate*100, bray_false_positives, n_datasets), "\n")
  cat("  Euclidean:", sprintf("%.1f%% (%d/%d)", euclidean_type1_rate*100, euclidean_false_positives, n_datasets), "\n")
  cat("  Jaccard:", sprintf("%.1f%% (%d/%d)", jaccard_type1_rate*100, jaccard_false_positives, n_datasets), "\n")
  
  # Statistical test for Type I error control
  cat("\n📈 P-VALUE DISTRIBUTION\n")
  cat("=======================\n")
  
  # Check if p-values are uniformly distributed (should be for null data)
  melsi_p_values <- results$melsi_p_value
  bray_p_values <- results$bray_p_value
  euclidean_p_values <- results$euclidean_p_value
  jaccard_p_values <- results$jaccard_p_value
  
  cat("P-value statistics:\n")
  cat("  MeLSI: mean =", round(mean(melsi_p_values), 3), 
      "median =", round(median(melsi_p_values), 3), "\n")
  cat("  Bray-Curtis: mean =", round(mean(bray_p_values), 3), 
      "median =", round(median(bray_p_values), 3), "\n")
  cat("  Euclidean: mean =", round(mean(euclidean_p_values), 3), 
      "median =", round(median(euclidean_p_values), 3), "\n")
  cat("  Jaccard: mean =", round(mean(jaccard_p_values), 3), 
      "median =", round(median(jaccard_p_values), 3), "\n")
  
  # Assessment
  cat("\n🎯 TYPE I ERROR ASSESSMENT\n")
  cat("==========================\n")
  
  # Check if MeLSI has acceptable Type I error rate (5-10% is acceptable)
  if (melsi_type1_rate <= 0.10) {
    cat("✅ MeLSI Type I error rate is ACCEPTABLE (", sprintf("%.1f%%", melsi_type1_rate*100), ")\n")
  } else {
    cat("❌ MeLSI Type I error rate is TOO HIGH (", sprintf("%.1f%%", melsi_type1_rate*100), ")\n")
  }
  
  # Compare to other methods
  if (abs(melsi_type1_rate - bray_type1_rate) <= 0.05) {
    cat("✅ MeLSI Type I error rate is SIMILAR to Bray-Curtis\n")
  } else {
    cat("⚠️  MeLSI Type I error rate DIFFERS from Bray-Curtis\n")
  }
  
  # Save results
  write.csv(results, "null_hypothesis_results.csv", row.names = FALSE)
  cat("\n📁 Results saved to: null_hypothesis_results.csv\n")
  
  return(list(
    results = results,
    type1_rates = c(
      MeLSI = melsi_type1_rate,
      Bray_Curtis = bray_type1_rate,
      Euclidean = euclidean_type1_rate,
      Jaccard = jaccard_type1_rate
    ),
    assessment = melsi_type1_rate <= 0.10
  ))
}

# Run the null hypothesis tests
if (interactive()) {
  null_results <- run_null_hypothesis_tests(n_datasets = 20, n_samples = 100, n_taxa = 200)
}
