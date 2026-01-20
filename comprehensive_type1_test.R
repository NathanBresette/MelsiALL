# Comprehensive Type I Error Testing for Improved MeLSI
# Tests multiple scenarios to ensure proper false positive control

library(vegan)
source("melsi_improved.R")

# Function to create different types of null datasets
create_null_datasets <- function(n_datasets_per_type = 10) {
  datasets <- list()
  
  # Type 1: Random groups, log-normal data
  cat("Creating Type 1 null datasets (random groups, log-normal)...\n")
  for (i in 1:n_datasets_per_type) {
    set.seed(i)
    counts <- matrix(rlnorm(100 * 200, meanlog = 2, sdlog = 1), nrow = 100, ncol = 200)
    metadata <- data.frame(Group = sample(c("Group1", "Group2"), 100, replace = TRUE))
    datasets[[paste0("Type1_", i)]] <- list(counts = counts, metadata = metadata, type = "random_groups")
  }
  
  # Type 2: Balanced groups, log-normal data  
  cat("Creating Type 2 null datasets (balanced groups, log-normal)...\n")
  for (i in 1:n_datasets_per_type) {
    set.seed(i + 100)
    counts <- matrix(rlnorm(100 * 200, meanlog = 2, sdlog = 1), nrow = 100, ncol = 200)
    metadata <- data.frame(Group = rep(c("Group1", "Group2"), each = 50))
    # Shuffle to break any potential structure
    counts <- counts[sample(1:100), ]
    datasets[[paste0("Type2_", i)]] <- list(counts = counts, metadata = metadata, type = "balanced_groups")
  }
  
  # Type 3: Random groups, normal data
  cat("Creating Type 3 null datasets (random groups, normal)...\n")
  for (i in 1:n_datasets_per_type) {
    set.seed(i + 200)
    counts <- matrix(rnorm(100 * 200, mean = 10, sd = 3), nrow = 100, ncol = 200)
    counts <- pmax(1, counts)  # Ensure positive values
    metadata <- data.frame(Group = sample(c("Group1", "Group2"), 100, replace = TRUE))
    datasets[[paste0("Type3_", i)]] <- list(counts = counts, metadata = metadata, type = "normal_data")
  }
  
  # Type 4: Balanced groups, sparse data (many zeros)
  cat("Creating Type 4 null datasets (balanced groups, sparse)...\n")
  for (i in 1:n_datasets_per_type) {
    set.seed(i + 300)
    # Create sparse data (many zeros)
    counts <- matrix(0, nrow = 100, ncol = 200)
    # Randomly assign some non-zero values
    n_nonzero <- rbinom(1, 100 * 200, 0.3)  # 30% non-zero
    nonzero_indices <- sample(1:(100 * 200), n_nonzero)
    counts[nonzero_indices] <- rlnorm(n_nonzero, meanlog = 1, sdlog = 1)
    metadata <- data.frame(Group = rep(c("Group1", "Group2"), each = 50))
    datasets[[paste0("Type4_", i)]] <- list(counts = counts, metadata = metadata, type = "sparse_data")
  }
  
  # Type 5: Small sample size
  cat("Creating Type 5 null datasets (small samples)...\n")
  for (i in 1:n_datasets_per_type) {
    set.seed(i + 400)
    counts <- matrix(rlnorm(40 * 100, meanlog = 2, sdlog = 1), nrow = 40, ncol = 100)
    metadata <- data.frame(Group = sample(c("Group1", "Group2"), 40, replace = TRUE))
    datasets[[paste0("Type5_", i)]] <- list(counts = counts, metadata = metadata, type = "small_samples")
  }
  
  return(datasets)
}

# Function to run comprehensive Type I error tests
run_comprehensive_type1_tests <- function() {
  cat("🚀 COMPREHENSIVE TYPE I ERROR TESTING\n")
  cat("=====================================\n")
  cat("Testing improved MeLSI across multiple null scenarios\n")
  cat("Goal: Ensure Type I error rate is ~5% across all scenarios\n\n")
  
  # Create diverse null datasets
  datasets <- create_null_datasets(n_datasets_per_type = 15)
  n_datasets <- length(datasets)
  
  cat("Created", n_datasets, "null datasets across 5 different scenarios\n\n")
  
  # Initialize results
  results <- data.frame(
    dataset_name = character(),
    dataset_type = character(),
    n_samples = integer(),
    n_taxa = integer(),
    melsi_p_value = numeric(),
    melsi_significant = logical(),
    bray_p_value = numeric(),
    bray_significant = logical(),
    euclidean_p_value = numeric(),
    euclidean_significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # Test each dataset
  for (i in 1:n_datasets) {
    dataset_name <- names(datasets)[i]
    data <- datasets[[i]]
    
    cat("Testing", dataset_name, "(", data$type, ")...\n")
    cat("  Samples:", nrow(data$counts), "Taxa:", ncol(data$counts), "\n")
    
    # CLR transformation
    counts_clr <- log(data$counts + 1)
    counts_clr <- counts_clr - rowMeans(counts_clr)
    
    # Test improved MeLSI
    melsi_results <- run_melsi_improved(
      counts_clr, data$metadata$Group,
      n_perms = 99, B = 30, m_frac = 0.8,
      show_progress = FALSE
    )
    
    # Test traditional methods for comparison
    # Bray-Curtis
    dist_bray <- vegdist(data$counts, method = "bray")
    permanova_bray <- adonis2(dist_bray ~ Group, data = data$metadata, permutations = 99)
    
    # Euclidean
    dist_euclidean <- vegdist(counts_clr, method = "euclidean")
    permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = data$metadata, permutations = 99)
    
    # Store results
    results <- rbind(results, data.frame(
      dataset_name = dataset_name,
      dataset_type = data$type,
      n_samples = nrow(data$counts),
      n_taxa = ncol(data$counts),
      melsi_p_value = melsi_results$p_value,
      melsi_significant = melsi_results$p_value < 0.05,
      bray_p_value = permanova_bray$`Pr(>F)`[1],
      bray_significant = permanova_bray$`Pr(>F)`[1] < 0.05,
      euclidean_p_value = permanova_euclidean$`Pr(>F)`[1],
      euclidean_significant = permanova_euclidean$`Pr(>F)`[1] < 0.05
    ))
    
    # Show progress
    if (i %% 10 == 0) {
      cat("Completed", i, "of", n_datasets, "datasets\n")
    }
  }
  
  # Analyze results
  cat("\n📊 TYPE I ERROR ANALYSIS\n")
  cat("========================\n")
  
  # Overall Type I error rates
  melsi_false_positives <- sum(results$melsi_significant)
  bray_false_positives <- sum(results$bray_significant)
  euclidean_false_positives <- sum(results$euclidean_significant)
  
  melsi_type1_rate <- melsi_false_positives / n_datasets
  bray_type1_rate <- bray_false_positives / n_datasets
  euclidean_type1_rate <- euclidean_false_positives / n_datasets
  
  cat("Overall Type I Error Rates:\n")
  cat("  MeLSI:", sprintf("%.1f%% (%d/%d)", melsi_type1_rate*100, melsi_false_positives, n_datasets), "\n")
  cat("  Bray-Curtis:", sprintf("%.1f%% (%d/%d)", bray_type1_rate*100, bray_false_positives, n_datasets), "\n")
  cat("  Euclidean:", sprintf("%.1f%% (%d/%d)", euclidean_type1_rate*100, euclidean_false_positives, n_datasets), "\n")
  
  # Type I error rates by scenario
  cat("\nType I Error Rates by Scenario:\n")
  for (scenario in unique(results$dataset_type)) {
    scenario_results <- results[results$dataset_type == scenario, ]
    scenario_melsi_rate <- sum(scenario_results$melsi_significant) / nrow(scenario_results)
    scenario_bray_rate <- sum(scenario_results$bray_significant) / nrow(scenario_results)
    scenario_euclidean_rate <- sum(scenario_results$euclidean_significant) / nrow(scenario_results)
    
    cat("\n", scenario, ":\n")
    cat("  MeLSI:", sprintf("%.1f%% (%d/%d)", scenario_melsi_rate*100, 
                           sum(scenario_results$melsi_significant), nrow(scenario_results)), "\n")
    cat("  Bray-Curtis:", sprintf("%.1f%% (%d/%d)", scenario_bray_rate*100, 
                                 sum(scenario_results$bray_significant), nrow(scenario_results)), "\n")
    cat("  Euclidean:", sprintf("%.1f%% (%d/%d)", scenario_euclidean_rate*100, 
                               sum(scenario_results$euclidean_significant), nrow(scenario_results)), "\n")
  }
  
  # P-value distribution analysis
  cat("\n📈 P-VALUE DISTRIBUTION ANALYSIS\n")
  cat("================================\n")
  
  cat("P-value statistics (should be uniform for null data):\n")
  cat("  MeLSI: mean =", round(mean(results$melsi_p_value), 3), 
      "median =", round(median(results$melsi_p_value), 3), 
      "sd =", round(sd(results$melsi_p_value), 3), "\n")
  cat("  Bray-Curtis: mean =", round(mean(results$bray_p_value), 3), 
      "median =", round(median(results$bray_p_value), 3), 
      "sd =", round(sd(results$bray_p_value), 3), "\n")
  cat("  Euclidean: mean =", round(mean(results$euclidean_p_value), 3), 
      "median =", round(median(results$euclidean_p_value), 3), 
      "sd =", round(sd(results$euclidean_p_value), 3), "\n")
  
  # Statistical assessment
  cat("\n🎯 STATISTICAL ASSESSMENT\n")
  cat("=========================\n")
  
  # Check if Type I error rates are acceptable (5-10%)
  melsi_acceptable <- melsi_type1_rate >= 0.02 && melsi_type1_rate <= 0.10
  bray_acceptable <- bray_type1_rate >= 0.02 && bray_type1_rate <= 0.10
  euclidean_acceptable <- euclidean_type1_rate >= 0.02 && euclidean_type1_rate <= 0.10
  
  cat("Type I Error Rate Assessment (acceptable range: 2-10%):\n")
  cat("  MeLSI:", ifelse(melsi_acceptable, "✅ ACCEPTABLE", "❌ OUTSIDE ACCEPTABLE RANGE"), "\n")
  cat("  Bray-Curtis:", ifelse(bray_acceptable, "✅ ACCEPTABLE", "❌ OUTSIDE ACCEPTABLE RANGE"), "\n")
  cat("  Euclidean:", ifelse(euclidean_acceptable, "✅ ACCEPTABLE", "❌ OUTSIDE ACCEPTABLE RANGE"), "\n")
  
  # Compare MeLSI to traditional methods
  melsi_similar_to_bray <- abs(melsi_type1_rate - bray_type1_rate) <= 0.05
  melsi_similar_to_euclidean <- abs(melsi_type1_rate - euclidean_type1_rate) <= 0.05
  
  cat("\nComparison to Traditional Methods:\n")
  cat("  MeLSI vs Bray-Curtis:", ifelse(melsi_similar_to_bray, "✅ SIMILAR", "⚠️ DIFFERENT"), "\n")
  cat("  MeLSI vs Euclidean:", ifelse(melsi_similar_to_euclidean, "✅ SIMILAR", "⚠️ DIFFERENT"), "\n")
  
  # Overall assessment
  overall_success <- melsi_acceptable && melsi_similar_to_bray && melsi_similar_to_euclidean
  
  cat("\n🏆 OVERALL ASSESSMENT\n")
  cat("====================\n")
  if (overall_success) {
    cat("🎉 IMPROVED MeLSI PASSES COMPREHENSIVE TYPE I ERROR TESTS!\n")
    cat("✅ Type I error rate is acceptable\n")
    cat("✅ Similar to traditional methods\n")
    cat("✅ Statistically sound methodology\n")
  } else {
    cat("❌ IMPROVED MeLSI FAILS TYPE I ERROR TESTS\n")
    if (!melsi_acceptable) cat("❌ Type I error rate outside acceptable range\n")
    if (!melsi_similar_to_bray) cat("❌ Significantly different from Bray-Curtis\n")
    if (!melsi_similar_to_euclidean) cat("❌ Significantly different from Euclidean\n")
  }
  
  # Save detailed results
  write.csv(results, "comprehensive_type1_results.csv", row.names = FALSE)
  cat("\n📁 Detailed results saved to: comprehensive_type1_results.csv\n")
  
  return(list(
    results = results,
    type1_rates = c(
      MeLSI = melsi_type1_rate,
      Bray_Curtis = bray_type1_rate,
      Euclidean = euclidean_type1_rate
    ),
    assessment = overall_success,
    summary = list(
      n_datasets = n_datasets,
      melsi_false_positives = melsi_false_positives,
      bray_false_positives = bray_false_positives,
      euclidean_false_positives = euclidean_false_positives
    )
  ))
}

# Run the comprehensive tests
if (interactive()) {
  comprehensive_results <- run_comprehensive_type1_tests()
}
