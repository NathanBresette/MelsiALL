# Simple Type I Error Test - Just 5 datasets
# Tests the improved MeLSI (with permutation bug fixes) for proper Type I error control

library(vegan)
source("melsi_improved.R")

# Function to create 5 diverse null datasets
create_five_null_datasets <- function() {
  datasets <- list()
  
  cat("Creating 5 null datasets...\n")
  
  # Dataset 1: Random groups, log-normal data
  set.seed(1)
  counts1 <- matrix(rlnorm(100 * 200, meanlog = 2, sdlog = 1), nrow = 100, ncol = 200)
  metadata1 <- data.frame(Group = sample(c("Group1", "Group2"), 100, replace = TRUE))
  datasets[["Random_Groups"]] <- list(counts = counts1, metadata = metadata1, description = "Random groups, log-normal")
  
  # Dataset 2: Balanced groups, log-normal data
  set.seed(2)
  counts2 <- matrix(rlnorm(100 * 200, meanlog = 2, sdlog = 1), nrow = 100, ncol = 200)
  metadata2 <- data.frame(Group = rep(c("Group1", "Group2"), each = 50))
  counts2 <- counts2[sample(1:100), ]  # Shuffle to break structure
  datasets[["Balanced_Groups"]] <- list(counts = counts2, metadata = metadata2, description = "Balanced groups, shuffled")
  
  # Dataset 3: Random groups, normal data
  set.seed(3)
  counts3 <- matrix(rnorm(100 * 200, mean = 10, sd = 3), nrow = 100, ncol = 200)
  counts3 <- pmax(1, counts3)  # Ensure positive
  metadata3 <- data.frame(Group = sample(c("Group1", "Group2"), 100, replace = TRUE))
  datasets[["Normal_Data"]] <- list(counts = counts3, metadata = metadata3, description = "Random groups, normal data")
  
  # Dataset 4: Sparse data (many zeros)
  set.seed(4)
  counts4 <- matrix(0, nrow = 100, ncol = 200)
  n_nonzero <- rbinom(1, 100 * 200, 0.3)  # 30% non-zero
  nonzero_indices <- sample(1:(100 * 200), n_nonzero)
  counts4[nonzero_indices] <- rlnorm(n_nonzero, meanlog = 1, sdlog = 1)
  metadata4 <- data.frame(Group = sample(c("Group1", "Group2"), 100, replace = TRUE))
  datasets[["Sparse_Data"]] <- list(counts = counts4, metadata = metadata4, description = "Sparse data, many zeros")
  
  # Dataset 5: Small sample size
  set.seed(5)
  counts5 <- matrix(rlnorm(50 * 100, meanlog = 2, sdlog = 1), nrow = 50, ncol = 100)
  metadata5 <- data.frame(Group = sample(c("Group1", "Group2"), 50, replace = TRUE))
  datasets[["Small_Sample"]] <- list(counts = counts5, metadata = metadata5, description = "Small sample size (50)")
  
  return(datasets)
}

# Function to run simple Type I error test
run_simple_type1_test <- function() {
  cat("🚀 SIMPLE TYPE I ERROR TEST\n")
  cat("===========================\n")
  cat("Testing improved MeLSI on 5 null datasets\n")
  cat("Goal: Verify Type I error rate is ~5%\n\n")
  
  # Create 5 null datasets
  datasets <- create_five_null_datasets()
  
  # Initialize results
  results <- data.frame(
    dataset = character(),
    description = character(),
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
  for (i in 1:length(datasets)) {
    dataset_name <- names(datasets)[i]
    data <- datasets[[i]]
    
    cat("Testing", dataset_name, "-", data$description, "\n")
    cat("  Samples:", nrow(data$counts), "Taxa:", ncol(data$counts), "\n")
    
    # CLR transformation
    counts_clr <- log(data$counts + 1)
    counts_clr <- counts_clr - rowMeans(counts_clr)
    
    # Test improved MeLSI
    cat("  Running MeLSI...\n")
    melsi_results <- run_melsi_improved(
      counts_clr, data$metadata$Group,
      n_perms = 99, B = 30, m_frac = 0.8,
      show_progress = FALSE
    )
    
    # Test traditional methods for comparison
    cat("  Running traditional methods...\n")
    
    # Bray-Curtis
    dist_bray <- vegdist(data$counts, method = "bray")
    permanova_bray <- adonis2(dist_bray ~ Group, data = data$metadata, permutations = 99)
    
    # Euclidean
    dist_euclidean <- vegdist(counts_clr, method = "euclidean")
    permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = data$metadata, permutations = 99)
    
    # Store results
    results <- rbind(results, data.frame(
      dataset = dataset_name,
      description = data$description,
      n_samples = nrow(data$counts),
      n_taxa = ncol(data$counts),
      melsi_p_value = melsi_results$p_value,
      melsi_significant = melsi_results$p_value < 0.05,
      bray_p_value = permanova_bray$`Pr(>F)`[1],
      bray_significant = permanova_bray$`Pr(>F)`[1] < 0.05,
      euclidean_p_value = permanova_euclidean$`Pr(>F)`[1],
      euclidean_significant = permanova_euclidean$`Pr(>F)`[1] < 0.05
    ))
    
    # Show current results
    cat("  Results:\n")
    cat("    MeLSI: p =", round(melsi_results$p_value, 4), 
        ifelse(melsi_results$p_value < 0.05, "*** (FALSE POSITIVE)", ""), "\n")
    cat("    Bray-Curtis: p =", round(permanova_bray$`Pr(>F)`[1], 4), 
        ifelse(permanova_bray$`Pr(>F)`[1] < 0.05, "*** (FALSE POSITIVE)", ""), "\n")
    cat("    Euclidean: p =", round(permanova_euclidean$`Pr(>F)`[1], 4), 
        ifelse(permanova_euclidean$`Pr(>F)`[1] < 0.05, "*** (FALSE POSITIVE)", ""), "\n")
    cat("\n")
  }
  
  # Summary analysis
  cat("📊 TYPE I ERROR SUMMARY\n")
  cat("=======================\n")
  
  n_datasets <- nrow(results)
  melsi_false_positives <- sum(results$melsi_significant)
  bray_false_positives <- sum(results$bray_significant)
  euclidean_false_positives <- sum(results$euclidean_significant)
  
  melsi_type1_rate <- melsi_false_positives / n_datasets
  bray_type1_rate <- bray_false_positives / n_datasets
  euclidean_type1_rate <- euclidean_false_positives / n_datasets
  
  cat("False Positive Rates (should be ~5%):\n")
  cat("  MeLSI:", sprintf("%.1f%% (%d/%d)", melsi_type1_rate*100, melsi_false_positives, n_datasets), "\n")
  cat("  Bray-Curtis:", sprintf("%.1f%% (%d/%d)", bray_type1_rate*100, bray_false_positives, n_datasets), "\n")
  cat("  Euclidean:", sprintf("%.1f%% (%d/%d)", euclidean_type1_rate*100, euclidean_false_positives, n_datasets), "\n")
  
  # P-value statistics
  cat("\nP-value Statistics:\n")
  cat("  MeLSI: mean =", round(mean(results$melsi_p_value), 3), 
      "median =", round(median(results$melsi_p_value), 3), "\n")
  cat("  Bray-Curtis: mean =", round(mean(results$bray_p_value), 3), 
      "median =", round(median(results$bray_p_value), 3), "\n")
  cat("  Euclidean: mean =", round(mean(results$euclidean_p_value), 3), 
      "median =", round(median(results$euclidean_p_value), 3), "\n")
  
  # Assessment
  cat("\n🎯 ASSESSMENT\n")
  cat("=============\n")
  
  # Check if Type I error rates are acceptable (0-20% for 5 tests)
  melsi_acceptable <- melsi_type1_rate <= 0.20  # Allow up to 1/5 false positive
  bray_acceptable <- bray_type1_rate <= 0.20
  euclidean_acceptable <- euclidean_type1_rate <= 0.20
  
  cat("Type I Error Rate Assessment (acceptable: ≤20% for 5 tests):\n")
  cat("  MeLSI:", ifelse(melsi_acceptable, "✅ ACCEPTABLE", "❌ TOO HIGH"), "\n")
  cat("  Bray-Curtis:", ifelse(bray_acceptable, "✅ ACCEPTABLE", "❌ TOO HIGH"), "\n")
  cat("  Euclidean:", ifelse(euclidean_acceptable, "✅ ACCEPTABLE", "❌ TOO HIGH"), "\n")
  
  # Compare to traditional methods
  melsi_similar <- abs(melsi_type1_rate - bray_type1_rate) <= 0.20
  
  cat("\nComparison to Traditional Methods:\n")
  cat("  MeLSI vs Bray-Curtis:", ifelse(melsi_similar, "✅ SIMILAR", "⚠️ DIFFERENT"), "\n")
  
  # Overall assessment
  overall_success <- melsi_acceptable && melsi_similar
  
  cat("\n🏆 OVERALL ASSESSMENT\n")
  cat("====================\n")
  if (overall_success) {
    cat("🎉 IMPROVED MeLSI PASSES TYPE I ERROR TEST!\n")
    cat("✅ Type I error rate is acceptable\n")
    cat("✅ Similar to traditional methods\n")
    cat("✅ Ready for publication!\n")
  } else {
    cat("❌ IMPROVED MeLSI NEEDS MORE WORK\n")
    if (!melsi_acceptable) cat("❌ Type I error rate too high\n")
    if (!melsi_similar) cat("❌ Too different from traditional methods\n")
  }
  
  # Save results
  write.csv(results, "simple_type1_results.csv", row.names = FALSE)
  cat("\n📁 Results saved to: simple_type1_results.csv\n")
  
  return(list(
    results = results,
    type1_rates = c(
      MeLSI = melsi_type1_rate,
      Bray_Curtis = bray_type1_rate,
      Euclidean = euclidean_type1_rate
    ),
    assessment = overall_success
  ))
}

# Run the simple test
if (interactive()) {
  simple_results <- run_simple_type1_test()
}
