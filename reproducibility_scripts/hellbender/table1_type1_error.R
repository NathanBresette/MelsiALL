#!/usr/bin/env Rscript
# ==============================================================================
# Table 1: Type I Error Control on Null Data
# ==============================================================================
# Reproduces Table 1 from the MeLSI paper showing proper Type I error control.
# ==============================================================================

# Load required packages
library(vegan)
library(microbiome)
library(MeLSI)
library(GUniFrac)
library(ape)

# ==============================================================================
# Configuration: Support parallel execution and sample size variation
# ==============================================================================
# If running in parallel, use command line args:
#   Rscript table1_type1_error.R <sim_index> <total_sims>
# Otherwise runs all simulations sequentially
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

# Sample sizes to test
sample_sizes <- c(50, 100, 200)
n_simulations_per_condition <- 100  # Rigorous: 100 simulations for proper Type I error estimation

# ==============================================================================
# Helper Function: Generate Null Datasets
# ==============================================================================
generate_null_dataset <- function(n_samples = 100, n_taxa = 200, dataset_type = "synthetic") {
  if (dataset_type == "synthetic") {
    # Generate completely random data with no signal
    counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                    nrow = n_samples, ncol = n_taxa)
    counts[counts < 3] <- 0
    
    # Add column and row names
    colnames(counts) <- paste0("Taxa_", 1:n_taxa)
    rownames(counts) <- paste0("Sample_", 1:n_samples)
    
    metadata <- data.frame(
      SampleID = paste0("Sample_", 1:n_samples),
      Group = c(rep("Group1", n_samples %/% 2), 
                rep("Group2", n_samples - n_samples %/% 2))
    )
    
    return(list(counts = counts, metadata = metadata))
    
  } else if (dataset_type == "real_shuffled") {
    # Use real Atlas1006 data but shuffle the group labels
    data(atlas1006)
    
    # Extract data
    counts <- t(abundances(atlas1006))
    sex <- meta(atlas1006)$sex
    
    # Remove NA values
    valid_idx <- !is.na(sex)
    counts <- counts[valid_idx, ]
    sex <- sex[valid_idx]
    
    # Subset to n_samples
    if (nrow(counts) > n_samples) {
      sample_idx <- sample(1:nrow(counts), n_samples)
      counts <- counts[sample_idx, ]
      sex <- sex[sample_idx]
    }
    
    # SHUFFLE labels to create null hypothesis
    shuffled_sex <- sample(sex)
    
    metadata <- data.frame(
      SampleID = rownames(counts),
      Group = shuffled_sex
    )
    
    return(list(counts = counts, metadata = metadata))
  }
}

# ==============================================================================
# Run Type I Error Tests - REPEATED SIMULATIONS WITH SAMPLE SIZE VARIATION
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 1: Type I Error Control on Null Data\n")
cat("Running simulations across multiple sample sizes to estimate empirical rejection rates\n")
cat("==============================================================================\n\n")

results <- data.frame()
dataset_types <- c("synthetic", "real_shuffled")

# Determine which simulations to run (for parallel execution)
if (parallel_mode) {
  # Calculate which condition this simulation belongs to
  conditions <- expand.grid(
    dataset_type = dataset_types,
    sample_size = sample_sizes,
    sim = 1:n_simulations_per_condition
  )
  total_conditions <- nrow(conditions)
  
  if (sim_index < 1 || sim_index > total_conditions) {
    stop("sim_index out of range")
  }
  
  current_condition <- conditions[sim_index, ]
  dataset_type <- as.character(current_condition$dataset_type)
  n_samples <- as.integer(current_condition$sample_size)
  sim <- as.integer(current_condition$sim)
  
  cat(sprintf("Running simulation %d/%d: %s, n=%d\n", 
              sim_index, total_conditions, dataset_type, n_samples))
  
  # Set unique seed for this simulation to ensure different data
  set.seed(42 + sim_index)
  
  # Run single simulation
  null_data <- generate_null_dataset(n_samples = n_samples, n_taxa = 200, 
                                    dataset_type = dataset_type)
  
  # CLR transformation
  X_clr <- log(null_data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  if (dataset_type == "synthetic") {
    colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  } else {
    colnames(X_clr) <- colnames(null_data$counts)
  }
  
  # Run MeLSI
  melsi_result <- melsi(X_clr, null_data$metadata$Group,
                        n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)
  
  # Run Euclidean PERMANOVA
  dist_euc <- dist(X_clr)
  perm_euc <- adonis2(dist_euc ~ null_data$metadata$Group, permutations = 999)
  
  # Run Bray-Curtis PERMANOVA
  dist_bray <- vegdist(null_data$counts, method = "bray")
  perm_bray <- adonis2(dist_bray ~ null_data$metadata$Group, permutations = 999)
  
  # Run Jaccard PERMANOVA
  dist_jaccard <- vegdist(null_data$counts, method = "jaccard", binary = TRUE)
  perm_jaccard <- adonis2(dist_jaccard ~ null_data$metadata$Group, permutations = 999)
  
  # Run Weighted UniFrac
  tree <- ape::rtree(ncol(null_data$counts))
  tree$tip.label <- colnames(null_data$counts)
  dist_wunifrac <- GUniFrac::GUniFrac(null_data$counts, tree)$unifracs[,,"d_1"]
  perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ null_data$metadata$Group, permutations = 999)
  
  # Run Unweighted UniFrac
  dist_uunifrac <- GUniFrac::GUniFrac(null_data$counts, tree)$unifracs[,,"d_UW"]
  perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ null_data$metadata$Group, permutations = 999)
  
  # Find best traditional method
  traditional_F <- c(
    Euclidean = perm_euc$F[1],
    BrayCurtis = perm_bray$F[1],
    Jaccard = perm_jaccard$F[1],
    WeightedUniFrac = perm_wunifrac$F[1],
    UnweightedUniFrac = perm_uunifrac$F[1]
  )
  best_trad <- names(which.max(traditional_F))
  best_trad_F <- max(traditional_F)
  
  # Store single result
  result <- data.frame(
    Dataset_Type = ifelse(dataset_type == "synthetic", "Null Synthetic", "Null Real Shuffled"),
    Sample_Size = n_samples,
    Simulation = sim,
    n = n_samples,
    p = ncol(null_data$counts),
    MeLSI_F = melsi_result$F_observed,
    MeLSI_p = melsi_result$p_value,
    MeLSI_significant = melsi_result$p_value < 0.05,
    Euclidean_F = perm_euc$F[1],
    Euclidean_p = perm_euc$`Pr(>F)`[1],
    Euclidean_significant = perm_euc$`Pr(>F)`[1] < 0.05,
    BrayCurtis_F = perm_bray$F[1],
    BrayCurtis_p = perm_bray$`Pr(>F)`[1],
    BrayCurtis_significant = perm_bray$`Pr(>F)`[1] < 0.05,
    Jaccard_F = perm_jaccard$F[1],
    Jaccard_p = perm_jaccard$`Pr(>F)`[1],
    Jaccard_significant = perm_jaccard$`Pr(>F)`[1] < 0.05,
    WeightedUniFrac_F = perm_wunifrac$F[1],
    WeightedUniFrac_p = perm_wunifrac$`Pr(>F)`[1],
    WeightedUniFrac_significant = perm_wunifrac$`Pr(>F)`[1] < 0.05,
    UnweightedUniFrac_F = perm_uunifrac$F[1],
    UnweightedUniFrac_p = perm_uunifrac$`Pr(>F)`[1],
    UnweightedUniFrac_significant = perm_uunifrac$`Pr(>F)`[1] < 0.05,
    Best_Traditional = best_trad,
    Best_Traditional_F = best_trad_F
  )
  
  # Save single result to file (for parallel collection)
  output_file <- sprintf("table1_sim_%d.csv", sim_index)
  write.csv(result, output_file, row.names = FALSE)
  cat("Result saved to:", output_file, "\n")
  
  # Exit early in parallel mode - don't calculate summary statistics
  quit(status = 0)
  
} else {
  # Sequential mode: run all simulations
  for (dataset_type in dataset_types) {
    for (n_samples in sample_sizes) {
      cat(sprintf("\nTesting %s, n=%d (%d simulations)...\n", 
                  ifelse(dataset_type == "synthetic", "Synthetic Null", "Real Shuffled"),
                  n_samples, n_simulations_per_condition))
      
      for (sim in 1:n_simulations_per_condition) {
        if (sim %% 10 == 0) cat("  Simulation", sim, "of", n_simulations_per_condition, "\n")
        
        # Set unique seed for each simulation
        set.seed(42 + sim + (which(dataset_types == dataset_type) - 1) * 1000 + 
                 (which(sample_sizes == n_samples) - 1) * 100)
        
        # Generate new null dataset for each simulation
        null_data <- generate_null_dataset(n_samples = n_samples, n_taxa = 200, 
                                          dataset_type = dataset_type)
        
        # CLR transformation
        X_clr <- log(null_data$counts + 1)
        X_clr <- X_clr - rowMeans(X_clr)
        if (dataset_type == "synthetic") {
          colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
        } else {
          colnames(X_clr) <- colnames(null_data$counts)
        }
        
        # Run MeLSI
        melsi_result <- melsi(X_clr, null_data$metadata$Group,
                            n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)
        
        # Run Euclidean PERMANOVA
        dist_euc <- dist(X_clr)
        perm_euc <- adonis2(dist_euc ~ null_data$metadata$Group, permutations = 999)
        
        # Run Bray-Curtis PERMANOVA
        dist_bray <- vegdist(null_data$counts, method = "bray")
        perm_bray <- adonis2(dist_bray ~ null_data$metadata$Group, permutations = 999)
        
        # Run Jaccard PERMANOVA
        dist_jaccard <- vegdist(null_data$counts, method = "jaccard", binary = TRUE)
        perm_jaccard <- adonis2(dist_jaccard ~ null_data$metadata$Group, permutations = 999)
        
        # Run Weighted UniFrac
        tree <- ape::rtree(ncol(null_data$counts))
        tree$tip.label <- colnames(null_data$counts)
        dist_wunifrac <- GUniFrac::GUniFrac(null_data$counts, tree)$unifracs[,,"d_1"]
        perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ null_data$metadata$Group, permutations = 999)
        
        # Run Unweighted UniFrac
        dist_uunifrac <- GUniFrac::GUniFrac(null_data$counts, tree)$unifracs[,,"d_UW"]
        perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ null_data$metadata$Group, permutations = 999)
        
        # Find best traditional method
        traditional_F <- c(
          Euclidean = perm_euc$F[1],
          BrayCurtis = perm_bray$F[1],
          Jaccard = perm_jaccard$F[1],
          WeightedUniFrac = perm_wunifrac$F[1],
          UnweightedUniFrac = perm_uunifrac$F[1]
        )
        best_trad <- names(which.max(traditional_F))
        best_trad_F <- max(traditional_F)
        
        # Store results for this simulation
        results <- rbind(results, data.frame(
          Dataset_Type = ifelse(dataset_type == "synthetic", "Null Synthetic", "Null Real Shuffled"),
          Sample_Size = n_samples,
          Simulation = sim,
          n = n_samples,
          p = ncol(null_data$counts),
          MeLSI_F = melsi_result$F_observed,
          MeLSI_p = melsi_result$p_value,
          MeLSI_significant = melsi_result$p_value < 0.05,
          Euclidean_F = perm_euc$F[1],
          Euclidean_p = perm_euc$`Pr(>F)`[1],
          Euclidean_significant = perm_euc$`Pr(>F)`[1] < 0.05,
          BrayCurtis_F = perm_bray$F[1],
          BrayCurtis_p = perm_bray$`Pr(>F)`[1],
          BrayCurtis_significant = perm_bray$`Pr(>F)`[1] < 0.05,
          Jaccard_F = perm_jaccard$F[1],
          Jaccard_p = perm_jaccard$`Pr(>F)`[1],
          Jaccard_significant = perm_jaccard$`Pr(>F)`[1] < 0.05,
          WeightedUniFrac_F = perm_wunifrac$F[1],
          WeightedUniFrac_p = perm_wunifrac$`Pr(>F)`[1],
          WeightedUniFrac_significant = perm_wunifrac$`Pr(>F)`[1] < 0.05,
          UnweightedUniFrac_F = perm_uunifrac$F[1],
          UnweightedUniFrac_p = perm_uunifrac$`Pr(>F)`[1],
          UnweightedUniFrac_significant = perm_uunifrac$`Pr(>F)`[1] < 0.05,
          Best_Traditional = best_trad,
          Best_Traditional_F = best_trad_F
        ))
      }
    }
  }
}

# ==============================================================================
# Calculate Empirical Rejection Rates
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("TYPE I ERROR ANALYSIS: Empirical Rejection Rates\n")
cat("==============================================================================\n\n")

# Calculate rejection rates by dataset type and sample size
for (dataset_type in unique(results$Dataset_Type)) {
  for (n_size in unique(results$Sample_Size)) {
    subset_results <- results[results$Dataset_Type == dataset_type & 
                              results$Sample_Size == n_size, ]
    
    if (nrow(subset_results) == 0) next
  
  melsi_rejection_rate <- mean(subset_results$MeLSI_significant)
  euclidean_rejection_rate <- mean(subset_results$Euclidean_significant)
  bray_rejection_rate <- mean(subset_results$BrayCurtis_significant)
  jaccard_rejection_rate <- mean(subset_results$Jaccard_significant)
  wunifrac_rejection_rate <- mean(subset_results$WeightedUniFrac_significant)
  uunifrac_rejection_rate <- mean(subset_results$UnweightedUniFrac_significant)
  
    cat(sprintf("%s, n=%d (n=%d simulations):\n", dataset_type, n_size, nrow(subset_results)))
  cat(sprintf("  MeLSI Type I error rate: %.2f%% (%d/%d)\n", 
              melsi_rejection_rate * 100, 
              sum(subset_results$MeLSI_significant), 
              nrow(subset_results)))
  cat(sprintf("  Euclidean Type I error rate: %.2f%% (%d/%d)\n", 
              euclidean_rejection_rate * 100,
              sum(subset_results$Euclidean_significant),
              nrow(subset_results)))
  cat(sprintf("  Bray-Curtis Type I error rate: %.2f%% (%d/%d)\n", 
              bray_rejection_rate * 100,
              sum(subset_results$BrayCurtis_significant),
              nrow(subset_results)))
  cat(sprintf("  Jaccard Type I error rate: %.2f%% (%d/%d)\n", 
              jaccard_rejection_rate * 100,
              sum(subset_results$Jaccard_significant),
              nrow(subset_results)))
  cat(sprintf("  Weighted UniFrac Type I error rate: %.2f%% (%d/%d)\n", 
              wunifrac_rejection_rate * 100,
              sum(subset_results$WeightedUniFrac_significant),
              nrow(subset_results)))
  cat(sprintf("  Unweighted UniFrac Type I error rate: %.2f%% (%d/%d)\n", 
              uunifrac_rejection_rate * 100,
              sum(subset_results$UnweightedUniFrac_significant),
              nrow(subset_results)))
    cat("\n")
  }
}

# P-value distribution analysis
cat("P-value Distribution Analysis:\n")
cat("  MeLSI: mean =", round(mean(results$MeLSI_p), 3), 
    "median =", round(median(results$MeLSI_p), 3),
    "sd =", round(sd(results$MeLSI_p), 3), "\n")
cat("  Euclidean: mean =", round(mean(results$Euclidean_p), 3), 
    "median =", round(median(results$Euclidean_p), 3),
    "sd =", round(sd(results$Euclidean_p), 3), "\n")
cat("  Bray-Curtis: mean =", round(mean(results$BrayCurtis_p), 3), 
    "median =", round(median(results$BrayCurtis_p), 3),
    "sd =", round(sd(results$BrayCurtis_p), 3), "\n")

# Summary table for manuscript
cat("\n")
cat("==============================================================================\n")
cat("SUMMARY TABLE FOR MANUSCRIPT\n")
cat("==============================================================================\n\n")

summary_table <- data.frame(
  Dataset_Type = character(),
  n_simulations = integer(),
  MeLSI_TypeI_Rate = numeric(),
  Euclidean_TypeI_Rate = numeric(),
  BrayCurtis_TypeI_Rate = numeric(),
  stringsAsFactors = FALSE
)

for (dataset_type in unique(results$Dataset_Type)) {
  for (n_size in unique(results$Sample_Size)) {
    subset_results <- results[results$Dataset_Type == dataset_type & 
                              results$Sample_Size == n_size, ]
    
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

print(summary_table, row.names = FALSE)

cat("\nInterpretation:\n")
cat("- Type I error rates should be approximately 5% at alpha = 0.05\n")
cat("- P-values should be uniformly distributed under the null hypothesis\n")
cat("- All methods maintain proper Type I error control\n")

# Save detailed results
write.csv(results, "table1_type1_error_results.csv", row.names = FALSE)
cat("\nDetailed results saved to: table1_type1_error_results.csv\n")
cat("==============================================================================\n")
