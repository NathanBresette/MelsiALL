#!/usr/bin/env Rscript
# ==============================================================================
# Table 6: Effect of Feature Correlation on MeLSI Performance
# ==============================================================================
# Validates MeLSI performance when features are correlated, addressing the
# reviewer concern that synthetic data assumed independent taxa.
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)
library(GUniFrac)
library(ape)
library(MASS)  # For mvrnorm to generate correlated data

# ==============================================================================
# Configuration: Support parallel execution
# ==============================================================================
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

# Simulation parameters
n_simulations_per_condition <- 50  # Rigorous: 50 simulations per condition

# Correlation levels to test
correlation_levels <- c(0.0, 0.3, 0.6)  # None, Low, Moderate
correlation_names <- c("None", "Low", "Moderate")

# Fixed conditions for focused analysis
n_samples <- 100  # Medium sample size
n_taxa <- 200
effect_size <- "medium"  # Medium effect size

# ==============================================================================
# Helper Function: Generate Correlated Microbiome Data
# ==============================================================================
generate_correlated_dataset <- function(n_samples = 100, n_taxa = 200, 
                                       effect_size = "medium", 
                                       correlation_level = 0.0) {
  # Generate base counts with independent structure
  base_counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                       nrow = n_samples, ncol = n_taxa)
  base_counts[base_counts < 3] <- 0
  
  # Add correlation structure if correlation_level > 0
  if (correlation_level > 0) {
    # Create correlation matrix: block structure where taxa within blocks are correlated
    # This mimics realistic microbiome correlation (e.g., co-occurring taxa)
    n_blocks <- 10  # 10 blocks of correlated taxa
    taxa_per_block <- n_taxa %/% n_blocks
    
    # Generate correlated log-abundances using multivariate normal
    # Then convert back to counts
    log_base <- log(base_counts + 1)
    
    for (block in 1:n_blocks) {
      block_start <- (block - 1) * taxa_per_block + 1
      block_end <- min(block * taxa_per_block, n_taxa)
      block_indices <- block_start:block_end
      
      if (length(block_indices) > 1) {
        # Create correlation matrix for this block
        block_cor <- matrix(correlation_level, 
                           nrow = length(block_indices), 
                           ncol = length(block_indices))
        diag(block_cor) <- 1
        
        # Generate correlated noise for this block
        # Use Cholesky decomposition to generate correlated multivariate normal
        tryCatch({
          block_sd <- apply(log_base[, block_indices, drop = FALSE], 2, sd)
          block_sd[block_sd == 0] <- 1  # Avoid division by zero
          
          # Generate correlated residuals
          correlated_noise <- mvrnorm(n_samples, 
                                     mu = rep(0, length(block_indices)),
                                     Sigma = block_cor * outer(block_sd, block_sd))
          
          # Add correlated noise to log abundances
          log_base[, block_indices] <- log_base[, block_indices] + 
            correlated_noise * 0.3  # Scale to avoid overwhelming signal
        }, error = function(e) {
          # If correlation matrix is not positive definite, skip correlation for this block
          # This can happen with high correlation and small blocks
        })
      }
    }
    
    # Convert back to counts (ensure non-negative and preserve matrix structure)
    base_counts <- matrix(pmax(0, round(exp(log_base) - 1)), 
                         nrow = n_samples, ncol = n_taxa)
  }
  
  # Add signal based on effect size
  group1_samples <- n_samples %/% 2
  group2_samples <- n_samples - group1_samples
  
  if (effect_size == "small") {
    signal_taxa <- 5
    fold_change <- 1.5
  } else if (effect_size == "medium") {
    signal_taxa <- 10
    fold_change <- 2.0
  } else if (effect_size == "large") {
    signal_taxa <- 20
    fold_change <- 3.0
  }
  
  # Add differential abundance to first group
  signal_indices <- sample(1:n_taxa, signal_taxa)
  for (i in 1:group1_samples) {
    base_counts[i, signal_indices] <- base_counts[i, signal_indices] * fold_change
  }
  
  # Add column names for taxa
  colnames(base_counts) <- paste0("Taxa_", 1:n_taxa)
  rownames(base_counts) <- paste0("Sample_", 1:n_samples)
  
  metadata <- data.frame(
    SampleID = paste0("Sample_", 1:n_samples),
    Group = c(rep("Group1", group1_samples), rep("Group2", group2_samples))
  )
  
  return(list(counts = base_counts, metadata = metadata, true_signal_taxa = signal_indices))
}

# ==============================================================================
# Helper Function: Calculate Recovery Metrics (reuse from table2)
# ==============================================================================
calculate_recovery_metrics <- function(melsi_result, true_signal_taxa, X_clr, groups) {
  feature_weights <- diag(melsi_result$metric_matrix)
  n_filtered <- length(feature_weights)
  n_total <- ncol(X_clr)
  
  # Re-apply pre-filtering to determine which features were kept
  group1_idx <- which(groups == unique(groups)[1])
  group2_idx <- which(groups == unique(groups)[2])
  
  importance_scores <- numeric(n_total)
  for (j in 1:n_total) {
    mu1 <- mean(X_clr[group1_idx, j])
    mu2 <- mean(X_clr[group2_idx, j])
    sigma1_sq <- var(X_clr[group1_idx, j])
    sigma2_sq <- var(X_clr[group2_idx, j])
    denominator <- sqrt(sigma1_sq + sigma2_sq + 1e-6)
    importance_scores[j] <- abs(mu1 - mu2) / denominator
  }
  
  n_keep <- max(10, floor(n_total * 0.7))
  kept_features <- order(importance_scores, decreasing = TRUE)[1:n_keep]
  
  signal_in_filtered <- match(true_signal_taxa, kept_features)
  signal_in_filtered <- signal_in_filtered[!is.na(signal_in_filtered)]
  
  ranked_filtered <- order(feature_weights, decreasing = TRUE)
  
  k_values <- c(5, 10, 20)
  metrics <- list()
  
  for (k in k_values) {
    if (k > n_filtered) k <- n_filtered
    
    signal_ranks <- match(signal_in_filtered, ranked_filtered)
    signal_ranks <- signal_ranks[!is.na(signal_ranks)]
    
    if (length(signal_ranks) > 0) {
      precision_k <- sum(signal_ranks <= k) / k
      recall_k <- sum(signal_ranks <= k) / length(true_signal_taxa)
    } else {
      precision_k <- 0
      recall_k <- 0
    }
    
    metrics[[paste0("precision_", k)]] <- precision_k
    metrics[[paste0("recall_", k)]] <- recall_k
  }
  
  signal_ranks <- match(signal_in_filtered, ranked_filtered)
  signal_ranks <- signal_ranks[!is.na(signal_ranks)]
  mean_rank <- if (length(signal_ranks) > 0) mean(signal_ranks) else NA
  
  is_signal <- rep(0, n_filtered)
  if (length(signal_in_filtered) > 0) {
    is_signal[signal_in_filtered] <- 1
  }
  
  if (sum(is_signal) > 0 && sum(is_signal) < length(is_signal)) {
    sorted_idx <- order(feature_weights, decreasing = TRUE)
    sorted_labels <- is_signal[sorted_idx]
    
    n_pos <- sum(sorted_labels)
    n_neg <- length(sorted_labels) - n_pos
    if (n_pos > 0 && n_neg > 0) {
      tp <- cumsum(sorted_labels)
      fp <- cumsum(1 - sorted_labels)
      tpr <- tp / n_pos
      fpr <- fp / n_neg
      auc_roc <- sum(diff(c(0, fpr)) * tpr[-length(tpr)])
    } else {
      auc_roc <- NA
    }
  } else {
    auc_roc <- NA
  }
  
  return(list(
    precision_5 = metrics$precision_5,
    precision_10 = metrics$precision_10,
    precision_20 = metrics$precision_20,
    recall_5 = metrics$recall_5,
    recall_10 = metrics$recall_10,
    recall_20 = metrics$recall_20,
    mean_rank = mean_rank,
    auc_roc = auc_roc
  ))
}

# ==============================================================================
# Run Correlation Analysis - REPEATED SIMULATIONS
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("TABLE 6: Effect of Feature Correlation on MeLSI Performance\n")
cat("Running multiple simulations per correlation level\n")
cat("==============================================================================\n\n")

results <- data.frame()

# Determine which simulation to run (for parallel execution)
if (parallel_mode) {
  # Calculate which condition this simulation belongs to
  conditions <- expand.grid(
    correlation_level = correlation_levels,
    sim = 1:n_simulations_per_condition
  )
  total_conditions <- nrow(conditions)
  
  if (sim_index < 1 || sim_index > total_conditions) {
    stop("sim_index out of range")
  }
  
  current_condition <- conditions[sim_index, ]
  corr_level <- as.numeric(current_condition$correlation_level)
  sim <- as.integer(current_condition$sim)
  corr_name <- correlation_names[which(correlation_levels == corr_level)]
  
  cat(sprintf("Running simulation %d/%d: %s correlation (%.1f)\n", 
              sim_index, total_conditions, corr_name, corr_level))
  
  # Set unique seed for this simulation
  set.seed(42 + sim_index)
  
  # Generate correlated data
  corr_data <- generate_correlated_dataset(
    n_samples = n_samples,
    n_taxa = n_taxa,
    effect_size = effect_size,
    correlation_level = corr_level
  )
  
  # CLR transformation
  X_clr <- log(corr_data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
  
  # Run MeLSI
  melsi_result <- melsi(X_clr, corr_data$metadata$Group,
                       n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)
  
  # Calculate recovery metrics
  recovery_metrics <- calculate_recovery_metrics(
    melsi_result,
    corr_data$true_signal_taxa,
    X_clr,
    corr_data$metadata$Group
  )
  
  # Run Euclidean PERMANOVA
  dist_euc <- dist(X_clr)
  perm_euc <- adonis2(dist_euc ~ corr_data$metadata$Group, permutations = 999)
  
  # Run Bray-Curtis PERMANOVA
  dist_bray <- vegdist(corr_data$counts, method = "bray")
  perm_bray <- adonis2(dist_bray ~ corr_data$metadata$Group, permutations = 999)
  
  # Run Jaccard PERMANOVA
  dist_jaccard <- vegdist(corr_data$counts, method = "jaccard", binary = TRUE)
  perm_jaccard <- adonis2(dist_jaccard ~ corr_data$metadata$Group, permutations = 999)
  
  # Run Weighted UniFrac (using random tree for synthetic data)
  tree <- ape::rtree(ncol(corr_data$counts))
  tree$tip.label <- colnames(corr_data$counts)
  dist_wunifrac <- GUniFrac::GUniFrac(corr_data$counts, tree)$unifracs[,,"d_1"]
  perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ corr_data$metadata$Group, permutations = 999)
  
  # Run Unweighted UniFrac
  dist_uunifrac <- GUniFrac::GUniFrac(corr_data$counts, tree)$unifracs[,,"d_UW"]
  perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ corr_data$metadata$Group, permutations = 999)
  
  # Find best traditional method
  traditional_F <- c(Euclidean = perm_euc$F[1], 
                     BrayCurtis = perm_bray$F[1],
                     Jaccard = perm_jaccard$F[1],
                     WeightedUniFrac = perm_wunifrac$F[1],
                     UnweightedUniFrac = perm_uunifrac$F[1])
  best_trad <- names(which.max(traditional_F))
  best_trad_F <- max(traditional_F)
  
  # Store single result
  result <- data.frame(
    Correlation_Level = corr_name,
    Correlation_Value = corr_level,
    Sample_Size = n_samples,
    Simulation = sim,
    MeLSI_F = melsi_result$F_observed,
    MeLSI_p = melsi_result$p_value,
    MeLSI_significant = melsi_result$p_value < 0.05,
    # Recovery metrics
    Precision_5 = recovery_metrics$precision_5,
    Precision_10 = recovery_metrics$precision_10,
    Precision_20 = recovery_metrics$precision_20,
    Recall_5 = recovery_metrics$recall_5,
    Recall_10 = recovery_metrics$recall_10,
    Recall_20 = recovery_metrics$recall_20,
    Mean_Rank = recovery_metrics$mean_rank,
    AUC_ROC = recovery_metrics$auc_roc,
    # Traditional methods
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
  output_file <- sprintf("table6_sim_%d.csv", sim_index)
  write.csv(result, output_file, row.names = FALSE)
  cat("Result saved to:", output_file, "\n")
  
  # Exit early in parallel mode
  quit(status = 0)
  
} else {
  # Sequential mode: run all simulations
  for (corr_level in correlation_levels) {
    corr_name <- correlation_names[which(correlation_levels == corr_level)]
    cat(sprintf("\nTesting %s correlation (%.1f) (%d simulations)...\n", 
                corr_name, corr_level, n_simulations_per_condition))
    
    for (sim in 1:n_simulations_per_condition) {
      if (sim %% 10 == 0) cat("  Simulation", sim, "of", n_simulations_per_condition, "\n")
      
      set.seed(42 + sim + (which(correlation_levels == corr_level) - 1) * 100)
      
      # Generate correlated data
      corr_data <- generate_correlated_dataset(
        n_samples = n_samples,
        n_taxa = n_taxa,
        effect_size = effect_size,
        correlation_level = corr_level
      )
      
      # CLR transformation
      X_clr <- log(corr_data$counts + 1)
      X_clr <- X_clr - rowMeans(X_clr)
      colnames(X_clr) <- paste0("Taxa_", 1:ncol(X_clr))
      
      # Run MeLSI
      melsi_result <- melsi(X_clr, corr_data$metadata$Group,
                           n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)
      
      # Calculate recovery metrics
      recovery_metrics <- calculate_recovery_metrics(
        melsi_result,
        corr_data$true_signal_taxa,
        X_clr,
        corr_data$metadata$Group
      )
      
      # Run traditional methods
      dist_euc <- dist(X_clr)
      perm_euc <- adonis2(dist_euc ~ corr_data$metadata$Group, permutations = 999)
      
      dist_bray <- vegdist(corr_data$counts, method = "bray")
      perm_bray <- adonis2(dist_bray ~ corr_data$metadata$Group, permutations = 999)
      
      dist_jaccard <- vegdist(corr_data$counts, method = "jaccard", binary = TRUE)
      perm_jaccard <- adonis2(dist_jaccard ~ corr_data$metadata$Group, permutations = 999)
      
      tree <- ape::rtree(ncol(corr_data$counts))
      tree$tip.label <- colnames(corr_data$counts)
      dist_wunifrac <- GUniFrac::GUniFrac(corr_data$counts, tree)$unifracs[,,"d_1"]
      perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ corr_data$metadata$Group, permutations = 999)
      
      dist_uunifrac <- GUniFrac::GUniFrac(corr_data$counts, tree)$unifracs[,,"d_UW"]
      perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ corr_data$metadata$Group, permutations = 999)
      
      traditional_F <- c(Euclidean = perm_euc$F[1], 
                         BrayCurtis = perm_bray$F[1],
                         Jaccard = perm_jaccard$F[1],
                         WeightedUniFrac = perm_wunifrac$F[1],
                         UnweightedUniFrac = perm_uunifrac$F[1])
      best_trad <- names(which.max(traditional_F))
      best_trad_F <- max(traditional_F)
      
      # Store results
      results <- rbind(results, data.frame(
        Correlation_Level = corr_name,
        Correlation_Value = corr_level,
        Sample_Size = n_samples,
        Simulation = sim,
        MeLSI_F = melsi_result$F_observed,
        MeLSI_p = melsi_result$p_value,
        MeLSI_significant = melsi_result$p_value < 0.05,
        Precision_5 = recovery_metrics$precision_5,
        Precision_10 = recovery_metrics$precision_10,
        Precision_20 = recovery_metrics$precision_20,
        Recall_5 = recovery_metrics$recall_5,
        Recall_10 = recovery_metrics$recall_10,
        Recall_20 = recovery_metrics$recall_20,
        Mean_Rank = recovery_metrics$mean_rank,
        AUC_ROC = recovery_metrics$auc_roc,
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

# ==============================================================================
# Summary Statistics (only in sequential mode)
# ==============================================================================
if (!parallel_mode && nrow(results) > 0) {
  cat("\n")
  cat("==============================================================================\n")
  cat("SUMMARY STATISTICS\n")
  cat("==============================================================================\n\n")
  
  # Summary by correlation level
  for (corr_name in correlation_names) {
    subset_results <- results[results$Correlation_Level == corr_name, ]
    
    if (nrow(subset_results) == 0) next
    
    cat(sprintf("\n%s Correlation:\n", corr_name))
    cat(sprintf("  MeLSI Power: %.1f%% (%d/%d)\n",
                mean(subset_results$MeLSI_significant) * 100,
                sum(subset_results$MeLSI_significant), nrow(subset_results)))
    cat(sprintf("  MeLSI Mean F: %.3f (SD = %.3f)\n",
                mean(subset_results$MeLSI_F), sd(subset_results$MeLSI_F)))
    cat(sprintf("  Precision@10: %.3f (SD = %.3f)\n",
                mean(subset_results$Precision_10, na.rm = TRUE), 
                sd(subset_results$Precision_10, na.rm = TRUE)))
    cat(sprintf("  Recall@10: %.3f (SD = %.3f)\n",
                mean(subset_results$Recall_10, na.rm = TRUE),
                sd(subset_results$Recall_10, na.rm = TRUE)))
    cat(sprintf("  Best Traditional: %s (Power: %.1f%%)\n",
                names(sort(table(subset_results$Best_Traditional), decreasing = TRUE))[1],
                mean(subset_results[[paste0(names(sort(table(subset_results$Best_Traditional), decreasing = TRUE))[1], "_significant")]], na.rm = TRUE) * 100))
  }
  
  # Save detailed results
  write.csv(results, "table6_feature_correlation_results.csv", row.names = FALSE)
  cat("\nDetailed results saved to: table6_feature_correlation_results.csv\n")
  cat("==============================================================================\n")
}
