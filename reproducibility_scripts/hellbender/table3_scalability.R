#!/usr/bin/env Rscript
# ==============================================================================
# Table 3: Scalability Across Sample Size and Dimensionality
# ==============================================================================
# Rigorous version with repeated simulations for variance estimation
# ==============================================================================

# Load required packages
library(vegan)
library(MeLSI)
library(GUniFrac)
library(ape)

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
n_simulations_per_condition <- 10  # Rigorous: 10 simulations for variance estimation

# Conditions to test
# Varying n (fixed p=200): n = 20, 50, 100, 200, 500
# Varying p (fixed n=100): p = 50, 100, 200, 500, 1000
sample_sizes <- c(20, 50, 100, 200, 500)
taxa_sizes <- c(50, 100, 200, 500, 1000)

# ==============================================================================
# Helper Function: Generate Scalability Test Data
# ==============================================================================
generate_scalability_data <- function(n_samples, n_taxa) {
  counts <- matrix(rnbinom(n_samples * n_taxa, mu = 30, size = 0.8), 
                  nrow = n_samples, ncol = n_taxa)
  counts[counts < 3] <- 0
  
  # Add medium effect signal
  group1_size <- n_samples %/% 2
  signal_taxa <- min(10, n_taxa %/% 10)
  signal_indices <- sample(1:n_taxa, signal_taxa)
  
  for (i in 1:group1_size) {
    counts[i, signal_indices] <- counts[i, signal_indices] * 2
  }
  
  group <- c(rep("Group1", group1_size), rep("Group2", n_samples - group1_size))
  
  return(list(counts = counts, group = group))
}

# ==============================================================================
# Run Scalability Tests - REPEATED SIMULATIONS
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("REPRODUCING TABLE 3: Scalability Analysis (Rigorous)\n")
cat("Running multiple simulations per condition for variance estimation\n")
cat("==============================================================================\n\n")

results <- data.frame()

# Determine which simulation to run (for parallel execution)
if (parallel_mode) {
  # Calculate which condition this simulation belongs to
  # Conditions: varying n (5 conditions) + varying p (5 conditions) = 10 total
  # Each condition has 10 simulations
  conditions <- expand.grid(
    condition_type = c("vary_n", "vary_p"),
    condition_value = c(sample_sizes, taxa_sizes),
    sim = 1:n_simulations_per_condition
  )
  
  # Filter to valid conditions
  valid_conditions <- data.frame()
  for (ct in c("vary_n", "vary_p")) {
    if (ct == "vary_n") {
      for (n_val in sample_sizes) {
        for (sim in 1:n_simulations_per_condition) {
          valid_conditions <- rbind(valid_conditions, data.frame(
            condition_type = ct,
            condition_value = n_val,
            sim = sim
          ))
        }
      }
    } else {
      for (p_val in taxa_sizes) {
        for (sim in 1:n_simulations_per_condition) {
          valid_conditions <- rbind(valid_conditions, data.frame(
            condition_type = ct,
            condition_value = p_val,
            sim = sim
          ))
        }
      }
    }
  }
  
  total_conditions <- nrow(valid_conditions)
  
  if (sim_index < 1 || sim_index > total_conditions) {
    stop("sim_index out of range")
  }
  
  current_condition <- valid_conditions[sim_index, ]
  cond_type <- as.character(current_condition$condition_type)
  cond_val <- as.integer(current_condition$condition_value)
  sim <- as.integer(current_condition$sim)
  
  cat(sprintf("Running simulation %d/%d: %s, value=%d\n", 
              sim_index, total_conditions, cond_type, cond_val))
  
  # Set unique seed for this simulation
  set.seed(42 + sim_index)
  
  # Generate data based on condition type
  if (cond_type == "vary_n") {
    n_samples <- cond_val
    n_taxa <- 200
  } else {
    n_samples <- 100
    n_taxa <- cond_val
  }
  
  data <- generate_scalability_data(n_samples, n_taxa)
  
  # Set column names for counts matrix (needed for UniFrac)
  colnames(data$counts) <- paste0("Taxa_", 1:ncol(data$counts))
  
  # CLR transformation
  X_clr <- log(data$counts + 1)
  X_clr <- X_clr - rowMeans(X_clr)
  colnames(X_clr) <- colnames(data$counts)  # Use same names
  
  # Time MeLSI
  start_time <- Sys.time()
  melsi_result <- melsi(X_clr, data$group, n_perms = 200, B = 30, 
                       show_progress = FALSE, plot_vip = FALSE)
  melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Euclidean PERMANOVA
  start_time <- Sys.time()
  dist_euc <- dist(X_clr)
  perm_euc <- adonis2(dist_euc ~ data$group, permutations = 999)
  euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Bray-Curtis PERMANOVA
  start_time <- Sys.time()
  dist_bray <- vegdist(data$counts, method = "bray")
  perm_bray <- adonis2(dist_bray ~ data$group, permutations = 999)
  bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Jaccard PERMANOVA
  start_time <- Sys.time()
  dist_jaccard <- vegdist(data$counts, method = "jaccard", binary = TRUE)
  perm_jaccard <- adonis2(dist_jaccard ~ data$group, permutations = 999)
  jaccard_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Weighted UniFrac PERMANOVA (using random tree for synthetic data)
  start_time <- Sys.time()
  tree <- ape::rtree(ncol(data$counts))
  tree$tip.label <- colnames(data$counts)
  dist_wunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_1"]
  perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ data$group, permutations = 999)
  wunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Time Unweighted UniFrac PERMANOVA
  start_time <- Sys.time()
  dist_uunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_UW"]
  perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ data$group, permutations = 999)
  uunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Find best traditional method
  traditional_F <- c(Euclidean = perm_euc$F[1], 
                     BrayCurtis = perm_bray$F[1],
                     Jaccard = perm_jaccard$F[1],
                     WeightedUniFrac = perm_wunifrac$F[1],
                     UnweightedUniFrac = perm_uunifrac$F[1])
  best_trad <- names(which.max(traditional_F))
  best_trad_F <- max(traditional_F)
  best_trad_time <- switch(best_trad,
                          "Euclidean" = euc_time,
                          "BrayCurtis" = bray_time,
                          "Jaccard" = jaccard_time,
                          "WeightedUniFrac" = wunifrac_time,
                          "UnweightedUniFrac" = uunifrac_time)
  
  # Store single result
  result <- data.frame(
    Condition_Type = cond_type,
    Condition_Value = cond_val,
    n_samples = n_samples,
    n_taxa = n_taxa,
    Simulation = sim,
    MeLSI_time_sec = melsi_time,
    MeLSI_F = melsi_result$F_observed,
    MeLSI_p = melsi_result$p_value,
    Euclidean_F = perm_euc$F[1],
    Euclidean_time_sec = euc_time,
    BrayCurtis_F = perm_bray$F[1],
    BrayCurtis_time_sec = bray_time,
    Jaccard_F = perm_jaccard$F[1],
    Jaccard_time_sec = jaccard_time,
    WeightedUniFrac_F = perm_wunifrac$F[1],
    WeightedUniFrac_time_sec = wunifrac_time,
    UnweightedUniFrac_F = perm_uunifrac$F[1],
    UnweightedUniFrac_time_sec = uunifrac_time,
    Best_Traditional = best_trad,
    Best_Traditional_F = best_trad_F,
    Best_Traditional_time_sec = best_trad_time
  )
  
  # Save single result to file (for parallel collection)
  output_file <- sprintf("table3_sim_%d.csv", sim_index)
  write.csv(result, output_file, row.names = FALSE)
  cat("Result saved to:", output_file, "\n")
  
  # Exit early in parallel mode - don't calculate summary statistics
  quit(status = 0)
  
} else {
  # Sequential mode: run all simulations
  # Varying sample sizes
  for (n in sample_sizes) {
    cat(sprintf("\nTesting n = %d samples (fixed p=200)...\n", n))
    for (sim in 1:n_simulations_per_condition) {
      if (sim %% 5 == 0) cat("  Simulation", sim, "of", n_simulations_per_condition, "\n")
      
      set.seed(42 + sim + (which(sample_sizes == n) - 1) * 100)
      
      data <- generate_scalability_data(n, 200)
      
      # Set column names for counts matrix (needed for UniFrac)
      colnames(data$counts) <- paste0("Taxa_", 1:ncol(data$counts))
      
      # CLR transformation
      X_clr <- log(data$counts + 1)
      X_clr <- X_clr - rowMeans(X_clr)
      colnames(X_clr) <- colnames(data$counts)  # Use same names
      
      # Time MeLSI
      start_time <- Sys.time()
      melsi_result <- melsi(X_clr, data$group, n_perms = 200, B = 30, 
                           show_progress = FALSE, plot_vip = FALSE)
      melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Euclidean PERMANOVA
      start_time <- Sys.time()
      dist_euc <- dist(X_clr)
      perm_euc <- adonis2(dist_euc ~ data$group, permutations = 999)
      euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Bray-Curtis PERMANOVA
      start_time <- Sys.time()
      dist_bray <- vegdist(data$counts, method = "bray")
      perm_bray <- adonis2(dist_bray ~ data$group, permutations = 999)
      bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Jaccard PERMANOVA
      start_time <- Sys.time()
      dist_jaccard <- vegdist(data$counts, method = "jaccard", binary = TRUE)
      perm_jaccard <- adonis2(dist_jaccard ~ data$group, permutations = 999)
      jaccard_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Weighted UniFrac PERMANOVA (using random tree for synthetic data)
      start_time <- Sys.time()
      tree <- ape::rtree(ncol(data$counts))
      tree$tip.label <- colnames(data$counts)
      dist_wunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_1"]
      perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ data$group, permutations = 999)
      wunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Unweighted UniFrac PERMANOVA
      start_time <- Sys.time()
      dist_uunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_UW"]
      perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ data$group, permutations = 999)
      uunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Find best traditional method
      traditional_F <- c(Euclidean = perm_euc$F[1], 
                         BrayCurtis = perm_bray$F[1],
                         Jaccard = perm_jaccard$F[1],
                         WeightedUniFrac = perm_wunifrac$F[1],
                         UnweightedUniFrac = perm_uunifrac$F[1])
      best_trad <- names(which.max(traditional_F))
      best_trad_F <- max(traditional_F)
      best_trad_time <- switch(best_trad,
                              "Euclidean" = euc_time,
                              "BrayCurtis" = bray_time,
                              "Jaccard" = jaccard_time,
                              "WeightedUniFrac" = wunifrac_time,
                              "UnweightedUniFrac" = uunifrac_time)
      
      results <- rbind(results, data.frame(
        Condition_Type = "vary_n",
        Condition_Value = n,
        n_samples = n,
        n_taxa = 200,
        Simulation = sim,
        MeLSI_time_sec = melsi_time,
        MeLSI_F = melsi_result$F_observed,
        MeLSI_p = melsi_result$p_value,
        Euclidean_F = perm_euc$F[1],
        Euclidean_time_sec = euc_time,
        BrayCurtis_F = perm_bray$F[1],
        BrayCurtis_time_sec = bray_time,
        Jaccard_F = perm_jaccard$F[1],
        Jaccard_time_sec = jaccard_time,
        WeightedUniFrac_F = perm_wunifrac$F[1],
        WeightedUniFrac_time_sec = wunifrac_time,
        UnweightedUniFrac_F = perm_uunifrac$F[1],
        UnweightedUniFrac_time_sec = uunifrac_time,
        Best_Traditional = best_trad,
        Best_Traditional_F = best_trad_F,
        Best_Traditional_time_sec = best_trad_time
      ))
    }
  }
  
  # Varying taxa sizes
  for (p in taxa_sizes) {
    cat(sprintf("\nTesting p = %d taxa (fixed n=100)...\n", p))
    for (sim in 1:n_simulations_per_condition) {
      if (sim %% 5 == 0) cat("  Simulation", sim, "of", n_simulations_per_condition, "\n")
      
      set.seed(42 + sim + (which(taxa_sizes == p) - 1) * 100 + 1000)
      
      data <- generate_scalability_data(100, p)
      
      # Set column names for counts matrix (needed for UniFrac)
      colnames(data$counts) <- paste0("Taxa_", 1:ncol(data$counts))
      
      # CLR transformation
      X_clr <- log(data$counts + 1)
      X_clr <- X_clr - rowMeans(X_clr)
      colnames(X_clr) <- colnames(data$counts)  # Use same names
      
      # Time MeLSI
      start_time <- Sys.time()
      melsi_result <- melsi(X_clr, data$group, n_perms = 200, B = 30, 
                           show_progress = FALSE, plot_vip = FALSE)
      melsi_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Euclidean PERMANOVA
      start_time <- Sys.time()
      dist_euc <- dist(X_clr)
      perm_euc <- adonis2(dist_euc ~ data$group, permutations = 999)
      euc_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Bray-Curtis PERMANOVA
      start_time <- Sys.time()
      dist_bray <- vegdist(data$counts, method = "bray")
      perm_bray <- adonis2(dist_bray ~ data$group, permutations = 999)
      bray_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Jaccard PERMANOVA
      start_time <- Sys.time()
      dist_jaccard <- vegdist(data$counts, method = "jaccard", binary = TRUE)
      perm_jaccard <- adonis2(dist_jaccard ~ data$group, permutations = 999)
      jaccard_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Weighted UniFrac PERMANOVA (using random tree for synthetic data)
      start_time <- Sys.time()
      tree <- ape::rtree(ncol(data$counts))
      tree$tip.label <- colnames(data$counts)
      dist_wunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_1"]
      perm_wunifrac <- adonis2(as.dist(dist_wunifrac) ~ data$group, permutations = 999)
      wunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Time Unweighted UniFrac PERMANOVA
      start_time <- Sys.time()
      dist_uunifrac <- GUniFrac::GUniFrac(data$counts, tree)$unifracs[,,"d_UW"]
      perm_uunifrac <- adonis2(as.dist(dist_uunifrac) ~ data$group, permutations = 999)
      uunifrac_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      # Find best traditional method
      traditional_F <- c(Euclidean = perm_euc$F[1], 
                         BrayCurtis = perm_bray$F[1],
                         Jaccard = perm_jaccard$F[1],
                         WeightedUniFrac = perm_wunifrac$F[1],
                         UnweightedUniFrac = perm_uunifrac$F[1])
      best_trad <- names(which.max(traditional_F))
      best_trad_F <- max(traditional_F)
      best_trad_time <- switch(best_trad,
                              "Euclidean" = euc_time,
                              "BrayCurtis" = bray_time,
                              "Jaccard" = jaccard_time,
                              "WeightedUniFrac" = wunifrac_time,
                              "UnweightedUniFrac" = uunifrac_time)
      
      results <- rbind(results, data.frame(
        Condition_Type = "vary_p",
        Condition_Value = p,
        n_samples = 100,
        n_taxa = p,
        Simulation = sim,
        MeLSI_time_sec = melsi_time,
        MeLSI_F = melsi_result$F_observed,
        MeLSI_p = melsi_result$p_value,
        Euclidean_F = perm_euc$F[1],
        Euclidean_time_sec = euc_time,
        BrayCurtis_F = perm_bray$F[1],
        BrayCurtis_time_sec = bray_time,
        Jaccard_F = perm_jaccard$F[1],
        Jaccard_time_sec = jaccard_time,
        WeightedUniFrac_F = perm_wunifrac$F[1],
        WeightedUniFrac_time_sec = wunifrac_time,
        UnweightedUniFrac_F = perm_uunifrac$F[1],
        UnweightedUniFrac_time_sec = uunifrac_time,
        Best_Traditional = best_trad,
        Best_Traditional_F = best_trad_F,
        Best_Traditional_time_sec = best_trad_time
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
  
  # Summary by condition
  for (cond_type in c("vary_n", "vary_p")) {
    subset_results <- results[results$Condition_Type == cond_type, ]
    
    if (nrow(subset_results) == 0) next
    
    cat(sprintf("\n%s:\n", ifelse(cond_type == "vary_n", "Varying Sample Size (n)", "Varying Dimensionality (p)")))
    
    for (cond_val in unique(subset_results$Condition_Value)) {
      cond_subset <- subset_results[subset_results$Condition_Value == cond_val, ]
      
      cat(sprintf("  Value = %d:\n", cond_val))
      cat(sprintf("    MeLSI: Mean F = %.3f (SD = %.3f), Mean Time = %.1fs (SD = %.1fs)\n",
                  mean(cond_subset$MeLSI_F), sd(cond_subset$MeLSI_F),
                  mean(cond_subset$MeLSI_time_sec), sd(cond_subset$MeLSI_time_sec)))
      cat(sprintf("    Best Traditional: Mean F = %.3f (SD = %.3f), Mean Time = %.1fs (SD = %.1fs)\n",
                  mean(cond_subset$Best_Traditional_F), sd(cond_subset$Best_Traditional_F),
                  mean(cond_subset$Best_Traditional_time_sec), sd(cond_subset$Best_Traditional_time_sec)))
    }
  }
  
  # Save detailed results
  write.csv(results, "table3_scalability_results.csv", row.names = FALSE)
  cat("\nDetailed results saved to: table3_scalability_results.csv\n")
  cat("==============================================================================\n")
}
