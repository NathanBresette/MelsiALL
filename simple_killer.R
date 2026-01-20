# Simple Killer Experiment
library(vegan)
source("melsi_robust.R")

# Create simple test data
set.seed(123)
n_samples <- 100
n_taxa <- 200
n_signal_taxa <- 30

# Create base counts
counts <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
                nrow = n_samples, ncol = n_taxa)

# Add 10% increase to first 30 taxa in second half
group2_start <- n_samples/2 + 1
counts[group2_start:n_samples, 1:n_signal_taxa] <- 
  counts[group2_start:n_samples, 1:n_signal_taxa] * 1.1

# Create metadata
metadata <- data.frame(
  Group = c(rep("Group1", n_samples/2), rep("Group2", n_samples/2))
)

cat("Dataset created:\n")
cat("Samples:", nrow(counts), "\n")
cat("Taxa:", ncol(counts), "\n")
cat("Signal taxa:", n_signal_taxa, "\n")

# Test DA analysis
cat("\nRunning DA analysis...\n")
p_values <- numeric(ncol(counts))
for (i in 1:ncol(counts)) {
  group1 <- counts[metadata$Group == "Group1", i]
  group2 <- counts[metadata$Group == "Group2", i]
  p_values[i] <- t.test(group1, group2)$p.value
}

# Multiple testing correction
p_adjusted <- p.adjust(p_values, method = "bonferroni")
n_significant_raw <- sum(p_values < 0.05)
n_significant_adjusted <- sum(p_adjusted < 0.05)

cat("DA Results:\n")
cat("Raw significant taxa:", n_significant_raw, "\n")
cat("Adjusted significant taxa:", n_significant_adjusted, "\n")

# Test MeLSI
cat("\nRunning MeLSI...\n")
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)

melsi_results <- run_melsi_permutation_test(
  counts_clr, metadata$Group,
  n_perms = 99,
  B = 50,
  m_frac = 0.7,
  show_progress = TRUE
)

cat("MeLSI Results:\n")
cat("F-statistic:", round(melsi_results$F_observed, 4), "\n")
cat("P-value:", round(melsi_results$p_value, 4), "\n")
cat("Significant:", ifelse(melsi_results$p_value < 0.05, "YES", "NO"), "\n")

# Test traditional methods
cat("\nRunning traditional methods...\n")
dist_bray <- vegdist(counts, method = "bray")
permanova_bray <- adonis2(dist_bray ~ Group, data = metadata, permutations = 99)

dist_euclidean <- vegdist(counts_clr, method = "euclidean")
permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = metadata, permutations = 99)

cat("Traditional Results:\n")
cat("Bray-Curtis F:", round(permanova_bray$F[1], 4), 
    "p:", round(permanova_bray$`Pr(>F)`[1], 4), "\n")
cat("Euclidean F:", round(permanova_euclidean$F[1], 4), 
    "p:", round(permanova_euclidean$`Pr(>F)`[1], 4), "\n")

# Summary
cat("\nEXPERIMENT SUMMARY:\n")
cat("DA (Bonferroni):", n_significant_adjusted, "significant taxa\n")
cat("MeLSI:", ifelse(melsi_results$p_value < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("Bray-Curtis:", ifelse(permanova_bray$`Pr(>F)`[1] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("Euclidean:", ifelse(permanova_euclidean$`Pr(>F)`[1] < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")

success <- (n_significant_adjusted == 0) && (melsi_results$p_value < 0.05)
if (success) {
  cat("\n🎉 KILLER EXPERIMENT SUCCEEDED!\n")
} else {
  cat("\n❌ KILLER EXPERIMENT FAILED\n")
}
