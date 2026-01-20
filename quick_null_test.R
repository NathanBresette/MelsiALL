# Quick null hypothesis test
library(vegan)
source("melsi_robust.R")

# Create a simple null dataset
set.seed(123)
n_samples <- 100
n_taxa <- 200

# Create counts from same distribution (no group differences)
counts <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
                nrow = n_samples, ncol = n_taxa)

# Random group assignments (no real biological difference)
metadata <- data.frame(
  Group = sample(c("Group1", "Group2"), n_samples, replace = TRUE)
)

cat("Null dataset created:\n")
cat("Samples:", nrow(counts), "\n")
cat("Taxa:", ncol(counts), "\n")
cat("Groups randomly assigned (no real difference)\n")

# Test MeLSI
cat("\nTesting MeLSI on NULL data...\n")
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)

melsi_results <- run_melsi_permutation_test(
  counts_clr, metadata$Group,
  n_perms = 99,
  B = 50,
  m_frac = 0.7,
  show_progress = TRUE
)

cat("\nMeLSI Results on NULL data:\n")
cat("F-statistic:", round(melsi_results$F_observed, 4), "\n")
cat("P-value:", round(melsi_results$p_value, 4), "\n")
cat("Significant:", ifelse(melsi_results$p_value < 0.05, "YES (FALSE POSITIVE!)", "NO (CORRECT)"), "\n")

# Test traditional methods
cat("\nTesting traditional methods on NULL data...\n")

# Bray-Curtis
dist_bray <- vegdist(counts, method = "bray")
permanova_bray <- adonis2(dist_bray ~ Group, data = metadata, permutations = 99)

# Euclidean
dist_euclidean <- vegdist(counts_clr, method = "euclidean")
permanova_euclidean <- adonis2(dist_euclidean ~ Group, data = metadata, permutations = 99)

cat("\nTraditional Results on NULL data:\n")
cat("Bray-Curtis: p =", round(permanova_bray$`Pr(>F)`[1], 4), 
    ifelse(permanova_bray$`Pr(>F)`[1] < 0.05, "*** (FALSE POSITIVE)", ""), "\n")
cat("Euclidean: p =", round(permanova_euclidean$`Pr(>F)`[1], 4), 
    ifelse(permanova_euclidean$`Pr(>F)`[1] < 0.05, "*** (FALSE POSITIVE)", ""), "\n")

# Assessment
cat("\nASSESSMENT:\n")
if (melsi_results$p_value < 0.05) {
  cat("❌ PROBLEM: MeLSI found significance on NULL data!\n")
  cat("   This suggests inflated Type I error rates\n")
} else {
  cat("✅ GOOD: MeLSI correctly found no significance on NULL data\n")
}

if (permanova_bray$`Pr(>F)`[1] < 0.05) {
  cat("❌ Bray-Curtis also found false positive\n")
} else {
  cat("✅ Bray-Curtis correctly found no significance\n")
}
