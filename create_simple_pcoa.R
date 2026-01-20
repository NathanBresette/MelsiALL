# Create Simple PCoA Plot with MeLSI Distances
library(ggplot2)
library(vegan)
library(plyr)

# Source MeLSI function
source("github/R/melsi_robust.R")

# Define colors
melsi_red <- "#E63946"
melsi_blue <- "#457B9D"
melsi_orange <- "#F77F00"
bg_gray <- "#F5F5F5"
text_gray <- "#333333"

# Generate synthetic data
cat("Generating synthetic microbiome data...\n")
set.seed(42)
n_samples <- 80
n_features <- 100
n_groups <- 2
n_per_group <- n_samples / n_groups

# Create base abundances
base_abundance <- matrix(rpois(n_samples * n_features, lambda = 8), 
                         nrow = n_samples, ncol = n_features)

# Add effect
group_effect <- 2.2
noise_level <- 0.3

for (i in 1:n_per_group) {
  signal_features <- c(1:10, 21:30, 41:50)
  base_abundance[i, signal_features] <- base_abundance[i, signal_features] * 
    (group_effect + rnorm(length(signal_features), 0, noise_level))
  
  signal_features2 <- c(11:20, 31:40, 51:60)
  base_abundance[i + n_per_group, signal_features2] <- 
    base_abundance[i + n_per_group, signal_features2] * 
    (group_effect + rnorm(length(signal_features2), 0, noise_level))
}

# Create group labels
groups <- factor(rep(c("Group 1", "Group 2"), each = n_per_group))

# Apply CLR transformation
counts_clr <- apply(base_abundance, 1, function(x) {
  log(x + 1e-10) - mean(log(x + 1e-10))
})
counts_clr <- t(counts_clr)

cat("Working with", nrow(counts_clr), "samples and", ncol(counts_clr), "features\n")

# Run MeLSI
cat("Running MeLSI...\n")
melsi_result <- melsi(counts_clr, groups, n_perms = 99, B = 30, 
                      m_frac = 0.8, show_progress = TRUE, plot_vip = FALSE)

# Extract learned metric
M_learned <- melsi_result$metric_matrix
n_features_used <- melsi_result$diagnostics$n_features_used

# Calculate MeLSI distances
if (n_features_used < ncol(counts_clr)) {
  X_filtered <- apply_conservative_prefiltering(counts_clr, groups, filter_frac = 0.7)
  dist_melsi <- calculate_mahalanobis_dist_robust(X_filtered, M_learned)
} else {
  dist_melsi <- calculate_mahalanobis_dist_robust(counts_clr, M_learned)
}

# Run PCoA
cat("Running PCoA...\n")
pcoa_result <- cmdscale(dist_melsi, k = 2, eig = TRUE)

# Extract variance explained
var_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100

# Prepare data for plotting
pcoa_data <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Group = groups
)

# Calculate ellipses
calculate_ellipse <- function(data, groups, level = 0.95) {
  results <- ddply(data, groups, function(df) {
    if (nrow(df) < 3) return(NULL)
    cov_mat <- cov(df[, c("x", "y")])
    mean_vec <- colMeans(df[, c("x", "y")])
    theta <- seq(0, 2*pi, length.out = 100)
    ellipse_points <- cbind(cos(theta), sin(theta))
    ellipse_points <- ellipse_points %*% chol(cov_mat * qchisq(level, 2))
    ellipse_points <- sweep(ellipse_points, 2, mean_vec, "+")
    return(data.frame(x = ellipse_points[, 1], y = ellipse_points[, 2]))
  })
  return(results)
}

ellipse_data <- pcoa_data
colnames(ellipse_data)[1:2] <- c("x", "y")
ellipses <- calculate_ellipse(ellipse_data, "Group", level = 0.95)

# Create PCoA plot
p <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_path(data = ellipses, aes(x = x, y = y, group = Group), 
            alpha = 0.25, linewidth = 1.3, inherit.aes = TRUE) +
  scale_color_manual(values = c("Group 1" = melsi_red, "Group 2" = melsi_blue)) +
  labs(
    title = "PCoA Visualization with MeLSI Learned Distances",
    subtitle = "Adaptive metric learning improves group separation",
    x = paste0("PCoA 1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PCoA 2 (", round(var_explained[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = text_gray),
    plot.title = element_text(size = 18, hjust = 0.5, color = text_gray, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = text_gray, margin = margin(b = 15)),
    axis.title = element_text(size = 14, color = text_gray, face = "bold"),
    axis.text = element_text(size = 12, color = text_gray),
    panel.grid = element_blank(),
    axis.line = element_line(color = text_gray, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 13, color = text_gray, face = "bold"),
    legend.text = element_text(size = 12, color = text_gray),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Save
ggsave("pcoa_plot.png", plot = p, width = 12, height = 8, dpi = 600)

cat("\n✅ PCoA plot saved as: pcoa_plot.png\n")
cat("F-statistic =", round(as.numeric(melsi_result$F_observed), 3), 
    ", p-value =", round(as.numeric(melsi_result$p_value), 4), "\n")
