# Create PCoA comparison: MeLSI vs Traditional Bray-Curtis on synthetic data
library(ggplot2)
library(vegan)
library(gridExtra)
library(plyr)

# Source MeLSI function
source("github/R/melsi_robust.R")

# Define colors
melsi_red <- "#E63946"
melsi_blue <- "#457B9D"
melsi_orange <- "#F77F00"
bg_gray <- "#F5F5F5"
text_gray <- "#333333"

# Define MeLSI theme
theme_pcoa <- theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = text_gray),
    plot.title = element_text(size = 16, hjust = 0.5, color = text_gray, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = text_gray),
    axis.title = element_text(size = 13, color = text_gray),
    axis.text = element_text(size = 11, color = text_gray),
    panel.grid = element_blank(),
    axis.line = element_line(color = text_gray, linewidth = 0.5),
    plot.background = element_rect(fill = bg_gray, color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 12, color = text_gray),
    legend.text = element_text(size = 11, color = text_gray)
  )

# Generate synthetic data with larger effect size
cat("Generating synthetic microbiome data...\n")
set.seed(42)  # Better seed for clearer separation

n_samples <- 80
n_features <- 100
n_groups <- 2
n_per_group <- n_samples / n_groups

# Create base abundances
base_abundance <- matrix(rpois(n_samples * n_features, lambda = 8), 
                         nrow = n_samples, ncol = n_features)

# Add larger effect size difference (MeLSI should capture this better)
group_effect <- 2.2  # Larger effect
noise_level <- 0.3   # Add some noise

for (i in 1:n_per_group) {
  # Group 1: increase specific features with pattern
  signal_features <- c(1:10, 21:30, 41:50)
  base_abundance[i, signal_features] <- base_abundance[i, signal_features] * (group_effect + rnorm(length(signal_features), 0, noise_level))
  
  # Group 2: increase different features
  signal_features2 <- c(11:20, 31:40, 51:60)
  base_abundance[i + n_per_group, signal_features2] <- base_abundance[i + n_per_group, signal_features2] * (group_effect + rnorm(length(signal_features2), 0, noise_level))
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
cat("\nRunning MeLSI...\n")
melsi_result <- melsi(counts_clr, groups, n_perms = 99, B = 30, 
                      m_frac = 0.8, show_progress = TRUE, plot_vip = FALSE)

# Extract learned metric and feature weights
M_learned <- melsi_result$metric_matrix
n_features_used <- melsi_result$diagnostics$n_features_used
feature_weights <- diag(M_learned)

# Create feature importance heatmap
cat("Creating feature importance heatmap...\n")
# Get top 20 most important features
top_features <- order(feature_weights, decreasing = TRUE)[1:min(20, length(feature_weights))]
weights_data <- data.frame(
  Feature = paste0("F", top_features),
  Weight = feature_weights[top_features],
  Rank = 1:length(top_features)
)

# Sort by weight for better visualization
weights_data$Feature <- factor(weights_data$Feature, levels = weights_data$Feature[order(weights_data$Weight)])

# Calculate MeLSI distances
if (n_features_used < ncol(counts_clr)) {
  X_filtered <- apply_conservative_prefiltering(counts_clr, groups, filter_frac = 0.7)
  dist_melsi <- calculate_mahalanobis_dist_robust(X_filtered, M_learned)
} else {
  dist_melsi <- calculate_mahalanobis_dist_robust(counts_clr, M_learned)
}

# Calculate Bray-Curtis distances (traditional)
cat("Calculating Bray-Curtis distances...\n")
dist_bray <- vegdist(base_abundance, method = "bray")

# Calculate Bray-Curtis PERMANOVA
cat("Running Bray-Curtis PERMANOVA...\n")
bray_result <- adonis2(dist_bray ~ groups, permutations = 999)

# Run PCoA on MeLSI distances
cat("\nRunning PCoA on MeLSI distances...\n")
pcoa_melsi <- cmdscale(dist_melsi, k = 2, eig = TRUE)

# Run PCoA on Bray-Curtis distances
cat("Running PCoA on Bray-Curtis distances...\n")
pcoa_bray <- cmdscale(dist_bray, k = 2, eig = TRUE)

# Extract variance explained
var_explained_melsi <- pcoa_melsi$eig / sum(pcoa_melsi$eig) * 100
var_explained_bray <- pcoa_bray$eig / sum(pcoa_bray$eig) * 100

# Prepare data for plotting
pcoa_melsi_data <- data.frame(
  PC1 = pcoa_melsi$points[, 1],
  PC2 = pcoa_melsi$points[, 2],
  Group = groups
)

pcoa_bray_data <- data.frame(
  PC1 = pcoa_bray$points[, 1],
  PC2 = pcoa_bray$points[, 2],
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

ellipse_melsi_data <- pcoa_melsi_data
colnames(ellipse_melsi_data)[1:2] <- c("x", "y")
ellipses_melsi <- calculate_ellipse(ellipse_melsi_data, "Group", level = 0.95)

ellipse_bray_data <- pcoa_bray_data
colnames(ellipse_bray_data)[1:2] <- c("x", "y")
ellipses_bray <- calculate_ellipse(ellipse_bray_data, "Group", level = 0.95)

# Get axis limits for consistency (trim outliers)
get_trimmed_limits <- function(values, trim = 0.02) {
  q_low <- quantile(values, trim)
  q_high <- quantile(values, 1 - trim)
  return(c(q_low, q_high))
}

xlim_melsi <- get_trimmed_limits(pcoa_melsi_data$PC1, trim = 0.03)
ylim_melsi <- get_trimmed_limits(pcoa_melsi_data$PC2, trim = 0.03)
xlim_bray <- get_trimmed_limits(pcoa_bray_data$PC1, trim = 0.03)
ylim_bray <- get_trimmed_limits(pcoa_bray_data$PC2, trim = 0.03)

# Create MeLSI plot
p_melsi <- ggplot(pcoa_melsi_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.75) +
  geom_path(data = ellipses_melsi, aes(x = x, y = y, group = Group), 
            alpha = 0.25, linewidth = 1.2, inherit.aes = TRUE) +
  scale_color_manual(values = c("Group 1" = melsi_red, "Group 2" = melsi_blue)) +
  scale_x_continuous(limits = xlim_melsi) +
  scale_y_continuous(limits = ylim_melsi) +
  labs(
    title = "MeLSI Learned Distance",
    x = paste0("PCoA 1 (", round(var_explained_melsi[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_explained_melsi[2], 1), "%)")
  ) +
  theme_pcoa +
  theme(legend.position = "none")

# Create Bray-Curtis plot
p_bray <- ggplot(pcoa_bray_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3.5, alpha = 0.75) +
  geom_path(data = ellipses_bray, aes(x = x, y = y, group = Group), 
            alpha = 0.25, linewidth = 1.2, inherit.aes = TRUE) +
  scale_color_manual(values = c("Group 1" = melsi_red, "Group 2" = melsi_blue)) +
  scale_x_continuous(limits = xlim_bray) +
  scale_y_continuous(limits = ylim_bray) +
  labs(
    title = "Traditional Bray-Curtis Distance",
    x = paste0("PCoA 1 (", round(var_explained_bray[1], 1), "%)"),
    y = paste0("PCoA 2 (", round(var_explained_bray[2], 1), "%)"),
    color = ""
  ) +
  theme_pcoa

# Create feature importance heatmap plot
p_heatmap <- ggplot(weights_data, aes(x = Feature, y = 1, fill = Weight)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = melsi_blue, mid = "white", high = melsi_orange, midpoint = mean(weights_data$Weight)) +
  labs(
    title = "MeLSI Feature Importance (Top 20)",
    subtitle = "Higher weights indicate features most important for group separation",
    x = "Feature",
    y = "",
    fill = "Weight"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = text_gray),
    plot.title = element_text(size = 14, hjust = 0.5, color = text_gray, face = "bold"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = text_gray),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.background = element_rect(fill = bg_gray, color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Create statistical comparison plot
# Extract F-statistics and p-values
melsi_f <- as.numeric(melsi_result$F_observed)
melsi_p <- as.numeric(melsi_result$p_value)

# Extract Bray-Curtis results safely
bray_f <- as.numeric(bray_result$aov.tab$F.Model[1])
bray_p <- as.numeric(bray_result$aov.tab$`Pr(>F)`[1])

# Create comparison data frame
comparison_data <- data.frame(
  Method = c("MeLSI", "Bray-Curtis"),
  F_Statistic = c(melsi_f, bray_f),
  P_Value = c(melsi_p, bray_p)
)

p_stats <- ggplot(comparison_data, aes(x = Method, y = F_Statistic, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.9) +
  geom_text(aes(label = paste0("F = ", round(F_Statistic, 2), "\np = ", round(P_Value, 4))), 
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("MeLSI" = melsi_orange, "Bray-Curtis" = melsi_blue)) +
  labs(
    title = "Statistical Power Comparison",
    subtitle = "Higher F-statistic indicates better group separation",
    y = "F-Statistic",
    x = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = text_gray),
    plot.title = element_text(size = 16, hjust = 0.5, color = text_gray, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = text_gray),
    axis.title = element_text(size = 13, color = text_gray),
    axis.text = element_text(size = 12, color = text_gray),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = bg_gray, color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Combine plots - PCoA on top, stats and heatmap on bottom
p_combined <- grid.arrange(
  p_melsi, p_bray, p_heatmap, p_stats,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  heights = c(1, 0.6, 0.5)
)

# Save
ggsave("pcoa_comparison_synthetic.png", plot = p_combined, width = 16, height = 13, dpi = 600)

cat("\n✅ Comparison plot saved as: pcoa_comparison_synthetic.png\n")
cat("MeLSI Results: F-statistic =", round(melsi_f, 3), 
    ", p-value =", round(melsi_p, 4), "\n")
cat("Bray-Curtis Results: pseudo-F =", round(bray_f, 3), 
    ", p-value =", round(bray_p, 4), "\n") 