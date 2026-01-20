# Create Feature Importance Heatmap for MeLSI
library(ggplot2)
library(RColorBrewer)

# Source MeLSI function
source("github/R/melsi_robust.R")

# Define colors
melsi_red <- "#E63946"
melsi_blue <- "#457B9D"
melsi_orange <- "#F77F00"
bg_gray <- "#F5F5F5"
text_gray <- "#333333"

# Generate synthetic data
cat("Generating synthetic data...\n")
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

# Add meaningful feature names (taxa names for real microbiome data)
# In practice, these would be actual bacterial names
feature_names <- paste0("Taxon_", sprintf("%02d", 1:n_features))

# Run MeLSI
cat("Running MeLSI...\n")
melsi_result <- melsi(counts_clr, groups, n_perms = 99, B = 30, 
                      m_frac = 0.8, show_progress = TRUE, plot_vip = FALSE)

# Extract feature weights
M_learned <- melsi_result$metric_matrix
feature_weights <- diag(M_learned)

# Get top features
top_features <- order(feature_weights, decreasing = TRUE)[1:30]  # Top 30

# Create data frame for plotting
weights_data <- data.frame(
  Taxon = feature_names[top_features],
  Weight = feature_weights[top_features],
  Rank = 1:length(top_features)
)

# Sort by weight for better visualization
weights_data$Taxon <- factor(weights_data$Taxon, 
                             levels = weights_data$Taxon[order(weights_data$Weight)])

# Create the heatmap
p <- ggplot(weights_data, aes(x = Taxon, y = 1, fill = Weight)) +
  geom_tile(color = "white", linewidth = 0.8, height = 0.9) +
  scale_fill_gradient2(
    low = melsi_blue, 
    mid = "white", 
    high = melsi_orange, 
    midpoint = median(weights_data$Weight),
    name = "Weight",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 15,
      barheight = 0.8
    )
  ) +
  labs(
    title = "MeLSI Feature Importance",
    subtitle = "Top 30 taxa driving group separation (ordered by importance)",
    x = "Taxa",
    y = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = text_gray),
    plot.title = element_text(size = 18, hjust = 0.5, color = text_gray, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = text_gray, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = text_gray),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14, color = text_gray, face = "bold"),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 13, color = text_gray, face = "bold"),
    legend.text = element_text(size = 11, color = text_gray),
    plot.background = element_rect(fill = bg_gray, color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Save
ggsave("feature_importance_heatmap.png", plot = p, width = 14, height = 4, dpi = 600)

cat("\n✅ Feature importance heatmap saved as: feature_importance_heatmap.png\n")
cat("Top 5 most important features:\n")
for (i in 1:5) {
  cat(sprintf("%d. %s (weight: %.4f)\n", i, 
             weights_data$Taxon[order(weights_data$Weight, decreasing = TRUE)[i]],
             sort(weights_data$Weight, decreasing = TRUE)[i]))
}









