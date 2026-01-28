#!/usr/bin/env Rscript
# Regenerate DietSwap figures with visible points

suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(vegan)
  library(gridExtra)
  library(grid)
  library(MeLSI)
})

set.seed(42)

# Load DietSwap dataset
cat("Loading DietSwap dataset...\n")
data(dietswap)

# Subset to baseline samples (timepoint.within.group == 1) and two groups (DI, HE)
dietswap_subset <- subset_samples(dietswap, timepoint.within.group == 1)
dietswap_subset <- subset_samples(dietswap_subset, group %in% c("DI", "HE"))
dietswap_subset <- prune_taxa(taxa_sums(dietswap_subset) > 0, dietswap_subset)

# Extract counts and metadata
counts <- as(otu_table(dietswap_subset), "matrix")
if (taxa_are_rows(dietswap_subset)) {
  counts <- t(counts)
}
metadata <- data.frame(sample_data(dietswap_subset))
metadata$group <- droplevels(metadata$group)
counts <- counts[rownames(metadata), ]

# Get group labels
groups <- metadata$group

# Remove samples with missing labels
valid_idx <- !is.na(groups)
counts <- counts[valid_idx, , drop = FALSE]
groups <- groups[valid_idx]

cat("Samples:", length(groups), "\n")
cat("Taxa:", ncol(counts), "\n")
cat("Groups:", paste(unique(groups), collapse = ", "), "\n")

# CLR transformation
cat("Applying CLR transformation...\n")
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)
colnames(counts_clr) <- colnames(counts)

# Run MeLSI
cat("Running MeLSI...\n")
melsi_res <- melsi(counts_clr, groups, n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)

# Group colors
group_colors <- c("DI" = "#1f77b4", "HE" = "#ff7f0e")
names(group_colors) <- c("DI", "HE")

# VIP Plot
cat("Creating VIP plot...\n")
vip_df <- data.frame(
  Feature = names(melsi_res$feature_weights),
  Weight = as.numeric(melsi_res$feature_weights)
) %>%
  arrange(desc(Weight)) %>%
  slice_head(n = 15) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

# VIP plot without directionality
vip_plot_no_dir <- ggplot(vip_df, aes(x = Weight, y = Feature)) +
  geom_col(fill = "#1f77b4") +
  geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 3) +
  labs(title = "Top DietSwap Taxa Weighted by MeLSI", x = "Feature Weight", y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

# VIP plot with directionality
has_directionality <- !is.null(melsi_res$directionality) && length(melsi_res$directionality) > 0
if (has_directionality) {
  vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$Feature)]
  vip_df$Directionality <- ifelse(is.na(vip_df$Directionality), "Unknown", vip_df$Directionality)
  
  unique_directions <- unique(vip_df$Directionality[vip_df$Directionality != "Unknown"])
  color_map <- setNames(group_colors[1:length(unique_directions)], unique_directions)
  if (any(vip_df$Directionality == "Unknown")) {
    color_map["Unknown"] <- "#95A5A6"
  }
  
  vip_plot_with_dir <- ggplot(vip_df, aes(x = Weight, y = Feature, fill = Directionality)) +
    geom_col() +
    scale_fill_manual(values = color_map, name = "Directionality") +
    geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 3) +
    labs(title = "Top DietSwap Taxa Weighted by MeLSI", x = "Feature Weight", y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
} else {
  vip_plot_with_dir <- vip_plot_no_dir
}

# Combined VIP plot
combined_vip <- arrangeGrob(vip_plot_no_dir, vip_plot_with_dir, ncol = 2, widths = c(1, 1.15))
ggsave("manuscript/figures/dietswap_vip_combined.png", combined_vip, width = 12, height = 4, dpi = 300)
cat("VIP plot saved\n")

# PCoA Plot
cat("Creating PCoA plot...\n")
dist_melsi <- MeLSI:::calculate_mahalanobis_dist_robust(counts_clr, melsi_res$metric_matrix)
pcoa_res <- ape::pcoa(dist_melsi)

explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
pcoa_df <- data.frame(
  PCoA1 = pcoa_res$vectors[, 1],
  PCoA2 = pcoa_res$vectors[, 2],
  Group = groups
)

# Create PCoA plot with LARGE, VISIBLE points
pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 5, alpha = 1.0, stroke = 0.8) +  # Large, solid, visible points
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "DietSwap PCoA with MeLSI Distance",
    x = paste0("PCoA1 (", explained[1], "%)"),
    y = paste0("PCoA2 (", explained[2], "%)"),
    color = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("manuscript/figures/dietswap_pcoa.png", pcoa_plot, width = 6, height = 4, dpi = 300)
cat("PCoA plot saved\n")

cat("\nAll figures regenerated successfully!\n")
