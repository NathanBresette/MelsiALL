#!/usr/bin/env Rscript
# Create SKIOME combined VIP + PCoA figure with perfect alignment
# Uses backup files to regenerate both plots together

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# Source MeLSI functions
source("../../github/R/melsi_robust.R")

# Load backup files
cat("Loading backup files...\n")
load("skiome_omnibus_results_backup.RData")

cat("Loaded backup files:\n")
cat("  - melsi_omnibus results object\n")
cat("  - X_clr: CLR-transformed data\n")
cat("  - groups: Group labels\n\n")

# Calculate directionality for 3 groups
cat("Calculating directionality for 3 groups...\n")
groups_unique <- unique(groups)
cat("  Groups found:", paste(groups_unique, collapse = ", "), "\n")

# Calculate mean abundance for each group
mean_by_group <- list()
for (g in groups_unique) {
  group_idx <- which(groups == g)
  mean_by_group[[as.character(g)]] <- colMeans(X_clr[group_idx, , drop = FALSE])
}

mean_matrix <- do.call(rbind, mean_by_group)
rownames(mean_matrix) <- groups_unique

directionality_calc <- apply(mean_matrix, 2, function(x) {
  max_idx <- which.max(x)
  return(groups_unique[max_idx])
})
names(directionality_calc) <- colnames(X_clr)

# Update results object
melsi_omnibus$directionality <- directionality_calc

# Clean taxon names
extract_taxon_name <- function(full_name) {
  name <- gsub("^sk__Bacteria;", "", full_name)
  if (grepl("g__[^;]+", name)) {
    genus <- gsub(".*g__([^;]+).*", "\\1", name)
    if (genus != "" && genus != name) return(genus)
  }
  if (grepl("f__[^;]+", name)) {
    family <- gsub(".*f__([^;]+).*", "\\1", name)
    if (family != "" && family != name) return(paste0(family, " (family)"))
  }
  if (grepl("o__[^;]+", name)) {
    order <- gsub(".*o__([^;]+).*", "\\1", name)
    if (order != "" && order != name) return(paste0(order, " (order)"))
  }
  return(gsub("^k__;", "", name))
}

cleaned_names <- sapply(names(melsi_omnibus$feature_weights), extract_taxon_name)
names(melsi_omnibus$feature_weights) <- cleaned_names
cleaned_dir_names <- sapply(names(melsi_omnibus$directionality), extract_taxon_name)
names(melsi_omnibus$directionality) <- cleaned_dir_names
if (!is.null(colnames(X_clr))) {
  colnames(X_clr) <- sapply(colnames(X_clr), extract_taxon_name)
}

# Generate colors matching PCoA
group_colors <- scales::hue_pal()(length(groups_unique))
names(group_colors) <- groups_unique

# Create VIP plot
cat("\nCreating VIP plot...\n")
vip_plot <- plot_vip(melsi_omnibus, top_n = 15, 
                     title = "SKIOME: Top 15 Taxa by MeLSI Feature Weights",
                     directionality = TRUE)

# Update VIP colors to match PCoA
plot_data <- vip_plot$data
if ("Directionality" %in% names(plot_data)) {
  directionality_values <- unique(plot_data$Directionality)
  color_map <- group_colors[directionality_values]
  names(color_map) <- directionality_values
  vip_plot <- vip_plot + 
    scale_fill_manual(values = color_map, name = "Group", drop = FALSE) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

# Create PCoA plot
cat("Creating PCoA plot...\n")
dist_matrix <- melsi_omnibus$distance_matrix
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
var_explained <- pcoa_result$eig[1:2] / sum(pcoa_result$eig[abs(pcoa_result$eig) > 1e-10]) * 100

plot_data_pcoa <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Group = as.factor(groups)
)

pcoa_plot <- ggplot(plot_data_pcoa, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 1.2) +
  scale_color_manual(values = group_colors, name = "Group") +
  labs(
    title = "SKIOME: PCoA using MeLSI Distance",
    x = sprintf("PCoA1 (%.1f%%)", var_explained[1]),
    y = sprintf("PCoA2 (%.1f%%)", var_explained[2])
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# Combine plots side-by-side with perfect alignment
cat("Combining plots side-by-side...\n")
combined_plot <- grid.arrange(
  vip_plot,
  pcoa_plot,
  ncol = 2,
  widths = c(1, 1)
)

# Save combined figure
output_dir <- "../../manuscript/figures"
output_path <- file.path(output_dir, "skiome_combined.png")
ggsave(output_path, combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
cat("Combined figure saved to:", output_path, "\n")

# Also save as PostScript for publication
output_path_ps <- file.path(output_dir, "skiome_combined.ps")
postscript(output_path_ps, width = 12, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special")
print(combined_plot)
dev.off()
cat("Combined figure (PostScript) saved to:", output_path_ps, "\n")

cat("\nDone!\n")
