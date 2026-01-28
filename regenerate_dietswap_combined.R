#!/usr/bin/env Rscript
# Regenerate DietSwap VIP and PCoA figures combined side-by-side with matching heights

suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(MeLSI)
  library(gridExtra)
  library(grid)
})

set.seed(42)

# Load DietSwap dataset
cat("Loading DietSwap dataset...\n")
data(dietswap)

# Subset to baseline samples
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

groups <- metadata$group
valid_idx <- !is.na(groups)
counts <- counts[valid_idx, , drop = FALSE]
groups <- groups[valid_idx]

cat("Samples:", length(groups), "\n")

# CLR transformation
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)
colnames(counts_clr) <- colnames(counts)

# Run MeLSI
cat("Running MeLSI...\n")
melsi_res <- melsi(counts_clr, groups, n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)

# Create PCoA plot first to extract colors
cat("Creating PCoA plot...\n")
pcoa_plot <- plot_pcoa(melsi_res, counts_clr, groups, title = "DietSwap PCoA with MeLSI")

# Extract colors from PCoA plot
pcoa_build <- ggplot_build(pcoa_plot)
group_levels <- levels(as.factor(groups))
# Get colors in order of group levels
pcoa_colors <- scales::hue_pal()(length(group_levels))
# Match colors to groups based on the plot
pcoa_data <- pcoa_build$data[[1]]
unique_colors <- unique(pcoa_data$colour[order(pcoa_data$group)])
if (length(unique_colors) == length(group_levels)) {
  pcoa_colors <- unique_colors[order(unique(pcoa_data$group))]
  names(pcoa_colors) <- group_levels
} else {
  # Fallback: use default ggplot2 colors
  pcoa_colors <- scales::hue_pal()(length(group_levels))
  names(pcoa_colors) <- group_levels
}

# Create VIP plot using MeLSI function
cat("Creating VIP plot...\n")
vip_plot <- plot_vip(melsi_res, top_n = 15, title = "DietSwap VIP with MeLSI")

# Update VIP plot colors to match PCoA
# Get directionality from melsi results
if (!is.null(melsi_res$directionality)) {
  vip_directionality <- melsi_res$directionality
  # Map directionality groups to PCoA colors
  vip_color_map <- character()
  for (group_name in group_levels) {
    # Find directionality entries that match this group
    matching_dirs <- grep(group_name, names(vip_directionality), value = TRUE, ignore.case = TRUE)
    if (length(matching_dirs) > 0) {
      for (dir in matching_dirs) {
        vip_color_map[dir] <- pcoa_colors[group_name]
      }
    }
  }
  # Also check for exact matches
  for (dir in unique(vip_directionality)) {
    if (dir %in% group_levels && !dir %in% names(vip_color_map)) {
      vip_color_map[dir] <- pcoa_colors[dir]
    }
  }
  
  # Update VIP plot with matching colors
  if (length(vip_color_map) > 0) {
    vip_plot <- vip_plot + scale_fill_manual(values = vip_color_map, name = "Higher in")
  }
}

# Ensure both plots have the same height
# Extract the gtable from both plots
vip_grob <- ggplotGrob(vip_plot)
pcoa_grob <- ggplotGrob(pcoa_plot)

# Get heights
vip_heights <- vip_grob$heights
pcoa_heights <- pcoa_grob$heights

# Use the maximum height for both
max_heights <- unit.pmax(vip_heights, pcoa_heights)
vip_grob$heights <- max_heights
pcoa_grob$heights <- max_heights

# Combine side-by-side
cat("Combining plots...\n")
combined_plot <- arrangeGrob(
  vip_grob,
  pcoa_grob,
  ncol = 2,
  widths = c(1, 1)
)

# Save combined plot
output_file <- "manuscript/figures/dietswap_combined.png"
ggsave(output_file, combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
cat("Combined plot saved to:", output_file, "\n")

# Also save individual plots for reference
ggsave("manuscript/figures/dietswap_vip.png", vip_plot, width = 6, height = 6, dpi = 300, bg = "white")
ggsave("manuscript/figures/dietswap_pcoa.png", pcoa_plot, width = 6, height = 6, dpi = 300, bg = "white")
cat("Individual plots saved\n")
