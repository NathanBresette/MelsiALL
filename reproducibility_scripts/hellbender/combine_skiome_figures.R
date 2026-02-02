#!/usr/bin/env Rscript
# Combine SKIOME VIP and PCoA plots into a single figure

library(ggplot2)
library(gridExtra)

# Set output directory
output_dir <- "../../manuscript/figures"

# Read the VIP plot (PNG)
vip_img <- png::readPNG("skiome_omnibus_vip_plot_fixed.png")
vip_grob <- grid::rasterGrob(vip_img, interpolate = TRUE)

# Read the PCoA plot (PNG - need to convert from PS or use existing)
# Check if we have a PNG version of PCoA
pcoa_png_path <- file.path(output_dir, "skiome_pcoa_plot.png")
if (file.exists(pcoa_png_path)) {
  pcoa_img <- png::readPNG(pcoa_png_path)
  pcoa_grob <- grid::rasterGrob(pcoa_img, interpolate = TRUE)
} else {
  # If no PNG, we'll need to load the PCoA from the backup and regenerate
  cat("PCoA PNG not found, will need to regenerate from backup\n")
  # For now, create a placeholder
  pcoa_grob <- NULL
}

# Create combined plot
if (!is.null(pcoa_grob)) {
  combined_plot <- grid.arrange(
    vip_grob,
    pcoa_grob,
    ncol = 2,
    widths = c(1, 1)
  )
  
  # Save combined figure
  output_path <- file.path(output_dir, "skiome_combined.png")
  ggsave(output_path, combined_plot, width = 12, height = 6, dpi = 300, bg = "white")
  cat("Combined figure saved to:", output_path, "\n")
} else {
  cat("Could not create combined figure - PCoA PNG missing\n")
}
