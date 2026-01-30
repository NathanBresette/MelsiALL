#!/usr/bin/env Rscript
# Combine SKIOME VIP and PCoA plots into a single figure

library(png)
library(grid)
library(gridExtra)

# Read the PNG files
vip_img <- readPNG("skiome_vip_plot.png")
pcoa_img <- readPNG("skiome_pcoa_plot.png")

# Convert to grobs
vip_grob <- rasterGrob(vip_img)
pcoa_grob <- rasterGrob(pcoa_img)

# Combine side by side
combined <- grid.arrange(vip_grob, pcoa_grob, ncol = 2, widths = c(1, 1))

# Save as PNG
png("skiome_combined.png", width = 16, height = 8, units = "in", res = 300)
grid.arrange(vip_grob, pcoa_grob, ncol = 2, widths = c(1, 1))
dev.off()

cat("Combined figure saved to: skiome_combined.png\n")
