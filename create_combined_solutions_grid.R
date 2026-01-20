# Combined Solutions Grid - 2x2 Layout
library(ggplot2)
library(gridExtra)
library(grid)

# Set theme to match your style
theme_combined <- theme_minimal() +
  theme(
    text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

# Colors - blue theme only
melsi_blue <- "#2E86AB"
text_gray <- "#2C3E50"

# Create individual panels
create_tech_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.9, 
             label = "Technical Innovations", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• First adaptive distance metric learning", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.75, 
             label = "• Ensemble approach prevents overfitting", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.7, 
             label = "• Conservative pre-filtering", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Rigorous permutation testing", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Automatic feature importance", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.55, yend = 0.85,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.5, 1) + theme_combined
  return(p)
}

create_impact_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.9, 
             label = "Research Impact", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• More powerful microbiome studies", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.75, 
             label = "• Interpretable biological insights", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.7, 
             label = "• Universal applicability", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Open-source accessibility", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Actionable biological knowledge", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.55, yend = 0.85,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.5, 1) + theme_combined
  return(p)
}

create_performance_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.9, 
             label = "Performance Advantages", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• 28% improvement in F-statistics", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.75, 
             label = "• Perfect Type I error control", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.7, 
             label = "• 28.7% computational speedup", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Validated on real datasets", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Superior to traditional methods", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.55, yend = 0.85,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.5, 1) + theme_combined
  return(p)
}

create_validation_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.9, 
             label = "Rigorous Validation", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• Synthetic data testing", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.75, 
             label = "• Real microbiome datasets", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.7, 
             label = "• Power analysis validation", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Multiple testing correction", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Computational efficiency", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.55, yend = 0.85,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.5, 1) + theme_combined
  return(p)
}

# Create the combined grid
create_combined_grid <- function() {
  tech_panel <- create_tech_panel()
  impact_panel <- create_impact_panel()
  performance_panel <- create_performance_panel()
  validation_panel <- create_validation_panel()
  
  # Arrange in 2x2 grid
  combined_plot <- grid.arrange(
    tech_panel, impact_panel,
    performance_panel, validation_panel,
    ncol = 2, nrow = 2,
    top = textGrob("MeLSI Solutions: Comprehensive Advantages", 
                   gp = gpar(fontsize = 16, fontface = "bold", col = melsi_blue))
  )
  
  return(combined_plot)
}

# Create and save the combined grid
cat("Creating combined solutions grid...\n")
combined_grid <- create_combined_grid()

# Save the figure
ggsave("figure6_combined_solutions_grid.png", combined_grid, 
       width = 12, height = 8, dpi = 300, bg = "white")

cat("Combined solutions grid saved as: figure6_combined_solutions_grid.png\n")










