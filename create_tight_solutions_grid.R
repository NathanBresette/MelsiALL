# Tight Solutions Grid - Reduced Spacing
library(ggplot2)
library(gridExtra)
library(grid)

# Set theme with tighter spacing
theme_tight <- theme_minimal() +
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
    plot.margin = margin(4, 4, 4, 4)  # Reduced margins
  )

# Colors - blue theme only
melsi_blue <- "#2E86AB"
text_gray <- "#2C3E50"

# Create individual panels with tighter spacing
create_tech_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, 
             label = "Technical Innovations", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.88, 
             label = "• First adaptive distance metric learning", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.84, 
             label = "• Ensemble approach prevents overfitting", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• Conservative pre-filtering", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.76, 
             label = "• Rigorous permutation testing", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.72, 
             label = "• Automatic feature importance", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.68, yend = 0.92,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.65, 1) + theme_tight
  return(p)
}

create_impact_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, 
             label = "Research Impact", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.88, 
             label = "• More powerful microbiome studies", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.84, 
             label = "• Interpretable biological insights", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• Universal applicability", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.76, 
             label = "• Open-source accessibility", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.72, 
             label = "• Actionable biological knowledge", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.68, yend = 0.92,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.65, 1) + theme_tight
  return(p)
}

create_performance_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, 
             label = "Performance Advantages", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.88, 
             label = "• 28% improvement in F-statistics", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.84, 
             label = "• Perfect Type I error control", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• 28.7% computational speedup", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.76, 
             label = "• Validated on real datasets", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.72, 
             label = "• Superior to traditional methods", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.68, yend = 0.92,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.65, 1) + theme_tight
  return(p)
}

create_validation_panel <- function() {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.95, 
             label = "Rigorous Validation", 
             size = 4, fontface = "bold", color = melsi_blue, hjust = 0) +
    annotate("text", x = 0.5, y = 0.88, 
             label = "• Synthetic data testing", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.84, 
             label = "• Real microbiome datasets", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.8, 
             label = "• Power analysis validation", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.76, 
             label = "• Multiple testing correction", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("text", x = 0.5, y = 0.72, 
             label = "• Computational efficiency", 
             size = 3, color = text_gray, hjust = 0) +
    annotate("segment", x = 0.45, xend = 0.45, y = 0.68, yend = 0.92,
             color = melsi_blue, linewidth = 2) +
    xlim(0, 1) + ylim(0.65, 1) + theme_tight
  return(p)
}

# Create the tight combined grid
create_tight_grid <- function() {
  tech_panel <- create_tech_panel()
  impact_panel <- create_impact_panel()
  performance_panel <- create_performance_panel()
  validation_panel <- create_validation_panel()
  
  # Arrange in 2x2 grid with tighter spacing
  combined_plot <- grid.arrange(
    tech_panel, impact_panel,
    performance_panel, validation_panel,
    ncol = 2, nrow = 2,
    top = textGrob("MeLSI Solutions: Comprehensive Advantages", 
                   gp = gpar(fontsize = 16, fontface = "bold", col = melsi_blue)),
    padding = unit(0.1, "line")  # Reduced padding between panels
  )
  
  return(combined_plot)
}

# Create and save the tight grid
cat("Creating tight solutions grid...\n")
tight_grid <- create_tight_grid()

# Save the figure
ggsave("figure6_tight_solutions_grid.png", tight_grid, 
       width = 12, height = 8, dpi = 300, bg = "white")

cat("Tight solutions grid saved as: figure6_tight_solutions_grid.png\n")










