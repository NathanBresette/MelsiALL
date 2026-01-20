# MeLSI Summary Figure - Key Results at a Glance
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Set theme
theme_summary <- theme_minimal() +
  theme(
    text = element_text(size = 11, family = "Arial"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Colors
melsi_red <- "#E74C3C"
melsi_blue <- "#2E86AB"
melsi_gray <- "#95A5A6"
melsi_green <- "#27AE60"
melsi_orange <- "#F39C12"

# Create summary figure with multiple panels
create_summary_figure <- function() {
  
  # Panel 1: Key Performance Metrics
  performance_data <- data.frame(
    Metric = c("F-statistic\nImprovement", "Type I Error\nControl", "Computational\nSpeedup", "Feature\nInterpretability"),
    Value = c("28%", "Perfect", "28.7%", "Automatic"),
    Color = c(melsi_red, melsi_green, melsi_blue, melsi_orange)
  )
  
  p1 <- ggplot(performance_data, aes(x = reorder(Metric, Value), y = 1, fill = Color)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Value), color = "white", fontface = "bold", size = 4) +
    scale_fill_identity() +
    labs(title = "Key Performance Metrics", x = "", y = "") +
    theme_summary +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
    coord_flip()
  
  # Panel 2: Real Data Results
  real_data <- data.frame(
    Dataset = c("Atlas1006\n(Human Gut)", "SoilRep\n(Environmental)"),
    MeLSI_F = c(4.85, 1.49),
    Best_Traditional_F = c(4.73, 0.98),
    Improvement = c("2.5%", "52%")
  )
  
  p2 <- ggplot(real_data, aes(x = Dataset)) +
    geom_col(aes(y = MeLSI_F), fill = melsi_red, alpha = 0.8, width = 0.4, position = position_nudge(x = -0.2)) +
    geom_col(aes(y = Best_Traditional_F), fill = melsi_gray, alpha = 0.8, width = 0.4, position = position_nudge(x = 0.2)) +
    geom_text(aes(y = MeLSI_F, label = paste("MeLSI:", MeLSI_F)), 
              position = position_nudge(x = -0.2), hjust = -0.1, size = 3) +
    geom_text(aes(y = Best_Traditional_F, label = paste("Best:", Best_Traditional_F)), 
              position = position_nudge(x = 0.2), hjust = -0.1, size = 3) +
    labs(title = "Real Data Performance", x = "Dataset", y = "F-statistic") +
    theme_summary +
    coord_flip()
  
  # Panel 3: Validation Summary
  validation_data <- data.frame(
    Test = c("Type I Error\n(Null Data)", "Power Analysis\n(Strong Signal)", "Scalability\n(Large Datasets)"),
    Result = c("Perfect Control", "Appropriate Power", "Efficient"),
    Status = c("✓", "✓", "✓")
  )
  
  p3 <- ggplot(validation_data, aes(x = Test, y = 1)) +
    geom_col(fill = melsi_green, alpha = 0.8) +
    geom_text(aes(label = paste(Status, Result)), color = "white", fontface = "bold", size = 3.5) +
    labs(title = "Statistical Validation", x = "", y = "") +
    theme_summary +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
    coord_flip()
  
  # Combine panels
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 3, 
                               top = textGrob("MeLSI: Key Results Summary", 
                                            gp = gpar(fontsize = 16, fontface = "bold")))
  
  return(combined_plot)
}

# Create and save summary figure
cat("Creating summary figure...\n")
summary_fig <- create_summary_figure()

# Save the figure
ggsave("figure_summary_key_results.png", summary_fig, 
       width = 12, height = 6, dpi = 300, bg = "white")

cat("Summary figure saved as: figure_summary_key_results.png\n")
