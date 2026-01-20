# Simple MeLSI Summary Figure
library(ggplot2)
library(gridExtra)

# Set theme
theme_simple <- theme_minimal() +
  theme(
    text = element_text(size = 11),
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

# Create simple summary figure
create_simple_summary <- function() {
  
  # Panel 1: Key Results
  results_data <- data.frame(
    Result = c("28% F-statistic\nImprovement", "Perfect Type I\nError Control", "28.7% Computational\nSpeedup", "Automatic Feature\nImportance"),
    Value = c("Atlas1006", "Null Data", "Pre-filtering", "VIP Plots"),
    fill_color = c(melsi_red, melsi_green, melsi_blue, melsi_gray)
  )
  
  p1 <- ggplot(results_data, aes(x = reorder(Result, Value), y = 1, fill = fill_color)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Value), color = "white", fontface = "bold", size = 3.5) +
    scale_fill_identity() +
    labs(title = "Key Performance Results", x = "", y = "") +
    theme_simple +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
    coord_flip()
  
  # Panel 2: Real Data Comparison
  comparison_data <- data.frame(
    Dataset = c("Atlas1006\n(Human Gut)", "SoilRep\n(Environmental)"),
    MeLSI = c(4.85, 1.49),
    Best_Traditional = c(4.73, 0.98)
  )
  
  p2 <- ggplot(comparison_data, aes(x = Dataset)) +
    geom_col(aes(y = MeLSI), fill = melsi_red, alpha = 0.8, width = 0.4, position = position_nudge(x = -0.2)) +
    geom_col(aes(y = Best_Traditional), fill = melsi_gray, alpha = 0.8, width = 0.4, position = position_nudge(x = 0.2)) +
    geom_text(aes(y = MeLSI, label = paste("MeLSI:", MeLSI)), 
              position = position_nudge(x = -0.2), hjust = -0.1, size = 3) +
    geom_text(aes(y = Best_Traditional, label = paste("Best:", Best_Traditional)), 
              position = position_nudge(x = 0.2), hjust = -0.1, size = 3) +
    labs(title = "Real Data Performance", x = "Dataset", y = "F-statistic") +
    theme_simple +
    coord_flip()
  
  # Combine panels
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  
  return(combined_plot)
}

# Create and save
cat("Creating simple summary figure...\n")
simple_summary <- create_simple_summary()

ggsave("figure_summary_simple.png", simple_summary, 
       width = 10, height = 5, dpi = 300, bg = "white")

cat("Simple summary figure saved as: figure_summary_simple.png\n")










