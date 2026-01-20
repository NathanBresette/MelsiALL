# Improved MeLSI Algorithm Flowchart for Poster
library(ggplot2)
library(grid)
library(gridExtra)

# Set theme
theme_flowchart <- theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Colors
melsi_blue <- "#2E86AB"
melsi_red <- "#E74C3C"
melsi_gray <- "#95A5A6"

# Create improved flowchart
create_improved_flowchart <- function() {
  # Define steps with better descriptions
  steps_data <- data.frame(
    step_num = 1:7,
    title = c("Conservative\nPre-filtering", 
              "Bootstrap\nSampling",
              "Feature\nSubsampling", 
              "Gradient\nOptimization",
              "Ensemble\nAveraging",
              "Robust Distance\nCalculation",
              "Permutation\nTesting"),
    details = c("Keep 70% of features\nwith highest variance",
                "Create multiple\ntraining sets",
                "Use 80% of features\nper weak learner", 
                "Learn optimal weights\nvia gradient descent",
                "Combine 30 weak\nlearners with weighting",
                "Calculate Mahalanobis\ndistances robustly",
                "Validate with 75\npermutations"),
    x = c(1, 2.5, 4, 5.5, 7, 8.5, 10),
    y = rep(1, 7)
  )
  
  # Create the plot
  p <- ggplot(steps_data, aes(x = x, y = y)) +
    # Draw main boxes
    geom_rect(aes(xmin = x - 0.4, xmax = x + 0.4, 
                  ymin = y - 0.3, ymax = y + 0.3),
              fill = melsi_blue, alpha = 0.8, color = "white", linewidth = 2) +
    # Add step numbers
    geom_text(aes(x = x, y = y + 0.15, label = paste0(step_num, ".")), 
              color = "white", size = 4, fontface = "bold") +
    # Add titles
    geom_text(aes(x = x, y = y, label = title), 
              color = "white", size = 3.5, fontface = "bold") +
    # Add details
    geom_text(aes(x = x, y = y - 0.15, label = details), 
              color = "white", size = 2.5) +
    # Add arrows
    geom_segment(aes(x = x + 0.4, xend = x + 1.1 - 0.4, y = y, yend = y),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5) +
    # Remove last arrow
    geom_segment(aes(x = 8.5 + 0.4, xend = 10 - 0.4, y = y, yend = y),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5) +
    labs(
      title = "MeLSI's 7-Step Process",
      subtitle = "Ensemble metric learning with rigorous statistical validation"
    ) +
    theme_flowchart +
    xlim(0.5, 10.5) +
    ylim(0.3, 1.7)
  
  return(p)
}

# Create and save the improved flowchart
cat("Creating improved algorithm flowchart...\n")
improved_flowchart <- create_improved_flowchart()

# Save the figure
ggsave("figure2_algorithm_flow_improved.png", improved_flowchart, 
       width = 14, height = 5, dpi = 300, bg = "white")

cat("Improved flowchart saved as: figure2_algorithm_flow_improved.png\n")










