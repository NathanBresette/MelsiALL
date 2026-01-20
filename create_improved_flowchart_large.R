# Improved MeLSI Algorithm Flowchart with Larger Boxes
library(ggplot2)

# Set theme with no grid lines
theme_clean <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Colors
melsi_blue <- "#2E86AB"
melsi_gray <- "#95A5A6"

# Create improved flowchart with larger boxes
create_improved_flowchart_large <- function() {
  # Define steps with better descriptions and larger boxes
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
  
  # Create the plot with larger boxes
  p <- ggplot(steps_data, aes(x = x, y = y)) +
    # Draw larger main boxes (increased from 0.4 to 0.6 width, 0.3 to 0.4 height)
    geom_rect(aes(xmin = x - 0.6, xmax = x + 0.6, 
                  ymin = y - 0.4, ymax = y + 0.4),
              fill = melsi_blue, alpha = 0.8, color = "white", linewidth = 2) +
    # Add step numbers (larger font)
    geom_text(aes(x = x, y = y + 0.2, label = paste0(step_num, ".")), 
              color = "white", size = 5, fontface = "bold") +
    # Add titles (larger font)
    geom_text(aes(x = x, y = y, label = title), 
              color = "white", size = 4, fontface = "bold") +
    # Add details (larger font)
    geom_text(aes(x = x, y = y - 0.2, label = details), 
              color = "white", size = 3) +
    # Add arrows (adjusted for larger boxes)
    geom_segment(aes(x = x + 0.6, xend = x + 1.3 - 0.6, y = y, yend = y),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5) +
    # Remove last arrow
    geom_segment(aes(x = 8.5 + 0.6, xend = 10 - 0.6, y = y, yend = y),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5) +
    labs(
      title = "MeLSI's 7-Step Process",
      subtitle = "Ensemble metric learning with rigorous statistical validation"
    ) +
    theme_clean +
    xlim(0.3, 10.7) +
    ylim(0.2, 1.8)
  
  return(p)
}

# Create and save the improved flowchart with larger boxes
cat("Creating improved algorithm flowchart with larger boxes...\n")
improved_flowchart_large <- create_improved_flowchart_large()

# Save the figure
ggsave("figure2_algorithm_flow_large_boxes.png", improved_flowchart_large, 
       width = 16, height = 6, dpi = 300, bg = "white")

cat("Improved flowchart with larger boxes saved as: figure2_algorithm_flow_large_boxes.png\n")










