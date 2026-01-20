# Final MeLSI Algorithm Flowchart - No Grids, Detailed Summaries
library(ggplot2)

# Set theme with no grid lines
theme_final <- theme_minimal() +
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
melsi_blue <- "#457B9D"
melsi_gray <- "#95A5A6"

# Create final flowchart with detailed summaries
create_final_flowchart <- function() {
  # Define steps with detailed summaries
  steps_data <- data.frame(
    step_num = 1:7,
    title = c("Conservative\nPre-filtering", 
              "Bootstrap\nSampling",
              "Feature\nSubsampling", 
              "Gradient\nOptimization",
              "Ensemble\nAveraging",
              "Robust Distance\nCalculation",
              "Permutation\nTesting"),
    details = c("Keep 70% of features\nwith highest variance\nto reduce noise",
                "Create multiple\ntraining sets by\nresampling data",
                "Use 80% of features\nper weak learner\nto prevent overfitting", 
                "Learn optimal weights\nvia gradient descent\nfor each feature subset",
                "Combine 30 weak\nlearners with\nperformance weighting",
                "Calculate Mahalanobis\ndistances robustly\nwith eigenvalue decomposition",
                "Validate significance\nwith 75 permutations\nfor reliable p-values"),
    x = c(1, 2.5, 4, 5.5, 7, 8.5, 10),
    y = rep(1, 7)
  )
  
  # Create the plot with larger boxes and detailed text
  p <- ggplot(steps_data, aes(x = x, y = y)) +
    # Draw larger main boxes with white fill
    geom_rect(aes(xmin = x - 0.7, xmax = x + 0.7, 
                  ymin = y - 0.45, ymax = y + 0.45),
              fill = melsi_blue, alpha = 0.9, color = melsi_blue, linewidth = 2) +
    # Add step numbers and titles together as one string
    geom_text(aes(x = x, y = y + 0.15, label = paste0(step_num, ". ", title)), 
              color = "black", size = 5, fontface = "bold") +
    # Add detailed summaries in black
    geom_text(aes(x = x, y = y - 0.1, label = details), 
              color = "black", size = 3.5) +
    # Add arrows for first 5 steps
    geom_segment(aes(x = x + 0.7, xend = x + 1.3 - 0.7, y = y, yend = y),
                 data = steps_data[1:5, ], 
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5, inherit.aes = FALSE) +
    # Add arrow between 6 and 7
    geom_segment(x = 8.5 + 0.7, xend = 10 - 0.7, y = 1, yend = 1,
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
                 color = melsi_gray, linewidth = 1.5) +
    labs(
      title = "MeLSI's 7-Step Process",
      subtitle = "Ensemble metric learning with rigorous statistical validation"
    ) +
    theme_final +
    xlim(0.3, 10.7) +
    ylim(0.2, 1.8)
  
  return(p)
}

# Create and save the final flowchart
cat("Creating final algorithm flowchart with detailed summaries...\n")
final_flowchart <- create_final_flowchart()

# Save the figure
ggsave("figure2_algorithm_flow_final.png", final_flowchart, 
       width = 16, height = 6, dpi = 300, bg = "white")

cat("Final flowchart saved as: figure2_algorithm_flow_final.png\n")

