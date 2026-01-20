# MeLSI Poster Figures
# Figures 1, 2, and 3 for the MeLSI poster presentation

library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Set consistent theme and colors
theme_poster <- theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", size = 0.5)
  )

# Color scheme
melsi_red <- "#E74C3C"
melsi_blue <- "#2E86AB"
melsi_gray <- "#95A5A6"
melsi_green <- "#27AE60"
melsi_orange <- "#F39C12"

# =============================================================================
# FIGURE 1: The Problem & Solution - Fixed vs Adaptive Distance Metrics
# =============================================================================

create_figure1 <- function() {
  # Create data for the comparison
  taxa <- c("Bacteroides", "Faecalibacterium", "Ruminococcus", "Bifidobacterium", 
            "Lactobacillus", "Akkermansia", "Prevotella", "Roseburia")
  
  # Traditional fixed weights (all equal)
  traditional_weights <- rep(1, length(taxa))
  
  # MeLSI learned weights (example from your results)
  melsi_weights <- c(2.1, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6)
  
  # Create comparison data
  comparison_data <- data.frame(
    Taxa = rep(taxa, 2),
    Method = rep(c("Traditional\n(Fixed Weights)", "MeLSI\n(Learned Weights)"), each = length(taxa)),
    Weight = c(traditional_weights, melsi_weights),
    Order = rep(1:length(taxa), 2)
  )
  
  # Create the plot
  p1 <- ggplot(comparison_data, aes(x = reorder(Taxa, Order), y = Weight, fill = Method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Traditional\n(Fixed Weights)" = melsi_gray, 
                                "MeLSI\n(Learned Weights)" = melsi_red)) +
    labs(
      title = "Fixed vs. Adaptive Distance Metrics",
      subtitle = "Traditional methods treat all taxa equally, MeLSI learns optimal weights",
      x = "Microbial Taxa",
      y = "Distance Metric Weight",
      fill = "Method"
    ) +
    theme_poster +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "bottom"
    ) +
    coord_flip()
  
  return(p1)
}

# =============================================================================
# FIGURE 2: MeLSI Algorithm Flow - 7-Step Process
# =============================================================================

create_figure2 <- function() {
  # Create flowchart data
  steps <- c("1. Conservative\nPre-filtering\n(70% features)", 
             "2. Bootstrap\nSampling", 
             "3. Feature\nSubsampling\n(80% per learner)",
             "4. Gradient\nOptimization", 
             "5. Ensemble\nAveraging\n(30 learners)",
             "6. Robust Distance\nCalculation", 
             "7. Permutation\nTesting")
  
  # Create positions for flowchart
  flowchart_data <- data.frame(
    Step = steps,
    x = c(1, 2, 3, 4, 5, 6, 7),
    y = rep(1, 7),
    width = 0.8,
    height = 0.6
  )
  
  # Create the flowchart
  p2 <- ggplot(flowchart_data, aes(x = x, y = y)) +
    # Draw boxes
    geom_rect(aes(xmin = x - width/2, xmax = x + width/2, 
                  ymin = y - height/2, ymax = y + height/2),
              fill = melsi_blue, alpha = 0.7, color = "white", size = 1) +
    # Add step numbers and text
    geom_text(aes(label = Step), color = "white", size = 3, fontface = "bold") +
    # Add arrows between steps
    geom_segment(aes(x = x + width/2, xend = x + 1 - width/2, y = y, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = melsi_gray, size = 1) +
    # Remove last arrow
    geom_segment(aes(x = 6 + width/2, xend = 7 - width/2, y = y, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = melsi_gray, size = 1) +
    labs(
      title = "MeLSI's 7-Step Process",
      subtitle = "Ensemble metric learning with rigorous validation"
    ) +
    theme_poster +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    xlim(0.5, 7.5) +
    ylim(0.5, 1.5)
  
  return(p2)
}

# =============================================================================
# FIGURE 3: Performance Comparison - MeLSI vs Traditional Methods
# =============================================================================

create_figure3 <- function() {
  # Performance data from your results
  methods <- c("MeLSI", "Euclidean", "Bray-Curtis", "Jaccard", "Weighted UniFrac")
  
  # Atlas1006 results
  atlas_F <- c(4.85, 4.73, 4.44, 3.41, 1.67)
  atlas_p <- c(0.013, 0.001, 0.001, 0.001, 0.18)
  
  # Synthetic Strong results  
  synthetic_F <- c(1.67, 1.60, 1.28, 1.15, 1.07)
  synthetic_p <- c(0.013, 0.001, 0.012, 0.012, 0.37)
  
  # Create comparison data
  comparison_data <- data.frame(
    Method = rep(methods, 2),
    Dataset = rep(c("Atlas1006 (Real Data)", "Synthetic Strong"), each = length(methods)),
    F_statistic = c(atlas_F, synthetic_F),
    P_value = c(atlas_p, synthetic_p),
    Significant = c(atlas_p < 0.05, synthetic_p < 0.05)
  )
  
  # Create the plot
  p3 <- ggplot(comparison_data, aes(x = reorder(Method, F_statistic), y = F_statistic, 
                                   fill = ifelse(Method == "MeLSI", "MeLSI", "Traditional"))) +
    geom_col(alpha = 0.8) +
    facet_wrap(~Dataset, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("MeLSI" = melsi_red, "Traditional" = melsi_gray)) +
    labs(
      title = "MeLSI Outperforms Traditional Methods",
      subtitle = "Higher F-statistics indicate better group separation",
      x = "Method",
      y = "F-statistic",
      fill = "Method Type"
    ) +
    theme_poster +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    ) +
    coord_flip()
  
  # Add significance indicators
  p3 <- p3 + 
    geom_text(aes(label = ifelse(Significant, "*", "")), 
              hjust = -0.1, vjust = 0.5, size = 4, color = melsi_green)
  
  return(p3)
}

# =============================================================================
# CREATE AND SAVE ALL FIGURES
# =============================================================================

# Create all figures
cat("Creating Figure 1: Problem & Solution...\n")
fig1 <- create_figure1()

cat("Creating Figure 2: Algorithm Flow...\n")
fig2 <- create_figure2()

cat("Creating Figure 3: Performance Comparison...\n")
fig3 <- create_figure3()

# Save figures
ggsave("figure1_problem_solution.png", fig1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("figure2_algorithm_flow.png", fig2, width = 12, height = 4, dpi = 300, bg = "white")
ggsave("figure3_performance_comparison.png", fig3, width = 10, height = 8, dpi = 300, bg = "white")

cat("\nAll figures saved successfully!\n")
cat("Files created:\n")
cat("- figure1_problem_solution.png\n")
cat("- figure2_algorithm_flow.png\n") 
cat("- figure3_performance_comparison.png\n")

# Display figures (optional - comment out if running in headless mode)
# print(fig1)
# print(fig2) 
# print(fig3)










