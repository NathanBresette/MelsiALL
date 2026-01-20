# Research Impact Graphic for MeLSI
library(ggplot2)
library(gridExtra)

# Set theme to match your style
theme_impact <- theme_minimal() +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# Colors matching your palette
melsi_blue <- "#2E86AB"
melsi_red <- "#E74C3C"
melsi_gray <- "#95A5A6"
melsi_green <- "#27AE60"
text_gray <- "#2C3E50"

# Create research impact graphic
create_research_impact_graphic <- function() {
  
  # Define impact areas with icons and descriptions
  impact_data <- data.frame(
    area = c("More Powerful\nStudies", "Biological\nInsights", "Universal\nApplicability", "Open Source\nAccess", "Actionable\nKnowledge"),
    description = c("Enables more sensitive\ndetection of patterns", "Interpretable feature\nimportance visualization", "Any microbiome dataset\n(human, environmental, clinical)", "R package available\nfor widespread use", "Transforms statistics\ninto biological knowledge"),
    x = c(1, 2, 3, 4, 5),
    y = rep(1, 5),
    color = c(melsi_blue, melsi_red, melsi_green, melsi_gray, melsi_blue)
  )
  
  # Create the main plot
  p <- ggplot(impact_data, aes(x = x, y = y)) +
    # Draw impact boxes
    geom_rect(aes(xmin = x - 0.4, xmax = x + 0.4, 
                  ymin = y - 0.3, ymax = y + 0.3,
                  fill = color),
              alpha = 0.8, color = "white", linewidth = 2) +
    # Add area titles
    geom_text(aes(x = x, y = y + 0.1, label = area), 
              color = "white", size = 3.5, fontface = "bold") +
    # Add descriptions
    geom_text(aes(x = x, y = y - 0.1, label = description), 
              color = "white", size = 2.8) +
    # Add impact icons/symbols
    geom_text(aes(x = x, y = y + 0.25, label = c("⚡", "🔬", "🌍", "📦", "💡")), 
              color = "white", size = 6) +
    # Set colors
    scale_fill_identity() +
    # Set limits
    xlim(0.5, 5.5) +
    ylim(0.4, 1.6) +
    labs(
      title = "Research Impact",
      subtitle = "How MeLSI Advances Microbiome Research"
    ) +
    theme_impact
  
  return(p)
}

# Create and save the research impact graphic
cat("Creating research impact graphic...\n")
impact_graphic <- create_research_impact_graphic()

# Save the figure
ggsave("figure6_research_impact.png", impact_graphic, 
       width = 12, height = 6, dpi = 300, bg = "white")

cat("Research impact graphic saved as: figure6_research_impact.png\n")










