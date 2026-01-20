# Research Impact Text Graphic - Blue Theme, No Boxes
library(ggplot2)

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

# Colors - blue theme only
melsi_blue <- "#2E86AB"
text_gray <- "#2C3E50"
light_blue <- "#EBF3FD"

# Create impact text graphic
create_impact_text_graphic <- function() {
  
  # Create a clean text-based layout
  p <- ggplot() +
    # Add title
    annotate("text", x = 0.5, y = 0.9, 
             label = "Research Impact", 
             size = 5, fontface = "bold", color = melsi_blue, hjust = 0) +
    
    # Add subtitle
    annotate("text", x = 0.5, y = 0.85, 
             label = "How MeLSI Advances Microbiome Research", 
             size = 4, color = text_gray, hjust = 0) +
    
    # Add impact points with bullet points
    annotate("text", x = 0.5, y = 0.75, 
             label = "• Enables more powerful microbiome studies", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.7, 
             label = "• Provides interpretable feature importance for biological insights", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Applicable to any microbiome dataset (human, environmental, clinical)", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Open-source R package available for widespread use", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.55, 
             label = "• Transforms statistical results into actionable biological knowledge", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    # Add a subtle blue accent line
    annotate("segment", x = 0.45, xend = 0.45, y = 0.5, yend = 0.8,
             color = melsi_blue, linewidth = 3) +
    
    # Set limits
    xlim(0, 1) +
    ylim(0.4, 1) +
    theme_impact
  
  return(p)
}

# Create and save the impact text graphic
cat("Creating research impact text graphic...\n")
impact_text <- create_impact_text_graphic()

# Save the figure
ggsave("figure6_research_impact_text.png", impact_text, 
       width = 10, height = 6, dpi = 300, bg = "white")

cat("Research impact text graphic saved as: figure6_research_impact_text.png\n")










