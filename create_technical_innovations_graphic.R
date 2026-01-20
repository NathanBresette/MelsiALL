# Technical Innovations Graphic - Blue Theme, No Boxes
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

# Create technical innovations graphic
create_technical_innovations_graphic <- function() {
  
  # Create a clean text-based layout
  p <- ggplot() +
    # Add title
    annotate("text", x = 0.5, y = 0.9, 
             label = "Key Innovations", 
             size = 5, fontface = "bold", color = melsi_blue, hjust = 0) +
    
    # Add subtitle
    annotate("text", x = 0.5, y = 0.85, 
             label = "Technical Advances in MeLSI", 
             size = 4, color = text_gray, hjust = 0) +
    
    # Add innovation points with bullet points
    annotate("text", x = 0.5, y = 0.75, 
             label = "• First adaptive distance metric learning for microbiome data", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.7, 
             label = "• Ensemble approach prevents overfitting and improves generalization", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.65, 
             label = "• Conservative pre-filtering balances statistical power with efficiency", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.6, 
             label = "• Rigorous permutation testing ensures statistical reliability", 
             size = 3.5, color = text_gray, hjust = 0) +
    
    annotate("text", x = 0.5, y = 0.55, 
             label = "• Automatic feature importance visualization provides biological insights", 
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

# Create and save the technical innovations graphic
cat("Creating technical innovations graphic...\n")
tech_innovations <- create_technical_innovations_graphic()

# Save the figure
ggsave("figure6_technical_innovations.png", tech_innovations, 
       width = 10, height = 6, dpi = 300, bg = "white")

cat("Technical innovations graphic saved as: figure6_technical_innovations.png\n")










