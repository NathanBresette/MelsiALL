# Updated Layout Diagram with Beta Diversity Introduction
library(ggplot2)
library(dplyr)

# Set theme
theme_layout <- theme_void() +
  theme(
    text = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Colors
primary_blue <- "#2E86AB"
accent_red <- "#E74C3C"
text_gray <- "#2C3E50"
light_gray <- "#95A5A6"
success_green <- "#27AE60"

# Create updated content layout
create_updated_layout <- function() {
  
  # Define content areas with new order
  content_data <- data.frame(
    area = c("Title", "Intro", "Fig1", "Fig4", "Fig2", "Fig6", "Fig3", "Summary", "Fig5", "Conclusions", "Contact"),
    x = c(24, 8.4, 8.4, 8.4, 24, 24, 39.6, 39.6, 39.6, 12, 36),
    y = c(33.3, 25.2, 15.8, 6.4, 25.2, 6.4, 25.2, 15.8, 6.4, 2.7, 2.7),
    width = c(48, 14, 14, 14, 12, 12, 14, 14, 14, 20, 20),
    height = c(5.4, 8, 6, 8, 15, 8, 8, 6, 8, 5.4, 5.4),
    color = c(primary_blue, primary_blue, accent_red, accent_red, primary_blue, primary_blue, 
              accent_red, accent_red, success_green, primary_blue, accent_red)
  )
  
  # Create content plot
  p <- ggplot(content_data, aes(x = x, y = y)) +
    # Draw content areas
    geom_rect(aes(xmin = x - width/2, xmax = x + width/2, 
                  ymin = y - height/2, ymax = y + height/2,
                  fill = color),
              alpha = 0.4, color = "black", linewidth = 0.5) +
    # Add labels
    geom_text(aes(label = area), color = "black", fontface = "bold", size = 3) +
    # Add detailed descriptions
    geom_text(aes(x = x, y = y + height/2 - 0.3, 
                  label = case_when(
                    area == "Intro" ~ "Beta Diversity & Distance Metrics\nIntroduction",
                    area == "Fig1" ~ "Figure 1: Problem & Solution",
                    area == "Fig2" ~ "Figure 2: Algorithm Flow",
                    area == "Fig3" ~ "Figure 3: Performance",
                    area == "Fig4" ~ "Figure 4: Feature Importance",
                    area == "Fig5" ~ "Figure 5: Validation",
                    area == "Fig6" ~ "Figure 6: Innovation Box",
                    TRUE ~ area
                  )),
              color = "black", size = 2.5) +
    # Set colors
    scale_fill_identity() +
    # Set limits
    xlim(0, 48) +
    ylim(0, 36) +
    labs(title = "Updated Layout: Beta Diversity Introduction → Figure 1") +
    theme_layout
  
  return(p)
}

# Create the updated layout
cat("Creating updated layout diagram...\n")
updated_layout <- create_updated_layout()

# Save the diagram
ggsave("poster_updated_layout.png", updated_layout, width = 12, height = 9, dpi = 300, bg = "white")

cat("Updated layout diagram saved as: poster_updated_layout.png\n")










