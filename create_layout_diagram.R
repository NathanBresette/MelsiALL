# Create Visual Layout Diagram for MeLSI Poster
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)

# Set theme
theme_layout <- theme_void() +
  theme(
    text = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Colors from our palette
bg_color <- "#FFFFFF"
primary_blue <- "#2E86AB"
accent_red <- "#E74C3C"
text_gray <- "#2C3E50"
light_gray <- "#95A5A6"
success_green <- "#27AE60"

# Create layout diagram
create_layout_diagram <- function() {
  
  # Define poster dimensions (48" x 36")
  poster_width <- 48
  poster_height <- 36
  
  # Define section dimensions
  header_height <- 5.4  # 15% of 36"
  main_height <- 25.2   # 70% of 36"
  bottom_height <- 5.4  # 15% of 36"
  
  left_width <- 16.8    # 35% of 48"
  center_width <- 14.4  # 30% of 48"
  right_width <- 16.8   # 35% of 48"
  
  # Create layout data
  layout_data <- data.frame(
    section = c("Header", "Left Col", "Center Col", "Right Col", "Bottom"),
    x = c(24, 8.4, 24, 39.6, 24),
    y = c(33.3, 18.9, 18.9, 18.9, 2.7),
    width = c(48, 16.8, 14.4, 16.8, 48),
    height = c(5.4, 25.2, 25.2, 25.2, 5.4),
    color = c(primary_blue, light_gray, accent_red, success_green, text_gray)
  )
  
  # Create the layout plot
  p <- ggplot(layout_data, aes(x = x, y = y)) +
    # Draw main sections
    geom_rect(aes(xmin = x - width/2, xmax = x + width/2, 
                  ymin = y - height/2, ymax = y + height/2,
                  fill = color),
              alpha = 0.3, color = "black", linewidth = 1) +
    # Add section labels
    geom_text(aes(label = section), color = "black", fontface = "bold", size = 4) +
    # Add content labels
    geom_text(aes(label = case_when(
      section == "Header" ~ "Title, Name, Contact",
      section == "Left Col" ~ "Fig 1: Problem\nIntroduction\nFig 4: VIP",
      section == "Center Col" ~ "Fig 2: Algorithm\nFig 6: Innovation",
      section == "Right Col" ~ "Fig 3: Results\nSummary\nFig 5: Validation",
      section == "Bottom" ~ "Conclusions\nFuture Work\nContact"
    )), color = "black", size = 3, vjust = 0.5) +
    # Add dimensions
    geom_text(aes(x = x, y = y + height/2 + 0.5, 
                  label = paste0(width, '" x ', height, '"')),
              color = "black", size = 2.5) +
    # Set colors
    scale_fill_identity() +
    # Set limits
    xlim(0, 48) +
    ylim(0, 36) +
    labs(title = "MeLSI Poster Layout (48\" x 36\")") +
    theme_layout
  
  return(p)
}

# Create detailed content layout
create_content_layout <- function() {
  
  # Define content areas
  content_data <- data.frame(
    area = c("Title", "Fig1", "Intro", "Fig4", "Fig2", "Fig6", "Fig3", "Summary", "Fig5", "Conclusions", "Contact"),
    x = c(24, 8.4, 8.4, 8.4, 24, 24, 39.6, 39.6, 39.6, 12, 36),
    y = c(33.3, 25.2, 15.8, 6.4, 25.2, 6.4, 25.2, 15.8, 6.4, 2.7, 2.7),
    width = c(48, 14, 14, 14, 12, 12, 14, 14, 14, 20, 20),
    height = c(5.4, 8, 6, 8, 15, 8, 8, 6, 8, 5.4, 5.4),
    color = c(primary_blue, accent_red, light_gray, accent_red, primary_blue, primary_blue, 
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
    # Add figure numbers
    geom_text(aes(x = x, y = y + height/2 - 0.3, 
                  label = case_when(
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
    labs(title = "Detailed Content Layout with Figure Placement") +
    theme_layout
  
  return(p)
}

# Create color palette reference
create_color_palette <- function() {
  
  color_data <- data.frame(
    color_name = c("Background", "Primary Blue", "Accent Red", "Text Gray", "Light Gray", "Success Green"),
    hex_code = c("#FFFFFF", "#2E86AB", "#E74C3C", "#2C3E50", "#95A5A6", "#27AE60"),
    usage = c("Main background", "Headers, algorithm boxes", "MeLSI results, highlights", "Main text", "Traditional methods", "Validation checkmarks"),
    x = 1:6,
    y = 1
  )
  
  p <- ggplot(color_data, aes(x = x, y = y)) +
    geom_rect(aes(xmin = x - 0.4, xmax = x + 0.4, ymin = y - 0.3, ymax = y + 0.3,
                  fill = hex_code),
              color = "black", linewidth = 1) +
    geom_text(aes(label = color_name, y = y - 0.5), fontface = "bold", size = 3) +
    geom_text(aes(label = hex_code, y = y + 0.5), size = 2.5) +
    geom_text(aes(label = usage, y = y + 0.8), size = 2) +
    scale_fill_identity() +
    xlim(0.5, 6.5) +
    ylim(0, 2) +
    labs(title = "MeLSI Color Palette") +
    theme_layout
  
  return(p)
}

# Create all layout diagrams
cat("Creating layout diagrams...\n")

layout_diagram <- create_layout_diagram()
content_layout <- create_content_layout()
color_palette <- create_color_palette()

# Save diagrams
ggsave("poster_layout_overview.png", layout_diagram, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("poster_content_layout.png", content_layout, width = 12, height = 9, dpi = 300, bg = "white")
ggsave("poster_color_palette.png", color_palette, width = 10, height = 4, dpi = 300, bg = "white")

cat("Layout diagrams saved:\n")
cat("- poster_layout_overview.png\n")
cat("- poster_content_layout.png\n")
cat("- poster_color_palette.png\n")
