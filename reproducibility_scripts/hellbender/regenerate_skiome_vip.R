#!/usr/bin/env Rscript
# Regenerate SKIOME VIP plot with proper 3-group coloring using backup files

# Load required libraries
library(ggplot2)

# Source the MeLSI functions
source("../../github/R/melsi_robust.R")

# Script should be run from the hellbender directory
# where the backup files are located

cat("Loading backup files...\n")

# Load the backup RData file
load("skiome_omnibus_results_backup.RData")

cat("Loaded backup files:\n")
cat("  - melsi_omnibus results object\n")
cat("  - X_clr: CLR-transformed data\n")
cat("  - groups: Group labels\n\n")

# Calculate directionality for 3 groups properly
cat("Calculating directionality for 3 groups...\n")
groups_unique <- unique(groups)
cat("  Groups found:", paste(groups_unique, collapse = ", "), "\n")

# Calculate mean abundance for each group
mean_by_group <- list()
for (g in groups_unique) {
  group_idx <- which(groups == g)
  mean_by_group[[as.character(g)]] <- colMeans(X_clr[group_idx, , drop = FALSE])
}

# Create matrix of means
mean_matrix <- do.call(rbind, mean_by_group)
rownames(mean_matrix) <- groups_unique

# For each feature, find which group has highest mean
directionality_calc <- apply(mean_matrix, 2, function(x) {
  max_idx <- which.max(x)
  return(groups_unique[max_idx])
})
names(directionality_calc) <- colnames(X_clr)

cat("  Directionality calculated for", length(directionality_calc), "features\n")
cat("  Group distribution:\n")
for (g in groups_unique) {
  count <- sum(directionality_calc == g)
  cat("    ", g, ":", count, "features\n")
}

# Update the results object with calculated directionality
melsi_omnibus$directionality <- directionality_calc

# Clean up taxon names - extract most specific level (genus > family > order)
cat("\nCleaning taxon names...\n")
if (!is.null(names(melsi_omnibus$feature_weights))) {
  extract_taxon_name <- function(full_name) {
    # Remove "sk__Bacteria;" prefix
    name <- gsub("^sk__Bacteria;", "", full_name)
    
    # Try to extract genus (g__)
    if (grepl("g__[^;]+", name)) {
      genus <- gsub(".*g__([^;]+).*", "\\1", name)
      if (genus != "" && genus != name) {
        return(genus)
      }
    }
    
    # If no genus, try family (f__)
    if (grepl("f__[^;]+", name)) {
      family <- gsub(".*f__([^;]+).*", "\\1", name)
      if (family != "" && family != name) {
        return(paste0(family, " (family)"))
      }
    }
    
    # If no family, try order (o__)
    if (grepl("o__[^;]+", name)) {
      order <- gsub(".*o__([^;]+).*", "\\1", name)
      if (order != "" && order != name) {
        return(paste0(order, " (order)"))
      }
    }
    
    # Fallback: return cleaned original
    return(gsub("^k__;", "", name))
  }
  
  # Apply to all feature names
  cleaned_names <- sapply(names(melsi_omnibus$feature_weights), extract_taxon_name)
  names(melsi_omnibus$feature_weights) <- cleaned_names
  
  # Also clean directionality names
  cleaned_dir_names <- sapply(names(melsi_omnibus$directionality), extract_taxon_name)
  names(melsi_omnibus$directionality) <- cleaned_dir_names
  
  # Also clean column names in X_clr if they exist
  if (!is.null(colnames(X_clr))) {
    colnames(X_clr) <- sapply(colnames(X_clr), extract_taxon_name)
  }
  
  cat("  Extracted most specific taxonomic level (genus > family > order)\n")
  cat("  Sample cleaned names:\n")
  top15_idx <- head(order(melsi_omnibus$feature_weights, decreasing = TRUE), 5)
  print(cleaned_names[top15_idx])
}

# Generate VIP plot with proper 3-group coloring matching PCoA colors
cat("\nGenerating VIP plot with 3-group directionality...\n")

# Use same color scheme as PCoA plot for consistency (matching Atlas1006/DietSwap approach)
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(scales)

# Generate colors using same palette as PCoA (hue_pal)
group_colors <- scales::hue_pal()(length(groups_unique))
names(group_colors) <- groups_unique
cat("  Group colors (matching PCoA):\n")
for (g in groups_unique) {
  cat("    ", g, ":", group_colors[g], "\n")
}

# Use the plot_vip function from MeLSI, then update colors to match PCoA
vip_plot <- plot_vip(melsi_omnibus, top_n = 15, 
                     title = "SKIOME: Top 15 Taxa by MeLSI Feature Weights",
                     directionality = TRUE)

# Update the plot colors to match PCoA
# Extract the plot data and update fill colors
plot_data <- vip_plot$data
if (!is.null(plot_data$fill) || !is.null(plot_data$Directionality)) {
  # Get current fill mapping
  if ("Directionality" %in% names(plot_data)) {
    # Map directionality to group colors
    directionality_values <- unique(plot_data$Directionality)
    color_map <- group_colors[directionality_values]
    names(color_map) <- directionality_values
    
    # Update the plot's scale_fill_manual
    vip_plot <- vip_plot + 
      scale_fill_manual(values = color_map, name = "Group", drop = FALSE)
  }
}

# Save as PostScript
postscript("skiome_omnibus_vip_plot_fixed.ps", width = 10, height = 8, 
           horizontal = FALSE, onefile = FALSE, paper = "special")
print(vip_plot)
dev.off()

cat("VIP plot saved to: skiome_omnibus_vip_plot_fixed.ps\n")

# Also save as PNG for easier viewing
png("skiome_omnibus_vip_plot_fixed.png", width = 10, height = 8, units = "in", res = 300)
print(vip_plot)
dev.off()

cat("VIP plot saved to: skiome_omnibus_vip_plot_fixed.png\n")
cat("\nDone!\n")
