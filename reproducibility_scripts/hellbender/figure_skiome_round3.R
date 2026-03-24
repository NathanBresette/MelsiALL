#!/usr/bin/env Rscript
# ==============================================================================
# Round 3 Figure 3: SKIOME (Atopic Dermatitis / Healthy / Psoriasis)
# Standalone script — loads all results from saved RData/CSV, no MeLSI rerun
# ==============================================================================

options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(vegan)
  library(scales)
  library(gridExtra)
  library(grid)
  library(MeLSI)
})

set.seed(42)

output_dir <- "figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Helper: PCoA plot
# ==============================================================================
make_pcoa_plot <- function(dist_mat, groups, group_colors, title, subtitle = NULL) {
  pcoa_res <- ape::pcoa(dist_mat)
  explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
  pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$Group <- groups

  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.4) +
    stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1.0) +
    scale_color_manual(values = group_colors, name = "Group") +
    labs(
      title = title, subtitle = subtitle,
      x = paste0("PCoA1 (", explained[1], "%)"),
      y = paste0("PCoA2 (", explained[2], "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
}

# ==============================================================================
# Helper: VIP plot
# ==============================================================================
clean_taxon_name <- function(x) {
  # Split on semicolon, take last non-empty element
  parts <- strsplit(as.character(x), ";")[[1]]
  parts <- trimws(parts)
  parts <- parts[nchar(parts) > 0]
  name <- if (length(parts) > 0) parts[length(parts)] else x
  # Strip rank prefixes like g__, f__, s__, k__, p__, c__, o__, sk__
  name <- gsub("^[a-z]+__", "", name)
  # Truncate to 35 chars
  if (nchar(name) > 35) name <- paste0(substr(name, 1, 32), "...")
  name
}

make_vip_plot <- function(melsi_res, group_colors, title, n_top = 15) {
  raw_names <- names(melsi_res$feature_weights)
  clean_names <- sapply(raw_names, clean_taxon_name, USE.NAMES = FALSE)

  vip_df <- data.frame(
    Feature = clean_names,
    RawFeature = raw_names,
    Weight = as.numeric(melsi_res$feature_weights)
  ) %>%
    arrange(desc(Weight)) %>%
    slice_head(n = n_top) %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)))

  has_dir <- !is.null(melsi_res$directionality) && length(melsi_res$directionality) > 0

  if (has_dir) {
    vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$RawFeature)]
    vip_df$Directionality <- ifelse(is.na(vip_df$Directionality), "Unknown", vip_df$Directionality)
    unique_dirs <- unique(vip_df$Directionality[vip_df$Directionality != "Unknown"])

    color_map <- character()
    for (dir in unique_dirs) {
      group_name <- gsub("Higher in ", "", dir)
      if (group_name %in% names(group_colors)) {
        color_map[dir] <- group_colors[group_name]
      } else {
        color_map[dir] <- "#95A5A6"
      }
    }
    if (any(vip_df$Directionality == "Unknown")) color_map["Unknown"] <- "#95A5A6"

    p <- ggplot(vip_df, aes(x = Weight, y = Feature, fill = Directionality)) +
      geom_col() +
      scale_fill_manual(values = color_map, name = "Direction", drop = FALSE) +
      geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 2.5) +
      labs(title = title, x = "Learned Metric Weight", y = NULL)
  } else {
    p <- ggplot(vip_df, aes(x = Weight, y = Feature)) +
      geom_col(fill = "#1f77b4") +
      geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 2.5) +
      labs(title = title, x = "Learned Metric Weight", y = NULL)
  }

  p + theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
}

# ==============================================================================
# Load SKIOME data and results
# ==============================================================================
cat("Loading SKIOME data and MeLSI results...\n")
load("skiome_data_loaded.RData")       # loads: counts, metadata
load("skiome_omnibus_results_backup.RData")  # loads: melsi_omnibus (or similar)

counts_skiome <- counts
group_skiome  <- metadata$Group
cat("Samples:", nrow(counts_skiome), " Groups:", paste(unique(group_skiome), collapse = ", "), "\n")

counts_skiome_clr <- log(counts_skiome + 1)
counts_skiome_clr <- counts_skiome_clr - rowMeans(counts_skiome_clr)

groups_skiome <- unique(group_skiome)
colors_skiome <- scales::hue_pal()(length(groups_skiome))
names(colors_skiome) <- groups_skiome

# Resolve MeLSI result object name
if (exists("melsi_omnibus")) {
  melsi_res_skiome <- melsi_omnibus
} else if (exists("melsi_result")) {
  melsi_res_skiome <- melsi_result
} else {
  cat("Loaded objects:", paste(ls(), collapse = ", "), "\n")
  stop("Could not identify MeLSI result object.")
}

# Load distance matrix from CSV
# The CSV stores the lower triangle as a flat vector (130305 values = 511*510/2)
cat("Loading SKIOME MeLSI distance matrix from CSV...\n")
n_skiome <- nrow(counts_skiome)  # 511
d_vec <- read.csv("skiome_omnibus_distance_matrix.csv")[[1]]
cat("Vector length:", length(d_vec), "  Expected:", n_skiome * (n_skiome - 1) / 2, "\n")
cat("Range:", range(d_vec, na.rm = TRUE), "  Non-finite:", sum(!is.finite(d_vec)), "\n")

dist_melsi_skiome <- structure(
  d_vec,
  Size   = n_skiome,
  Labels = rownames(counts_skiome),
  Diag   = FALSE,
  Upper  = FALSE,
  class  = "dist"
)
rm(d_vec); gc()

# Compute traditional distances (for PCoA ordination only)
cat("Computing Euclidean distance...\n")
dist_euc <- dist(counts_skiome_clr)

cat("Computing Bray-Curtis distance...\n")
dist_bray <- vegdist(counts_skiome, method = "bray")

# ==============================================================================
# Build panels — use verified F/p from Hellbender CSVs (skip adonis2 to save memory)
# ==============================================================================
gc()
cat("Memory before panels:", format(sum(gc()[,2]) * 8 / 1024, digits=3), "MB used\n")
cat("Building panels...\n")

# Validate distance matrices before PCoA
cat("dist_melsi: range =", range(as.numeric(dist_melsi_skiome), na.rm=TRUE), " NaN =", sum(is.nan(as.numeric(dist_melsi_skiome))), "\n")
cat("dist_euc:   range =", range(as.numeric(dist_euc), na.rm=TRUE), " NaN =", sum(is.nan(as.numeric(dist_euc))), "\n")
cat("dist_bray:  range =", range(as.numeric(dist_bray), na.rm=TRUE), " NaN =", sum(is.nan(as.numeric(dist_bray))), "\n")

cat("Making VIP plot...\n")
p_vip <- make_vip_plot(melsi_res_skiome, colors_skiome,
                       "(A) SKIOME: Learned Metric Weights")
cat("VIP done\n")

cat("Making MeLSI PCoA...\n")
p_melsi <- make_pcoa_plot(dist_melsi_skiome, group_skiome, colors_skiome,
                          "(B) MeLSI PCoA", "F=4.972, p=0.005")
cat("MeLSI PCoA done\n"); gc()

cat("Making Euclidean PCoA...\n")
p_euc <- make_pcoa_plot(dist_euc, group_skiome, colors_skiome,
                        "(C) Euclidean (CLR) PCoA", "F=4.897, p=0.001")
cat("Euclidean PCoA done\n"); gc()

cat("Making Bray-Curtis PCoA...\n")
p_bray <- make_pcoa_plot(dist_bray, group_skiome, colors_skiome,
                         "(D) Bray-Curtis PCoA", "F=16.275, p=0.001")
cat("Bray-Curtis PCoA done\n"); gc()

# Free distance matrices before rendering
rm(dist_euc, dist_bray, dist_melsi_skiome, counts_skiome, counts_skiome_clr)
gc()

cat("Combining panels...\n")
combined <- grid.arrange(p_vip, p_melsi, p_euc, p_bray, ncol = 2, widths = c(1.2, 1))

# Save PDF
pdf_path <- file.path(output_dir, "skiome_figure3.pdf")
ggsave(pdf_path, combined, width = 14, height = 10, device = "pdf")
cat("Saved:", pdf_path, "\n")

# Try PNG
png_path <- file.path(output_dir, "skiome_figure3.png")
tryCatch({
  grDevices::png(png_path, width = 14, height = 10, units = "in", res = 300, type = "cairo")
  grid.draw(combined)
  dev.off()
  cat("Saved:", png_path, "\n")
}, error = function(e) {
  cat("PNG save failed:", conditionMessage(e), "\n")
})

cat("\nSKIOME figure complete.\n")
