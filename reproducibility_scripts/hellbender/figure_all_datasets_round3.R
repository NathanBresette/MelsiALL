#!/usr/bin/env Rscript
# ==============================================================================
# Round 3 Figure Generation: 4-Panel Figures for All Datasets
# ==============================================================================
# Generates Figures 1-3 with layout:
#   (A) VIP bar plot with directionality
#   (B) PCoA using MeLSI learned distance
#   (C) PCoA using Euclidean distance (CLR)
#   (D) PCoA using Bray-Curtis dissimilarity
#
# SKIOME: loads existing MeLSI results from RData (avoids expensive rerun)
# Atlas1006/DietSwap: runs MeLSI then saves results as RData for future use
# ==============================================================================

suppressPackageStartupMessages({
  library(microbiome)
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

# Output directory
output_dir <- "figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Helper: Create PCoA plot for any distance matrix
# ==============================================================================
make_pcoa_plot <- function(dist_mat, groups, group_colors, title, subtitle = NULL) {
  pcoa_res <- ape::pcoa(dist_mat)
  explained <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
  pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$Group <- groups

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.4) +
    stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1.0) +
    scale_color_manual(values = group_colors, name = "Group") +
    labs(
      title = title,
      subtitle = subtitle,
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
  return(p)
}

# ==============================================================================
# Helper: Create VIP plot with directionality
# ==============================================================================
make_vip_plot <- function(melsi_res, group_colors, title, n_top = 15) {
  vip_df <- data.frame(
    Feature = names(melsi_res$feature_weights),
    Weight = as.numeric(melsi_res$feature_weights)
  ) %>%
    arrange(desc(Weight)) %>%
    slice_head(n = n_top) %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)))

  has_dir <- !is.null(melsi_res$directionality) && length(melsi_res$directionality) > 0

  if (has_dir) {
    vip_df$Directionality <- melsi_res$directionality[as.character(vip_df$Feature)]
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
      labs(title = title, x = "Learned Metric Weight", y = NULL) +
      theme_bw() +
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
  } else {
    p <- ggplot(vip_df, aes(x = Weight, y = Feature)) +
      geom_col(fill = "#1f77b4") +
      geom_text(aes(label = sprintf("%.3f", Weight)), hjust = -0.05, size = 2.5) +
      labs(title = title, x = "Learned Metric Weight", y = NULL) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 8)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
  }
  return(p)
}

# ==============================================================================
# Helper: Build 4-panel figure from MeLSI result + data
# ==============================================================================
build_figure <- function(melsi_res, dist_melsi, counts_raw, counts_clr, groups,
                         group_colors, dataset_name, output_filename) {

  cat("\n==============================================================================\n")
  cat("Building figure:", dataset_name, "\n")
  cat("==============================================================================\n")
  cat("  MeLSI: F =", round(melsi_res$F_observed, 3), " p =", round(melsi_res$p_value, 3), "\n")

  # Euclidean distance on CLR data (deterministic)
  dist_euc <- dist(counts_clr)
  perm_euc <- adonis2(dist_euc ~ groups, permutations = 999)
  cat("  Euclidean: F =", round(perm_euc$F[1], 3), " p =", round(perm_euc$`Pr(>F)`[1], 3), "\n")

  # Bray-Curtis distance on raw counts (deterministic)
  dist_bray <- vegdist(counts_raw, method = "bray")
  perm_bray <- adonis2(dist_bray ~ groups, permutations = 999)
  cat("  Bray-Curtis: F =", round(perm_bray$F[1], 3), " p =", round(perm_bray$`Pr(>F)`[1], 3), "\n")

  # Panel A: VIP
  p_vip <- make_vip_plot(melsi_res, group_colors,
                         paste0("(A) ", dataset_name, ": Learned Metric Weights"))

  # Panel B: MeLSI PCoA
  p_melsi <- make_pcoa_plot(dist_melsi, groups, group_colors,
                            "(B) MeLSI PCoA",
                            paste0("F=", round(melsi_res$F_observed, 3),
                                   ", p=", round(melsi_res$p_value, 3)))

  # Panel C: Euclidean PCoA
  p_euc <- make_pcoa_plot(dist_euc, groups, group_colors,
                          "(C) Euclidean (CLR) PCoA",
                          paste0("F=", round(perm_euc$F[1], 3),
                                 ", p=", round(perm_euc$`Pr(>F)`[1], 3)))

  # Panel D: Bray-Curtis PCoA
  p_bray <- make_pcoa_plot(dist_bray, groups, group_colors,
                           "(D) Bray-Curtis PCoA",
                           paste0("F=", round(perm_bray$F[1], 3),
                                  ", p=", round(perm_bray$`Pr(>F)`[1], 3)))

  # Combine into 2x2 layout
  combined <- grid.arrange(
    p_vip, p_melsi,
    p_euc, p_bray,
    ncol = 2,
    widths = c(1.2, 1)
  )

  # Save PNG
  png_path <- file.path(output_dir, paste0(output_filename, ".png"))
  ggsave(png_path, combined, width = 14, height = 10, dpi = 300)
  cat("Saved:", png_path, "\n")

  # Save TIFF for journal submission
  tiff_path <- file.path(output_dir, paste0(output_filename, ".tif"))
  ggsave(tiff_path, combined, width = 14, height = 10, dpi = 300, compression = "lzw")
  cat("Saved:", tiff_path, "\n")

  return(list(melsi_res = melsi_res, perm_euc = perm_euc, perm_bray = perm_bray))
}

# ==============================================================================
# FIGURE 1: Atlas1006 (Male vs Female)
# ==============================================================================
cat("\n*** FIGURE 1: Atlas1006 ***\n")

data(atlas1006)
counts_atlas <- t(abundances(atlas1006))
meta_atlas <- meta(atlas1006)
sex <- meta_atlas$sex
valid_idx <- !is.na(sex)
counts_atlas <- counts_atlas[valid_idx, , drop = FALSE]
sex <- sex[valid_idx]

counts_atlas_clr <- log(counts_atlas + 1)
counts_atlas_clr <- counts_atlas_clr - rowMeans(counts_atlas_clr)

groups_atlas <- unique(sex)
colors_atlas <- scales::hue_pal()(length(groups_atlas))
names(colors_atlas) <- groups_atlas

# Check for saved results
if (file.exists("atlas1006_melsi_results.RData")) {
  cat("Loading saved Atlas1006 MeLSI results...\n")
  load("atlas1006_melsi_results.RData")
} else {
  cat("Running MeLSI for Atlas1006 (no saved results found)...\n")
  set.seed(42)
  melsi_res_atlas <- melsi(counts_atlas_clr, sex,
                           n_perms = 200, B = 30, m_frac = 0.8,
                           show_progress = TRUE, plot_vip = FALSE)
  # Compute MeLSI distance
  counts_atlas_filtered <- MeLSI:::apply_conservative_prefiltering(counts_atlas_clr, sex, filter_frac = 0.7)
  dist_melsi_atlas <- MeLSI:::calculate_mahalanobis_dist_robust(counts_atlas_filtered, melsi_res_atlas$metric_matrix)
  # Save for future use
  save(melsi_res_atlas, dist_melsi_atlas, file = "atlas1006_melsi_results.RData")
  cat("Saved Atlas1006 MeLSI results to atlas1006_melsi_results.RData\n")
}

res1 <- build_figure(melsi_res_atlas, dist_melsi_atlas,
                     counts_atlas, counts_atlas_clr, sex, colors_atlas,
                     "Atlas1006", "atlas1006_figure1")

# ==============================================================================
# FIGURE 2: DietSwap (Western vs High-fiber)
# ==============================================================================
cat("\n*** FIGURE 2: DietSwap ***\n")

data(dietswap)
meta_diet <- meta(dietswap)
baseline_idx <- meta_diet$timepoint.within.group == 1
counts_diet <- t(abundances(dietswap))[baseline_idx, , drop = FALSE]
group_diet <- meta_diet$nationality[baseline_idx]

counts_diet_clr <- log(counts_diet + 1)
counts_diet_clr <- counts_diet_clr - rowMeans(counts_diet_clr)

groups_diet <- unique(group_diet)
colors_diet <- scales::hue_pal()(length(groups_diet))
names(colors_diet) <- groups_diet

# Check for saved results
if (file.exists("dietswap_melsi_results.RData")) {
  cat("Loading saved DietSwap MeLSI results...\n")
  load("dietswap_melsi_results.RData")
} else {
  cat("Running MeLSI for DietSwap (no saved results found)...\n")
  set.seed(42)
  melsi_res_diet <- melsi(counts_diet_clr, group_diet,
                          n_perms = 200, B = 30, m_frac = 0.8,
                          show_progress = TRUE, plot_vip = FALSE)
  # Compute MeLSI distance
  counts_diet_filtered <- MeLSI:::apply_conservative_prefiltering(counts_diet_clr, group_diet, filter_frac = 0.7)
  dist_melsi_diet <- MeLSI:::calculate_mahalanobis_dist_robust(counts_diet_filtered, melsi_res_diet$metric_matrix)
  # Save for future use
  save(melsi_res_diet, dist_melsi_diet, file = "dietswap_melsi_results.RData")
  cat("Saved DietSwap MeLSI results to dietswap_melsi_results.RData\n")
}

res2 <- build_figure(melsi_res_diet, dist_melsi_diet,
                     counts_diet, counts_diet_clr, group_diet, colors_diet,
                     "DietSwap", "dietswap_figure2")

# ==============================================================================
# FIGURE 3: SKIOME (Atopic Dermatitis / Healthy / Psoriasis) — ALL 3 GROUPS
# ==============================================================================
cat("\n*** FIGURE 3: SKIOME (omnibus, 3 groups) ***\n")

# Load pre-parsed SKIOME data
skiome_rdata <- "skiome_data_loaded.RData"
skiome_results_rdata <- "skiome_omnibus_results_backup.RData"
skiome_dist_csv <- "skiome_omnibus_distance_matrix.csv"

if (file.exists(skiome_rdata) && file.exists(skiome_results_rdata)) {
  load(skiome_rdata)  # loads counts_skiome, group_skiome
  load(skiome_results_rdata)  # loads melsi_omnibus (or similar)
  cat("Loaded existing SKIOME data and MeLSI results\n")

  counts_skiome_clr <- log(counts_skiome + 1)
  counts_skiome_clr <- counts_skiome_clr - rowMeans(counts_skiome_clr)

  groups_skiome <- unique(group_skiome)
  colors_skiome <- scales::hue_pal()(length(groups_skiome))
  names(colors_skiome) <- groups_skiome

  # Load MeLSI distance matrix from CSV if available
  if (file.exists(skiome_dist_csv)) {
    cat("Loading SKIOME MeLSI distance matrix from CSV...\n")
    dist_mat_raw <- as.matrix(read.csv(skiome_dist_csv, row.names = 1))
    dist_melsi_skiome <- as.dist(dist_mat_raw)
  } else {
    cat("Computing SKIOME MeLSI distance from saved metric matrix...\n")
    counts_skiome_filtered <- MeLSI:::apply_conservative_prefiltering(counts_skiome_clr, group_skiome, filter_frac = 0.7)
    dist_melsi_skiome <- MeLSI:::calculate_mahalanobis_dist_robust(counts_skiome_filtered, melsi_omnibus$metric_matrix)
  }

  # Figure out which variable name the RData used for the MeLSI result
  # Common names: melsi_omnibus, melsi_result, res
  if (exists("melsi_omnibus")) {
    melsi_res_skiome <- melsi_omnibus
  } else if (exists("melsi_result")) {
    melsi_res_skiome <- melsi_result
  } else {
    # List all objects loaded and pick the one that looks like a MeLSI result
    loaded_objs <- ls()
    cat("Loaded objects:", paste(loaded_objs, collapse = ", "), "\n")
    stop("Could not identify MeLSI result object from RData. Check variable names.")
  }

  res3 <- build_figure(melsi_res_skiome, dist_melsi_skiome,
                       counts_skiome, counts_skiome_clr, group_skiome, colors_skiome,
                       "SKIOME", "skiome_figure3")
} else {
  cat("SKIOME data or results not found. Skipping Figure 3.\n")
  cat("Expected files: skiome_data_loaded.RData, skiome_omnibus_results_backup.RData\n")
}

cat("\n==============================================================================\n")
cat("All figures generated successfully.\n")
cat("==============================================================================\n")
