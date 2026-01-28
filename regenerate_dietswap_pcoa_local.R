#!/usr/bin/env Rscript
# Regenerate DietSwap PCoA with larger visible points

suppressPackageStartupMessages({
  library(microbiome)
  library(ggplot2)
  library(MeLSI)
})

set.seed(42)

# Load DietSwap dataset
cat("Loading DietSwap dataset...\n")
data(dietswap)

# Subset to baseline samples
dietswap_subset <- subset_samples(dietswap, timepoint.within.group == 1)
dietswap_subset <- subset_samples(dietswap_subset, group %in% c("DI", "HE"))
dietswap_subset <- prune_taxa(taxa_sums(dietswap_subset) > 0, dietswap_subset)

# Extract counts and metadata
counts <- as(otu_table(dietswap_subset), "matrix")
if (taxa_are_rows(dietswap_subset)) {
  counts <- t(counts)
}
metadata <- data.frame(sample_data(dietswap_subset))
metadata$group <- droplevels(metadata$group)
counts <- counts[rownames(metadata), ]

groups <- metadata$group
valid_idx <- !is.na(groups)
counts <- counts[valid_idx, , drop = FALSE]
groups <- groups[valid_idx]

cat("Samples:", length(groups), "\n")

# CLR transformation
counts_clr <- log(counts + 1)
counts_clr <- counts_clr - rowMeans(counts_clr)
colnames(counts_clr) <- colnames(counts)

# Run MeLSI
cat("Running MeLSI...\n")
melsi_res <- melsi(counts_clr, groups, n_perms = 200, B = 30, show_progress = FALSE, plot_vip = FALSE)

# PCoA Plot - use plot_pcoa() from MeLSI (default settings)
cat("Creating PCoA plot...\n")
pcoa_plot <- plot_pcoa(melsi_res, counts_clr, groups, title = "DietSwap PCoA with MeLSI")

# Save as PNG directly (better quality than converting from PS)
ggsave("manuscript/figures/dietswap_pcoa.png", pcoa_plot, width = 6, height = 4, dpi = 300, bg = "white")
cat("PCoA plot saved with large visible points!\n")
