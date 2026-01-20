# Demo script with realistic taxonomic names for microbiome data
# Load MeLSI package
library(MeLSI)

# Generate synthetic microbiome data with realistic taxonomic names
set.seed(42)
n_samples <- 60
n_taxa <- 50  # Reduced for better demonstration

# Create realistic bacterial species names commonly found in microbiome studies
realistic_taxa <- c(
    # Bacteroides genus (common gut bacteria)
    "Bacteroides_vulgatus", "Bacteroides_thetaiotaomicron", "Bacteroides_fragilis",
    "Bacteroides_ovatus", "Bacteroides_uniformis", "Bacteroides_caccae",
    
    # Firmicutes phylum
    "Faecalibacterium_prausnitzii", "Ruminococcus_bromii", "Ruminococcus_gnavus",
    "Roseburia_intestinalis", "Coprococcus_comes", "Dorea_longicatena",
    "Blautia_wexlerae", "Clostridium_innocuum", "Clostridium_leptum",
    
    # Bifidobacterium genus
    "Bifidobacterium_longum", "Bifidobacterium_adolescentis", "Bifidobacterium_bifidum",
    "Bifidobacterium_breve", "Bifidobacterium_catenulatum",
    
    # Lactobacillus genus
    "Lactobacillus_acidophilus", "Lactobacillus_casei", "Lactobacillus_plantarum",
    "Lactobacillus_rhamnosus", "Lactobacillus_reuteri",
    
    # Other important genera
    "Akkermansia_muciniphila", "Prevotella_copri", "Prevotella_oralis",
    "Parabacteroides_distasonis", "Alistipes_putredinis", "Odoribacter_splanchnicus",
    "Barnesiella_intestinihominis", "Porphyromonas_gingivalis",
    
    # Proteobacteria
    "Escherichia_coli", "Klebsiella_pneumoniae", "Enterobacter_cloacae",
    
    # Actinobacteria
    "Collinsella_aerofaciens", "Eggerthella_lenta", "Slackia_equolifaciens",
    
    # Additional diverse species
    "Fusobacterium_nucleatum", "Veillonella_parvula", "Streptococcus_salivarius",
    "Enterococcus_faecalis", "Staphylococcus_epidermidis", "Propionibacterium_acnes",
    "Corynebacterium_striatum", "Micrococcus_luteus", "Sarcina_ventriculi",
    "Methanobrevibacter_smithii", "Desulfovibrio_piger", "Succinivibrio_dextrinosolvens"
)

# Use the realistic names (repeat if we have more taxa than names)
if (n_taxa <= length(realistic_taxa)) {
    taxa_names <- realistic_taxa[1:n_taxa]
} else {
    # If we need more names, create variations
    taxa_names <- c(realistic_taxa, paste0("Unknown_bacterium_", 1:(n_taxa - length(realistic_taxa))))
}

# Create base microbiome data (log-normal distribution)
X <- matrix(rlnorm(n_samples * n_taxa, meanlog = 2, sdlog = 1), 
            nrow = n_samples, ncol = n_taxa)

# IMPORTANT: Set column names to realistic taxonomic names
colnames(X) <- taxa_names

# Create group labels
y <- c(rep("Control", n_samples/2), rep("Treatment", n_samples/2))

# Add signal to first 10 taxa in Treatment group (simulate differential abundance)
# These will be the most important features
treatment_indices <- which(y == "Treatment")
X[treatment_indices, 1:10] <- X[treatment_indices, 1:10] * 2.0  # Strong signal

cat("=== MeLSI Demo with Realistic Taxonomic Names ===\n")
cat("Sample size:", n_samples, "samples\n")
cat("Number of taxa:", n_taxa, "taxa\n")
cat("Groups: Control (n=", n_samples/2, "), Treatment (n=", n_samples/2, ")\n\n")

# Check that column names are preserved
cat("First 10 taxonomic names:\n")
print(head(colnames(X), 10))
cat("\n")

# CLR transformation (recommended for microbiome data)
X_clr <- log(X + 1)
X_clr <- X_clr - rowMeans(X_clr)

# IMPORTANT: Preserve column names after transformation
colnames(X_clr) <- colnames(X)

# Run MeLSI analysis (automatically detects 2 groups)
cat("Running MeLSI analysis...\n")
results <- melsi(X_clr, y, n_perms = 50, B = 20, m_frac = 0.8, show_progress = TRUE, plot_vip = TRUE)

# Display results
cat("\n=== MeLSI Results ===\n")
cat(sprintf("F-statistic: %.4f\n", results$F_observed))
cat(sprintf("P-value: %.4f\n", results$p_value))
cat(sprintf("Significant: %s\n", ifelse(results$p_value < 0.05, "Yes", "No")))

# Show feature names in results
cat("\nTop 10 most important taxa:\n")
top_10 <- head(sort(results$feature_weights, decreasing = TRUE), 10)
for (i in 1:length(top_10)) {
    cat(sprintf("  %2d. %-30s (weight: %.4f)\n", i, names(top_10)[i], top_10[i]))
}

cat("\n=== Analysis Complete ===\n")
cat("Notice how the actual bacterial species names are now displayed\n")
cat("instead of generic 'Feature_X' names!\n")











