# Supplementary Tables for MeLSI Manuscript

## Supplementary Table S1: Recovery of True Signal Taxa

This table shows MeLSI's ability to recover true signal taxa using learned feature weights across varying effect sizes and sample sizes.

| Effect Size | Sample Size | n | Precision at 5 | Precision at 10 | Precision at 20 | Recall at 5 | Recall at 10 | Recall at 20 | Mean Rank | AUC-ROC |
|-------------|-------------|---|----------------|-----------------|-----------------|-------------|--------------|--------------|-----------|---------|
| Small | 50 | 50 | 0.104 | 0.084 | 0.060 | 0.104 | 0.168 | 0.240 | 50.3 | 0.641 |
| Small | 100 | 50 | 0.104 | 0.084 | 0.064 | 0.104 | 0.168 | 0.256 | 49.3 | 0.649 |
| Small | 200 | 50 | 0.148 | 0.122 | 0.086 | 0.148 | 0.244 | 0.344 | 46.2 | 0.673 |
| Medium | 50 | 50 | 0.356 | 0.262 | 0.197 | 0.178 | 0.262 | 0.394 | 38.9 | 0.733 |
| Medium | 100 | 50 | 0.520 | 0.364 | 0.259 | 0.260 | 0.364 | 0.518 | 31.2 | 0.794 |
| Medium | 200 | 50 | 0.660 | 0.462 | 0.303 | 0.330 | 0.462 | 0.606 | 25.1 | 0.842 |
| Large | 50 | 50 | 0.876 | 0.728 | 0.560 | 0.219 | 0.364 | 0.560 | 26.2 | 0.858 |
| Large | 100 | 50 | 0.976 | 0.912 | 0.705 | 0.244 | 0.456 | 0.705 | 19.7 | 0.914 |
| Large | 200 | 50 | 1.000 | 0.988 | 0.804 | 0.250 | 0.494 | 0.804 | 14.4 | 0.960 |

**Abbreviations:** n, number of simulations; Precision at k, proportion of top-k features that are true signals; Recall at k, proportion of true signals found in top-k features; Mean Rank, average rank of true signal features (lower is better); AUC-ROC, area under receiver operating characteristic curve for classifying signal vs. non-signal taxa based on learned feature weights.

---

## Supplementary Table S2: Individual Method Comparisons for Power Analysis

This table shows detailed comparisons between MeLSI and each of the five traditional methods individually across all effect sizes and sample sizes. These comparisons support the rank calculations shown in Table 2.

| Effect Size | Sample Size | Traditional Method | MeLSI Power (%) | MeLSI Mean F | Traditional Power (%) | Traditional Mean F | Power Difference (%) | F Difference |
|-------------|-------------|-------------------|-----------------|--------------|----------------------|-------------------|---------------------|--------------|
| Small | 50 | Euclidean | 6 | 1.230 | 8 | 1.014 | -2 | 0.216 |
| Small | 50 | Bray-Curtis | 6 | 1.230 | 20 | 1.059 | -14 | 0.171 |
| Small | 50 | Jaccard | 6 | 1.230 | 0 | 0.987 | 6 | 0.243 |
| Small | 50 | Weighted UniFrac | 6 | 1.230 | 10 | 1.069 | -4 | 0.161 |
| Small | 50 | Unweighted UniFrac | 6 | 1.230 | 6 | 1.005 | 0 | 0.226 |
| Small | 100 | Euclidean | 10 | 1.342 | 8 | 1.041 | 2 | 0.301 |
| Small | 100 | Bray-Curtis | 10 | 1.342 | 20 | 1.095 | -10 | 0.247 |
| Small | 100 | Jaccard | 10 | 1.342 | 6 | 1.012 | 4 | 0.330 |
| Small | 100 | Weighted UniFrac | 10 | 1.342 | 8 | 1.048 | 2 | 0.294 |
| Small | 100 | Unweighted UniFrac | 10 | 1.342 | 4 | 1.012 | 6 | 0.331 |
| Small | 200 | Euclidean | 16 | 1.432 | 16 | 1.074 | 0 | 0.358 |
| Small | 200 | Bray-Curtis | 16 | 1.432 | 54 | 1.182 | -38 | 0.250 |
| Small | 200 | Jaccard | 16 | 1.432 | 6 | 0.988 | 10 | 0.444 |
| Small | 200 | Weighted UniFrac | 16 | 1.432 | 20 | 1.200 | -4 | 0.232 |
| Small | 200 | Unweighted UniFrac | 16 | 1.432 | 4 | 0.982 | 12 | 0.450 |
| Medium | 50 | Euclidean | 16 | 1.307 | 32 | 1.106 | -16 | 0.201 |
| Medium | 50 | Bray-Curtis | 16 | 1.307 | 74 | 1.325 | -58 | -0.018 |
| Medium | 50 | Jaccard | 16 | 1.307 | 0 | 0.961 | 16 | 0.346 |
| Medium | 50 | Weighted UniFrac | 16 | 1.307 | 32 | 1.320 | -16 | -0.013 |
| Medium | 50 | Unweighted UniFrac | 16 | 1.307 | 0 | 0.978 | 16 | 0.329 |
| Medium | 100 | Euclidean | 50 | 1.504 | 50 | 1.156 | 0 | 0.348 |
| Medium | 100 | Bray-Curtis | 50 | 1.504 | 90 | 1.400 | -40 | 0.104 |
| Medium | 100 | Jaccard | 50 | 1.504 | 0 | 0.978 | 50 | 0.526 |
| Medium | 100 | Weighted UniFrac | 50 | 1.504 | 60 | 1.404 | -10 | 0.100 |
| Medium | 100 | Unweighted UniFrac | 50 | 1.504 | 0 | 0.978 | 50 | 0.526 |
| Medium | 200 | Euclidean | 96 | 1.780 | 96 | 1.250 | 0 | 0.530 |
| Medium | 200 | Bray-Curtis | 96 | 1.780 | 100 | 1.540 | -4 | 0.240 |
| Medium | 200 | Jaccard | 96 | 1.780 | 0 | 0.978 | 96 | 0.802 |
| Medium | 200 | Weighted UniFrac | 96 | 1.780 | 98 | 1.540 | -2 | 0.240 |
| Medium | 200 | Unweighted UniFrac | 96 | 1.780 | 0 | 0.978 | 96 | 0.802 |
| Large | 50 | Euclidean | 84 | 1.585 | 90 | 1.280 | -6 | 0.305 |
| Large | 50 | Bray-Curtis | 84 | 1.585 | 100 | 1.560 | -16 | 0.025 |
| Large | 50 | Jaccard | 84 | 1.585 | 0 | 0.978 | 84 | 0.607 |
| Large | 50 | Weighted UniFrac | 84 | 1.585 | 94 | 1.560 | -10 | 0.025 |
| Large | 50 | Unweighted UniFrac | 84 | 1.585 | 0 | 0.978 | 84 | 0.607 |
| Large | 100 | Euclidean | 100 | 2.129 | 100 | 1.480 | 0 | 0.649 |
| Large | 100 | Bray-Curtis | 100 | 2.129 | 100 | 1.800 | 0 | 0.329 |
| Large | 100 | Jaccard | 100 | 2.129 | 0 | 0.978 | 100 | 1.151 |
| Large | 100 | Weighted UniFrac | 100 | 2.129 | 100 | 1.800 | 0 | 0.329 |
| Large | 100 | Unweighted UniFrac | 100 | 2.129 | 0 | 0.978 | 100 | 1.151 |
| Large | 200 | Euclidean | 100 | 3.129 | 100 | 1.780 | 0 | 1.349 |
| Large | 200 | Bray-Curtis | 100 | 3.129 | 100 | 2.200 | 0 | 0.929 |
| Large | 200 | Jaccard | 100 | 3.129 | 0 | 0.978 | 100 | 2.151 |
| Large | 200 | Weighted UniFrac | 100 | 3.129 | 100 | 2.200 | 0 | 0.929 |
| Large | 200 | Unweighted UniFrac | 100 | 3.129 | 0 | 0.978 | 100 | 2.151 |

**Abbreviations:** Power, empirical statistical power (percentage of simulations with p < 0.05); F, PERMANOVA F-statistic (mean across 50 simulations per condition); Power Difference, MeLSI Power - Traditional Power (%); F Difference, MeLSI Mean F - Traditional Mean F. Results based on 50 simulations per condition.

**Note:** The five traditional methods are: (1) **Euclidean distance** - standard Euclidean distance on CLR-transformed data; (2) **Bray-Curtis dissimilarity** - count-based dissimilarity metric; (3) **Jaccard dissimilarity** - binary (presence/absence) dissimilarity; (4) **Weighted UniFrac** - phylogenetically-informed distance using abundance-weighted branch lengths; (5) **Unweighted UniFrac** - phylogenetically-informed distance using presence/absence of taxa.

---

## Supplementary Table S3: Individual Method Comparisons for Scalability Analysis

This table will be generated after Table 3 job completes and will show detailed comparisons between MeLSI and each traditional method for scalability analysis.

*[To be added after Table 3 simulations complete]*

---

## Supplementary Table S4: Individual Method Comparisons for Feature Correlation Analysis

This table shows detailed comparisons between MeLSI and each of the five traditional methods individually across all correlation levels. These comparisons support the rank calculations shown in Table 6.

| Correlation Level | Correlation Value | Traditional Method | MeLSI Power (%) | MeLSI Mean F | Traditional Power (%) | Traditional Mean F | Power Difference (%) | F Difference |
|-------------------|-------------------|-------------------|-----------------|--------------|----------------------|-------------------|---------------------|--------------|
| None | 0.0 | Euclidean | 50.0 | 1.512 | 68.0 | 1.218 | -18.0 | 0.294 |
| None | 0.0 | Bray-Curtis | 50.0 | 1.512 | 100.0 | 1.619 | -50.0 | -0.107 |
| None | 0.0 | Jaccard | 50.0 | 1.512 | 8.0 | 1.024 | 42.0 | 0.488 |
| None | 0.0 | Weighted UniFrac | 50.0 | 1.512 | 64.0 | 1.671 | -14.0 | -0.158 |
| None | 0.0 | Unweighted UniFrac | 50.0 | 1.512 | 6.0 | 1.052 | 44.0 | 0.461 |
| Low | 0.3 | Euclidean | 42.0 | 1.481 | 58.0 | 1.190 | -16.0 | 0.291 |
| Low | 0.3 | Bray-Curtis | 42.0 | 1.481 | 96.0 | 1.500 | -54.0 | -0.019 |
| Low | 0.3 | Jaccard | 42.0 | 1.481 | 6.0 | 0.998 | 36.0 | 0.483 |
| Low | 0.3 | Weighted UniFrac | 42.0 | 1.481 | 40.0 | 1.498 | 2.0 | -0.018 |
| Low | 0.3 | Unweighted UniFrac | 42.0 | 1.481 | 6.0 | 0.992 | 36.0 | 0.489 |
| Moderate | 0.6 | Euclidean | 46.0 | 1.498 | 54.0 | 1.205 | -8.0 | 0.293 |
| Moderate | 0.6 | Bray-Curtis | 46.0 | 1.498 | 100.0 | 1.513 | -54.0 | -0.015 |
| Moderate | 0.6 | Jaccard | 46.0 | 1.498 | 6.0 | 1.000 | 40.0 | 0.498 |
| Moderate | 0.6 | Weighted UniFrac | 46.0 | 1.498 | 32.0 | 1.385 | 14.0 | 0.113 |
| Moderate | 0.6 | Unweighted UniFrac | 46.0 | 1.498 | 10.0 | 0.983 | 36.0 | 0.515 |
| High | 0.8 | Euclidean | 44.0 | 1.507 | 52.0 | 1.208 | -8.0 | 0.299 |
| High | 0.8 | Bray-Curtis | 44.0 | 1.507 | 96.0 | 1.492 | -52.0 | 0.015 |
| High | 0.8 | Jaccard | 44.0 | 1.507 | 6.0 | 1.051 | 38.0 | 0.456 |
| High | 0.8 | Weighted UniFrac | 44.0 | 1.507 | 20.0 | 1.376 | 24.0 | 0.131 |
| High | 0.8 | Unweighted UniFrac | 44.0 | 1.507 | 8.0 | 1.021 | 36.0 | 0.486 |

**Abbreviations:** Power, empirical statistical power (percentage of simulations with p < 0.05); F, PERMANOVA F-statistic (mean across 50 simulations per correlation level); Power Difference, MeLSI Power - Traditional Power (%); F Difference, MeLSI Mean F - Traditional Mean F. Results based on 50 simulations per correlation level.

**Note:** The five traditional methods are: (1) **Euclidean distance** - standard Euclidean distance on CLR-transformed data; (2) **Bray-Curtis dissimilarity** - count-based dissimilarity metric; (3) **Jaccard dissimilarity** - binary (presence/absence) dissimilarity; (4) **Weighted UniFrac** - phylogenetically-informed distance using abundance-weighted branch lengths; (5) **Unweighted UniFrac** - phylogenetically-informed distance using presence/absence of taxa.
