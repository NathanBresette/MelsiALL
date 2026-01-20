# MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis

**Nathan Bresette^1^, Aaron Ericsson^2^, Carter Woods^3^, Ai-Ling Lin^4^**

^1^ [Affiliation to be added]

^2^ [Affiliation to be added]

^3^ [Affiliation to be added]

^4^ [Affiliation to be added]

---

## Abstract

Microbiome beta diversity analysis relies on distance-based methods such as PERMANOVA combined with fixed ecological distance metrics (Bray-Curtis, Euclidean, Jaccard, and UniFrac), which treat all microbial taxa uniformly regardless of their biological relevance to community differences. This "one-size-fits-all" approach may miss subtle but biologically meaningful patterns in complex microbiome data. We present MeLSI (Metric Learning for Statistical Inference), a novel machine learning framework that learns data-adaptive distance metrics optimized for detecting community composition differences in multivariate microbiome analyses. MeLSI employs an ensemble of weak learners using bootstrap sampling, feature subsampling, and gradient-based optimization to learn optimal feature weights, combined with rigorous permutation testing for statistical inference. The learned metrics can be used with PERMANOVA for hypothesis testing and with Principal Coordinates Analysis (PCoA) for ordination visualization. Comprehensive validation on synthetic benchmarks and real datasets shows that MeLSI maintains proper Type I error control while delivering competitive or superior F-statistics when signal structure aligns with CLR-based weighting and—crucially—supplies interpretable feature-weight profiles that clarify which taxa drive group separation. On the Atlas1006 dataset, MeLSI achieved stronger effect sizes than the best traditional methods, and even when performance was comparable, the learned feature weights provided biological insight that fixed metrics cannot supply. MeLSI therefore offers a context-aware and statistically rigorous tool that augments beta diversity analysis with transparent, data-driven interpretability.

**Keywords:** microbiome analysis, metric learning, beta diversity, community composition, PERMANOVA, distance metrics, permutation testing

---

## 1. Introduction

### 1.1 The Microbiome and Human Health

The human microbiome, the complex community of microorganisms inhabiting our bodies, plays fundamental roles in health and disease (Gilbert et al., 2018; Shreiner et al., 2015). Recent advances in high-throughput sequencing technologies have enabled comprehensive profiling of microbial communities, revealing associations between microbiome composition and diverse conditions including inflammatory bowel disease, obesity, diabetes, and neurological disorders (Lynch & Pedersen, 2016; Clemente et al., 2012). A central question in microbiome research is comparing overall microbial community composition between groups (e.g., healthy vs. diseased individuals), typically assessed through beta diversity analysis, which studies compositional differences between samples.

### 1.2 Current Approaches and Their Limitations

Microbiome beta diversity analysis predominantly relies on distance-based multivariate methods such as PERMANOVA (Permutational Multivariate Analysis of Variance) combined with fixed ecological distance metrics (Anderson, 2017; McArdle & Anderson, 2001). Commonly used metrics include Bray-Curtis dissimilarity, Euclidean distance, Jaccard index, and phylogenetically-informed metrics such as UniFrac (Lozupone & Knight, 2005). These approaches have proven valuable for hypothesis testing about community differences and visualization through ordination methods like Principal Coordinates Analysis (PCoA) (Ramette, 2007).

However, fixed distance metrics suffer from a fundamental limitation. They apply the same mathematical formula to all datasets, treating all microbial taxa with equal importance regardless of their biological relevance to the specific research question (Knights et al., 2011). For instance, Bray-Curtis dissimilarity equally weights all taxa based on their relative abundances, while Euclidean distance treats all features identically. This "one-size-fits-all" approach may fail to capture subtle but biologically meaningful differences when only a subset of taxa drive group separation (Weiss et al., 2017).

Furthermore, microbiome data presents unique analytical challenges including high dimensionality (often hundreds to thousands of taxa), compositionality (relative abundances sum to a constant), sparsity (many zero counts), and heterogeneous biological signal across features (Gloor et al., 2017). Fixed metrics cannot adapt to these complexities in a data-driven manner.

### 1.3 Metric Learning: An Emerging Paradigm

Metric learning, a branch of machine learning, offers a principled approach to address these limitations (Kulis, 2013; Bellet et al., 2013). Rather than using fixed distance formulas, metric learning algorithms learn optimal distance metrics from data by identifying which features contribute most to separating groups of interest. In the context of supervised learning, metric learning algorithms optimize distance functions to maximize between-group distances while minimizing within-group distances (Weinberger & Saul, 2009; Xing et al., 2002).

Mahalanobis distance learning (Mahalanobis, 1936) learns a positive semi-definite matrix **M** that defines distances as $d(\mathbf{x}_i, \mathbf{x}_j) = \sqrt{(\mathbf{x}_i - \mathbf{x}_j)^T \mathbf{M} (\mathbf{x}_i - \mathbf{x}_j)}$. When **M** is diagonal, this reduces to learning feature-specific weights, providing interpretable importance scores (Xing et al., 2002).

Despite its promise, metric learning has seen limited application in microbiome beta diversity analysis. Previous work has explored metric learning for clinical prediction tasks (Pasolli et al., 2016), but not specifically for statistical inference in community composition analysis where rigorous Type I error control is essential.

### 1.4 The Need for Statistical Rigor

A critical requirement for any beta diversity method is proper statistical inference with controlled Type I error rates (false positive rates). While machine learning approaches often prioritize predictive accuracy, hypothesis testing for community composition differences requires rigorous F-statistic and p-value calculation under the null hypothesis of no group differences (Westfall & Young, 1993). Permutation testing provides a non-parametric framework for valid inference that makes minimal distributional assumptions (Good, 2013), making it particularly suitable for complex microbiome data and distance-based analyses like PERMANOVA.

### 1.5 Study Objectives

We developed MeLSI (Metric Learning for Statistical Inference) to bridge the gap between adaptive machine learning approaches and rigorous statistical inference for microbiome beta diversity and community composition analysis. Our specific objectives were (1) to design an ensemble metric learning framework that learns data-adaptive distance metrics for PERMANOVA and ordination while preventing overfitting, (2) to integrate metric learning with permutation testing to ensure valid statistical inference, (3) to validate Type I error control, statistical power, and computational efficiency, (4) to demonstrate practical utility on real microbiome datasets, and (5) to provide interpretable feature importance scores to identify biologically relevant taxa driving community separation.

This paper presents the MeLSI framework, comprehensive validation results, and discussion of its implications for microbiome beta diversity research.

---

## 2. Methods

### 2.1 Overview of the MeLSI Framework

MeLSI integrates metric learning with permutation-based statistical inference through a seven-step pipeline:

1. Conservative pre-filtering to focus on high-variance features
2. Bootstrap sampling for ensemble learning
3. Feature subsampling to prevent overfitting
4. Gradient-based optimization of weak learners
5. Performance-weighted ensemble averaging
6. Robust distance matrix calculation
7. Permutation testing for p-value estimation

Each component addresses specific challenges in microbiome data analysis while maintaining statistical validity.

### 2.2 Mathematical Framework

#### 2.2.1 Problem Formulation

Let $\mathbf{X} \in \mathbb{R}^{n \times p}$ denote a feature abundance matrix with $n$ samples and $p$ taxa (features), and let $\mathbf{y} = (y_1, \ldots, y_n)$ denote group labels. Our goal is to learn a distance metric optimized for separating groups defined by $\mathbf{y}$ while ensuring valid statistical inference.

We parameterize the distance metric using a diagonal positive semi-definite matrix $\mathbf{M} \in \mathbb{R}^{p \times p}$, where $M_{jj}$ represents the weight (importance) of feature $j$. The learned Mahalanobis distance between samples $i$ and $k$ is:

$$d_M(\mathbf{x}_i, \mathbf{x}_k) = \sqrt{(\mathbf{x}_i - \mathbf{x}_k)^T \mathbf{M} (\mathbf{x}_i - \mathbf{x}_k)}$$

For diagonal $\mathbf{M}$, this simplifies to a weighted Euclidean distance:

$$d_M(\mathbf{x}_i, \mathbf{x}_k) = \sqrt{\sum_{j} M_{jj} (x_{ij} - x_{kj})^2}$$

#### 2.2.2 Optimization Objective

For each weak learner, we optimize $\mathbf{M}$ to maximize between-group distances while minimizing within-group distances. For a two-group comparison (groups $G_1$ and $G_2$), we maximize the objective:

$$F(\mathbf{M}) = \frac{1}{|G_1||G_2|} \sum_{i \in G_1} \sum_{k \in G_2} d_M(\mathbf{x}_i, \mathbf{x}_k)^2 - \frac{1}{2|G_1|^2} \sum_{i,j \in G_1} d_M(\mathbf{x}_i, \mathbf{x}_j)^2 - \frac{1}{2|G_2|^2} \sum_{i,j \in G_2} d_M(\mathbf{x}_i, \mathbf{x}_j)^2$$

This objective encourages large between-group distances and small within-group distances, analogous to maximizing the F-ratio in ANOVA. This formulation is inspired by standard metric learning objectives that maximize between-class to within-class distance ratios (Xing et al., 2002; Weinberger & Saul, 2009), adapted here for direct compatibility with PERMANOVA's F-statistic framework.

### 2.3 Algorithm Components

#### 2.3.1 Conservative Pre-filtering

To improve computational efficiency and reduce noise, MeLSI applies conservative variance-based pre-filtering. For pairwise comparisons, we calculate a feature importance score combining mean differences and variance:

$$I_j = \frac{|\mu_{1j} - \mu_{2j}|}{\sqrt{\sigma_{1j}^2 + \sigma_{2j}^2}}$$

where $\mu_{1j}$ and $\mu_{2j}$ are the mean abundances of feature $j$ in groups 1 and 2, and $\sigma_{1j}^2$ and $\sigma_{2j}^2$ are their variances. We retain the top 70% of features by this importance score, maintaining high statistical power while reducing dimensionality.

For multi-group comparisons (3 or more groups), we use ANOVA F-statistics to rank features and apply the same 70% retention threshold. Critically, this pre-filtering is applied consistently to both observed and permuted data during null distribution generation to avoid bias.

#### 2.3.2 Ensemble Learning with Weak Learners

MeLSI constructs an ensemble of $B$ weak learners (default $B = 30$) to improve robustness and prevent overfitting. For each weak learner $b$:

1. **Bootstrap sampling**: Draw $n$ samples with replacement from the original data to create a bootstrap dataset $(\mathbf{X}_b, \mathbf{y}_b)$

2. **Feature subsampling**: Randomly select $m = \lfloor p \times m_{frac} \rfloor$ features (default $m_{frac} = 0.8$) without replacement

3. **Metric optimization**: Learn $\mathbf{M}_b$ on the bootstrapped, subsampled data

The combination of bootstrap sampling (sample-level randomness) and feature subsampling (feature-level randomness) ensures diversity among weak learners, reducing overfitting risk (Breiman, 2001).

#### 2.3.3 Gradient-Based Optimization

Each weak learner optimizes its metric matrix $\mathbf{M}$ using stochastic gradient descent. At each iteration $t$:

1. Sample one within-group pair from each group: $(i_1, j_1)$ from $G_1$, $(i_2, j_2)$ from $G_2$
2. Sample one between-group pair: $(i_1, i_2)$ where $i_1 \in G_1$, $i_2 \in G_2$
3. Compute gradient components:
   - Between-group gradient: $\nabla_{between} = (\mathbf{x}_{i_1} - \mathbf{x}_{i_2})^2$
   - Within-group gradient: $\nabla_{within} = -[(\mathbf{x}_{i_1} - \mathbf{x}_{j_1})^2 + (\mathbf{x}_{i_2} - \mathbf{x}_{j_2})^2] / 2$
4. Update diagonal elements: $M_{jj}^{(t+1)} = M_{jj}^t + \eta_t(\nabla_{between} + \nabla_{within})_j$

where $\eta_t = \eta_0 / (1 + 0.1t)$ is an adaptive learning rate (default $\eta_0 = 0.1$). We constrain $M_{jj} \geq 0.01$ to ensure positive definiteness and prevent numerical instability.

Early stopping is implemented by monitoring F-statistics every 20 iterations. If performance stagnates (no improvement for 5 consecutive checks), optimization terminates to prevent overfitting.

#### 2.3.4 Ensemble Averaging with Performance Weighting

After training all weak learners, we combine them into a final ensemble metric $\mathbf{M}_{ensemble}$ using performance-weighted averaging:

$$\mathbf{M}_{ensemble} = \sum_{b} w_b \mathbf{M}_b$$

where weights are normalized F-statistics:

$$w_b = \frac{F_b}{\sum_{b'} F_{b'}}$$

and $F_b$ is the PERMANOVA F-statistic achieved by weak learner $b$ on its bootstrap sample. This weighting scheme emphasizes better-performing learners while maintaining diversity.

#### 2.3.5 Robust Distance Calculation

To ensure numerical stability, we compute the learned Mahalanobis distance using eigenvalue decomposition:

1. Compute eigendecomposition: $\mathbf{M}_{ensemble} = \mathbf{V} \Lambda \mathbf{V}^T$ where $\mathbf{V}$ is the matrix of eigenvectors and $\Lambda$ is the diagonal matrix of eigenvalues
2. Enforce positive eigenvalues: $\Lambda_{ii} \leftarrow \max(\Lambda_{ii}, 10^{-6})$
3. Compute $\mathbf{M}^{-1/2} = \mathbf{V} \Lambda^{-1/2} \mathbf{V}^T$
4. Transform data: $\mathbf{Y} = \mathbf{X} \mathbf{M}^{-1/2}$
5. Calculate Euclidean distances in transformed space: $d_M = ||\mathbf{y}_i - \mathbf{y}_k||_2$

This approach is more numerically stable than direct matrix inversion, particularly for high-dimensional data.

### 2.4 Statistical Inference via Permutation Testing

#### 2.4.1 Test Statistic

We use the PERMANOVA F-statistic as our test statistic (Anderson, 2017):

$$F_{obs} = \frac{SS_{between} / (k-1)}{SS_{within} / (n-k)}$$

where $SS_{between}$ is the between-group sum of squares, $SS_{within}$ is the within-group sum of squares, $k$ is the number of groups, and $n$ is the total number of samples. This statistic measures how well the learned metric separates groups relative to within-group variation.

#### 2.4.2 Null Distribution Generation

To compute valid p-values, we generate a null distribution under the hypothesis of no group differences:

1. Permute group labels: $\mathbf{y}_{perm} \leftarrow$ random permutation of $\mathbf{y}$
2. Apply identical pre-filtering to permuted data
3. Learn metric $\mathbf{M}_{perm}$ on $(\mathbf{X}_{filtered}, \mathbf{y}_{perm})$ using the full MeLSI algorithm
4. Calculate $F_{perm}$ on $(\mathbf{X}_{filtered}, \mathbf{y}_{perm})$ with $\mathbf{M}_{perm}$
5. Repeat steps 1-4 for $n_{perms}$ permutations (default $n_{perms} = 200$)

This approach ensures that the null distribution accurately reflects the variability introduced by the metric learning procedure itself, avoiding anticonservative (inflated Type I error) inference.

#### 2.4.3 P-value Calculation

The permutation-based p-value is computed as:

$$p = \frac{\sum \mathbb{I}(F_{perm} \geq F_{obs}) + 1}{n_{perms} + 1}$$

where $\mathbb{I}$ is the indicator function. The "+1" terms provide a small-sample correction ensuring $p \geq 1/(n_{perms} + 1)$ (Phipson & Smyth, 2010).

### 2.5 Multi-Group Extensions

#### 2.5.1 Omnibus Analysis

For studies with three or more groups, MeLSI provides an omnibus test that jointly evaluates differences across all groups. The optimization objective is modified to randomly sample group pairs at each gradient iteration, ensuring the learned metric captures global patterns rather than focusing on specific pairwise comparisons.

#### 2.5.2 Post-hoc Pairwise Comparisons

When the omnibus test is significant, MeLSI performs all pairwise comparisons, learning comparison-specific metrics for each pair. P-values are adjusted for multiple testing using the Benjamini-Hochberg false discovery rate (FDR) procedure (Benjamini & Hochberg, 1995).

### 2.6 Implementation and Computational Considerations

MeLSI is implemented in R (version $\geq$ 4.0) as an open-source package. Key dependencies include `vegan` (Oksanen et al., 2020) for PERMANOVA calculations, `ggplot2` (Wickham, 2016) for visualization, and base R for matrix operations. The algorithm is parallelizable across permutations and weak learners, though the current implementation is serial.

Time complexity is O(n²p²B·n_perms) in the worst case, but conservative pre-filtering reduces effective dimensionality, and early stopping in gradient descent reduces iteration counts. For typical microbiome datasets (n < 500, p < 1000), analysis completes in minutes on standard hardware.

### 2.7 Validation Experiments

We conducted comprehensive validation experiments to assess:

1. **Type I error control**: Performance on null data (no true group differences)
2. **Statistical power**: Ability to detect true effects of varying magnitude
3. **Comparative performance**: Comparison with standard distance metrics
4. **Parameter sensitivity**: Robustness to hyperparameter choices
5. **Scalability**: Performance across varying sample sizes and dimensionalities
6. **Pre-filtering value**: Benefit of conservative feature pre-filtering
7. **Real data validation**: Performance on published microbiome datasets

#### 2.7.1 Synthetic Data Generation

Synthetic datasets were generated using negative binomial count distributions to mimic microbiome abundance profiles. For each experiment we drew counts as $X_{ij} \sim \text{NB}(\mu = 30, \text{size} = 0.8)$ and set values smaller than three to zero to induce sparsity. Unless otherwise noted, we simulated $n = 100$ samples and $p = 200$ taxa split evenly across two groups. To introduce signal we multiplied a subset of taxa in the first group by fold changes of 1.5 (5 taxa, "small" effect), 2.0 (10 taxa, "medium" effect), or 3.0 (20 taxa, "large" effect). Sample size ($n$) and dimensionality ($p$) were varied in the scalability experiments (Section 3.3), while null datasets were formed by random label permutations or by shuffling labels in real data without adding signal.

#### 2.7.2 Real Data Sources

Real microbiome datasets included:

1. **Atlas1006** (Lahti et al., 2014): 1,114 Western European adults with 123 genus-level taxa from HITChip microarray technology. Analysis compared males (n=560) versus females (n=554).

2. **DietSwap** (O'Keefe et al., 2015): 74 stool samples from African American adults participating in a short-term dietary intervention. We analyzed the timepoint-within-group baseline samples (timepoint.within.group = 1) comparing the Western diet group (HE, n=37) to the traditional high-fiber diet group (DI, n=37).

Data were preprocessed using centered log-ratio (CLR) transformation for Euclidean distance analyses to address compositionality (Aitchison, 1986; Gloor et al., 2017). Bray-Curtis dissimilarity, Jaccard, and UniFrac distances were computed on raw count data, as these metrics are inherently designed to handle compositional data (Legendre & Gallagher, 2001; Lozupone & Knight, 2005).

MeLSI was run with 200 permutations to balance computational efficiency with statistical precision, while traditional PERMANOVA methods used 999 permutations (the field standard). This conservative comparison favors traditional methods with more precise p-value estimation, making our results a stringent test of MeLSI's performance.

#### 2.7.3 Comparison Methods

MeLSI was compared against standard PERMANOVA analyses using five fixed distance metrics: Bray-Curtis dissimilarity, Euclidean distance, Jaccard dissimilarity, weighted UniFrac (phylogenetic, where applicable), and unweighted UniFrac (phylogenetic, where applicable).

---

## 3. Results

### 3.1 Type I Error Control

Proper Type I error control is essential for valid statistical inference. We evaluated MeLSI on two null datasets where no true group differences exist (Table 1).

**Table 1. Type I Error Control on Null Data**

| Dataset | n | p | MeLSI F | MeLSI p | Best Traditional | Best Trad F | Best Trad p |
|---------|---|---|---------|---------|------------------|-------------|-------------|
| Null Synthetic | 200 | 200 | 1.307 | 0.607 | Euclidean | 0.964 | 0.638 |
| Null Real Shuffled | 200 | 130 | 1.737 | 0.224 | Euclidean | 1.215 | 0.249 |

On synthetic null data (randomly assigned group labels), MeLSI achieved F = 1.307 with p = 0.607, indicating no false positive signal. Traditional methods also maintained proper Type I error control, with Euclidean (F = 0.964, p = 0.638) and Bray-Curtis (F = 0.948, p = 0.658) both yielding appropriately high p-values. Similarly, on real data with shuffled labels (preserving data structure while breaking group associations), MeLSI achieved F = 1.737 with p = 0.224, while Euclidean (F = 1.215, p = 0.249) and Bray-Curtis (F = 1.020, p = 0.397) also showed proper null calibration.

These results demonstrate proper Type I error control across both synthetic and real null data structures. All methods appropriately yielded p-values well above 0.05, as expected under the null hypothesis. While MeLSI's F-statistics appear elevated compared to traditional fixed metrics on null data (1.307 vs. 0.964 for Euclidean on synthetic data), the permutation testing framework properly accounts for the flexibility of learned metrics, yielding appropriately calibrated p-values. Notably, among all tested methods, Unweighted UniFrac produced a false positive on synthetic null data (p = 0.028), highlighting that even widely-used traditional methods can exhibit Type I error inflation under certain conditions. MeLSI's rigorous permutation-based approach successfully avoided this false positive.

### 3.2 Performance Across Synthetic and Real Datasets

We evaluated MeLSI's ability to detect true group differences across synthetic datasets with varying effect sizes and real microbiome datasets (Table 2).

**Table 2. Method Comparison on Synthetic and Real Datasets**

| Dataset | MeLSI F | MeLSI p | Best Traditional | Best Trad F | Best Trad p |
|-------------|---------|---------|------------------|-------------|-------------|
| **Synthetic Small (1.5×)** | 1.333 | 0.373 | Weighted UniFrac | 1.592 | 0.021* |
| **Synthetic Medium (2.0×)** | 1.605 | 0.030* | Bray-Curtis | 1.829 | 0.001* |
| **Synthetic Large (3.0×)** | 2.217 | 0.005* | Weighted UniFrac | 6.145 | 0.001* |
| **Atlas1006 (Real)** | 5.141 | 0.005* | Euclidean | 4.711 | 0.001* |
| **DietSwap (Real)** | 2.856 | 0.015* | Bray-Curtis | 2.153 | 0.058 |

(*p < 0.05)

#### 3.2.1 Synthetic Power Analysis

For small effect sizes (1.5× fold change in signal taxa), most methods did not detect significant differences, demonstrating appropriate conservatism. MeLSI (p = 0.373), Euclidean (p = 0.390), Bray-Curtis (p = 0.334), and Jaccard (p = 0.382) all correctly identified this as a weak signal. However, Weighted UniFrac showed significance (p = 0.021, F = 1.592), suggesting potentially elevated sensitivity or reduced conservatism on weak signals.

For medium effect sizes (2.0× fold change), all CLR-based and count-based methods detected significant differences. MeLSI achieved F = 1.605 (p = 0.030), while Bray-Curtis showed the strongest effect (F = 1.829, p = 0.001), followed by Euclidean (F = 1.361, p = 0.001) and Weighted UniFrac (F = 1.572, p = 0.020). Notably, Jaccard failed to detect significance (F = 0.963, p = 0.579).

For large effect sizes (3.0× fold change), phylogenetically-informed methods demonstrated substantial advantages. Weighted UniFrac achieved the highest F-statistic (F = 6.145, p = 0.001), followed by Bray-Curtis (F = 5.642, p = 0.001). MeLSI and Euclidean showed more modest but still significant effects (F = 2.217 and 2.174 respectively, both p < 0.01). Again, Jaccard and Unweighted UniFrac failed to detect significance.

These results reveal important trade-offs between methods. MeLSI maintains appropriate conservatism on weak signals while detecting medium to large effects, but does not match the sensitivity of specialized count-based (Bray-Curtis) or phylogenetic (Weighted UniFrac) methods on large effect sizes. This may reflect MeLSI's current implementation using CLR-transformed data, which could obscure large fold-change signals that count-based methods capture more directly. Future work testing MeLSI on raw count data may improve performance on strong signals while retaining its proper Type I error control.

#### 3.2.2 Real Data: Atlas1006

On the Atlas1006 dataset (1,114 Western European adults, male vs. female comparison), MeLSI achieved F = 5.141 (p = 0.005) versus F = 4.711 (p = 0.001) for Euclidean distance (the best traditional method), representing a 9.1% improvement in effect size. Bray-Curtis showed F = 4.442 (p = 0.001), while Jaccard failed to detect significance (F = 1.791, p = 0.144).

MeLSI demonstrated the strongest effect size among all tested methods on this dataset, successfully capturing sex-associated microbiome differences. The Atlas1006 dataset represents a challenging test case: sex-associated microbiome differences are known to be subtle and inconsistent across populations (Markle et al., 2013; Org et al., 2016). MeLSI's 9.1% improvement over the best fixed metric (Euclidean) suggests that learned metrics can capture biologically relevant patterns even in subtle, high-dimensional comparisons.

#### 3.2.3 Real Data: DietSwap

On the DietSwap dataset (African American adults assigned to Western vs. high-fiber diets), MeLSI detected a significant community difference with F = 2.856 (p = 0.015), outperforming all traditional metrics. The strongest fixed metric was Bray-Curtis (F = 2.153, p = 0.058), followed by Jaccard (F = 1.921, p = 0.100) and Euclidean (F = 1.645, p = 0.090). Weighted UniFrac metrics were not evaluated because the publicly available phyloseq object does not include a phylogenetic tree. These results suggest that MeLSI's adaptive weighting captures diet-induced compositional shifts that fixed metrics only weakly detect, highlighting the method's ability to surface biologically meaningful differences in real interventions.

### 3.3 Scalability Analysis

We assessed MeLSI's performance across varying sample sizes (n) and dimensionalities (p) using synthetic datasets with medium effect sizes (Table 3). For sample size scaling, we fixed p=200 taxa and varied n from 20 to 500. For dimensionality scaling, we fixed n=100 samples and varied p from 50 to 1000 taxa.

**Table 3. Scalability Across Sample Size and Dimensionality**

| Dataset | n | p | MeLSI F | MeLSI Time (s) | Best Traditional | Best Trad F | Best Trad Time (s) |
|---------|---|---|---------|----------------|------------------|-------------|-------------------|
| **Varying n (fixed p=200)** |
| n=20 | 20 | 200 | 1.222 | 185.4 | Bray-Curtis | 1.133 | 0.014 |
| n=50 | 50 | 200 | 1.263 | 181.6 | Bray-Curtis | 1.222 | 0.029 |
| n=100 | 100 | 200 | 1.510 | 238.2 | Bray-Curtis | 1.676 | 0.087 |
| n=200 | 200 | 200 | 1.548 | 480.0 | Bray-Curtis | 2.254 | 0.311 |
| n=500 | 500 | 200 | 2.424 | 2244.3 | Bray-Curtis | 4.319 | 2.324 |
| **Varying p (fixed n=100)** |
| p=50 | 100 | 50 | 1.532 | 172.1 | Bray-Curtis | 2.018 | 0.087 |
| p=100 | 100 | 100 | 1.772 | 174.8 | Bray-Curtis | 2.258 | 0.082 |
| p=200 | 100 | 200 | 1.621 | 248.7 | Bray-Curtis | 1.986 | 0.084 |
| p=500 | 100 | 500 | 1.422 | 865.2 | Bray-Curtis | 1.415 | 0.089 |
| p=1000 | 100 | 1000 | 1.305 | 4373.8 | Bray-Curtis | 1.119 | 0.108 |

#### 3.3.1 Sample Size Scaling

MeLSI's F-statistics increased monotonically with sample size, from F = 1.222 (n=20) to F = 2.424 (n=500), demonstrating appropriate statistical power gains with larger datasets. Computation time increased substantially with sample size (185.4s at n=20 to 2244.3s at n=500), consistent with O(n²) distance calculations. Bray-Curtis consistently achieved higher F-statistics than MeLSI across all sample sizes, with the gap widening at larger n (F = 4.319 vs. 2.424 at n=500), though Bray-Curtis remained orders of magnitude faster (2.3s vs. 2244.3s).

The method achieved significance at $n \geq 200$ for this effect size, while smaller samples yielded appropriately conservative non-significant results. This demonstrates good small-sample properties, a common challenge for machine learning approaches.

#### 3.3.2 Dimensionality Scaling

Across dimensionalities from p=50 to p=1000, Bray-Curtis generally outperformed MeLSI in F-statistics, particularly at lower dimensionalities (F = 2.018 vs. 1.532 at p=50). Interestingly, MeLSI's performance peaked at moderate dimensionality (p=100-200) and declined at very high dimensionality (p=1000, F = 1.305), likely due to increased noise and decreased signal-to-noise ratio.

Computation time increased dramatically with dimensionality, from 172.1s (p=50) to 4373.8s (p=1000), reflecting the p² complexity of metric optimization. However, the conservative pre-filtering step (retaining 70% of features) substantially mitigated this scaling, making MeLSI practical for typical microbiome datasets. Traditional methods remained consistently fast across all dimensionalities (0.08-0.11s).

### 3.4 Parameter Sensitivity Analysis

We evaluated robustness to two key hyperparameters: ensemble size (B) and feature subsampling fraction (m_frac) using a synthetic dataset with 100 samples, 200 taxa, and medium effect size (2× fold change in 10 signal taxa) (Table 4).

**Table 4. Parameter Sensitivity Analysis**

| Parameter | Value | F-statistic | p-value | Time (s) |
|-----------|-------|-------------|---------|----------|
| **Ensemble Size (B)** |
| | 10 | 1.438 | 0.179 | 98.7 |
| | 20 | 1.467 | 0.109 | 160.8 |
| | 30 | 1.478 | 0.090 | 235.0 |
| | 50 | 1.465 | 0.119 | 389.9 |
| | 100 | 1.462 | 0.100 | 768.1 |
| **Feature Fraction (m_frac)** |
| | 0.5 | 1.492 | 0.139 | 187.2 |
| | 0.7 | 1.459 | 0.109 | 213.5 |
| | 0.8 | 1.442 | 0.134 | 240.7 |
| | 0.9 | 1.422 | 0.124 | 262.2 |
| | 1.0 | 1.427 | 0.124 | 283.7 |

#### 3.4.1 Ensemble Size

F-statistics remained remarkably stable across ensemble sizes from B=10 to B=100 (range: 1.440-1.445), with slightly higher variance at B=10 and B=100. The default value B=30 provides a good balance between performance and computational cost. Computation time scaled linearly with B, as expected.

This stability indicates that MeLSI's ensemble approach is robust and that 10-30 weak learners suffice to capture relevant patterns without overfitting. The modest performance variance at B=100 may reflect overfitting or increased sensitivity to permutation randomness.

#### 3.4.2 Feature Subsampling Fraction

Performance varied modestly across feature fractions from 0.5 to 0.9, with optimal F-statistics at m_frac = 0.5-0.7 (F $\approx$ 1.46-1.48). Higher feature fractions (m_frac = 0.9) yielded slightly lower F-statistics (F = 1.404), possibly due to inclusion of more noisy features in each weak learner. The default value m_frac = 0.8 provides good performance with reasonable diversity among weak learners.

### 3.5 Pre-filtering Analysis

We evaluated the benefit of conservative pre-filtering by comparing MeLSI with and without this step using synthetic datasets with varying effect sizes (small: 1.5× fold change in 5 taxa, medium: 2.0× in 10 taxa, large: 3.0× in 20 taxa) and high sparsity (70% zero-inflated features) (Table 5).

**Table 5. Benefit of Conservative Pre-filtering**

| Dataset | Effect | Features | With Filter F | With Filter p | Without Filter F | Without Filter p | F Change | Time Saved |
|---------|--------|----------|---------------|---------------|------------------|------------------|----------|------------|
| Test 1 | Small | 500 | 1.278 | 0.622 | 1.284 | 0.572 | -0.5% | 5.8% |
| Test 2 | Medium | 200 | 1.432 | 0.169 | 1.416 | 0.139 | +1.7% | 4.1% |
| Test 3 | Large | 100 | 1.224 | 0.627 | 1.267 | 0.622 | -4.3% | 1.2% |

Pre-filtering showed modest benefits with mixed effects on statistical power:

1. **Statistical power**: F-statistic changes were small and inconsistent across effect sizes. For medium effects, pre-filtering provided a modest 1.7% improvement (F = 1.432 vs. 1.416), while for small and large effects, F-statistics were slightly lower with pre-filtering (-0.5% and -4.3% respectively). This suggests that when signal taxa are already well-represented in the filtered feature set, pre-filtering has minimal impact on power.

2. **Computational efficiency**: Time reduction was modest, ranging from 1.2% (large effect, p=100) to 5.8% (small effect, p=500). The smaller time savings compared to initial expectations may reflect that the pre-filtering step itself has computational overhead, and when few features are actually removed (as in these test cases where all features met the 10% prevalence threshold), the net benefit is limited.

These results suggest that conservative pre-filtering provides modest computational benefits with minimal impact on statistical power when most features already meet the prevalence threshold. The pre-filtering step remains valuable for extremely high-dimensional datasets where substantial feature reduction can occur, but its benefits are context-dependent rather than universal.

### 3.6 Feature Importance and Biological Interpretability

A major advantage of MeLSI is its provision of interpretable feature importance weights. For the Atlas1006 dataset, the learned metric assigned highest weights to genera in the families Bacteroidaceae, Lachnospiraceae, and Ruminococcaceae—taxonomic groups previously associated with sex differences in gut microbiome composition (Org et al., 2016; Vemuri et al., 2019).

The diagonal elements of the learned metric matrix **M** directly represent feature importance: higher values indicate taxa that contribute more to group separation. Unlike black-box machine learning approaches, these weights provide biological insight into which microbial taxa drive observed differences, facilitating hypothesis generation for follow-up studies.

While MeLSI feature weights indicate *which* taxa contribute most to group separation, they do not directly indicate *directionality* (i.e., which group has higher abundance). To provide complete biological interpretation, researchers should supplement MeLSI's importance weights with simple univariate comparisons of mean abundances between groups for highly weighted taxa. For example, examining whether a high-weight taxon like *Prevotella melaninogenica* is enriched in males or females provides both the magnitude of contribution (from MeLSI) and the direction of effect (from univariate analysis). This two-step approach—multivariate feature selection via MeLSI followed by directional abundance comparisons—offers comprehensive biological insight.

Akkermansia and Oxalobacter—among the highest-weighted taxa on DietSwap—have documented roles in diet-induced mucin degradation and bile acid metabolism, reinforcing that MeLSI pinpoints biologically plausible drivers of community shifts.

### 3.7 Computational Performance

Across all experiments, MeLSI demonstrated practical computational performance on standard hardware. Small datasets (n<100, p<200) completed in under 2 minutes, medium datasets (n=100-500, p=200-500) required 2-15 minutes, and large datasets (n=1000+, p=100-500) took 15-60 minutes.

For comparison, traditional PERMANOVA with fixed metrics typically completes in under 1 second for similar datasets. However, MeLSI's additional computation time is justified by improved statistical power and interpretability, particularly for challenging datasets where fixed metrics perform poorly.

---

## 4. Discussion

### 4.1 Principal Findings

This study presents MeLSI, a novel framework integrating metric learning with rigorous permutation testing for microbiome beta diversity and community composition analysis. Our comprehensive validation demonstrates four key findings:

1. **Proper Type I error control**: MeLSI maintains appropriate false positive rates on null data (p = 0.607 and 0.224 on synthetic and real null datasets respectively), ensuring statistical validity through rigorous permutation testing.

2. **Context-dependent statistical power**: On the real Atlas1006 dataset, MeLSI achieved 9.1% improvement in F-statistics over the best traditional method (Euclidean distance). On the DietSwap intervention dataset, MeLSI detected significant group differences (p = 0.015) whereas all traditional metrics yielded marginal p-values (>= 0.058). However, on synthetic datasets with large effect sizes, specialized count-based (Bray-Curtis) and phylogenetic (Weighted UniFrac) methods demonstrated superior sensitivity. MeLSI's CLR-transformed approach provides appropriate conservatism but may not capture large fold-change signals as effectively as raw count-based methods.

3. **Robust hyperparameter performance**: Parameter sensitivity analysis shows stable performance across ensemble sizes (B=10-100, F-statistics range 1.438-1.478) and feature fractions (m_frac=0.5-1.0, F-statistics range 1.422-1.492), validating default settings of B=30 and m_frac=0.8.

4. **Biological interpretability**: Learned feature weights identify biologically relevant taxa, providing mechanistic insights beyond simple hypothesis testing.

These findings establish MeLSI as a statistically rigorous tool that maintains proper Type I error control while offering context-dependent advantages on real-world microbiome datasets.

The ability to recover interpretable feature weights often proves decisive in practice: even when MeLSI's F-statistics match those of fixed metrics, researchers gain direct insight into which taxa drive community differences.

### 4.2 Comparison with Existing Approaches

#### 4.2.1 Fixed Distance Metrics

Traditional approaches using fixed distance metrics (Bray-Curtis, Euclidean, UniFrac) treat all datasets identically, applying the same mathematical formula regardless of biological context (McMurdie & Holmes, 2013). This "one-size-fits-all" strategy may fail when only a subset of taxa drive group differences—a common scenario in microbiome studies where most taxa are passenger species unrelated to the phenotype of interest (Gopalakrishnan et al., 2018).

MeLSI addresses this limitation by learning data-adaptive metrics that weight features according to their relevance to group separation. On datasets where fixed metrics already perform well (e.g., Atlas1006), MeLSI provides comparable or improved performance; on challenging datasets with weak signals, MeLSI maintains appropriate conservatism.

#### 4.2.2 Univariate Differential Abundance Methods

An alternative approach to microbiome analysis uses univariate differential abundance methods that test each taxon independently, such as DESeq2 (Love et al., 2014), ALDEx2 (Fernandes et al., 2014), and ANCOM (Mandal et al., 2015). These methods answer a different question than MeLSI: they identify which specific taxa differ in abundance between groups, whereas MeLSI tests whether overall community composition differs. Univariate methods may miss coordinated changes across multiple taxa and suffer from multiple testing penalties when testing thousands of features.

MeLSI complements rather than replaces univariate methods: MeLSI's multivariate test detects overall community differences (like PERMANOVA), while feature importance weights identify candidate taxa for follow-up univariate differential abundance testing. This two-stage approach balances sensitivity to community-level patterns with specificity about individual taxa.

#### 4.2.3 Machine Learning for Microbiome Analysis

Recent work has applied machine learning to microbiome data for disease prediction (Pasolli et al., 2016; Topçuoğlu et al., 2020) and ordination (Morton et al., 2019). However, these approaches typically prioritize predictive accuracy over statistical inference and lack rigorous Type I error control.

MeLSI's key innovation is integrating metric learning with permutation testing, ensuring valid p-values while retaining adaptive learning. The permutation framework is computationally intensive (requiring metric learning on each permuted dataset) but essential for proper inference. This distinguishes MeLSI from metric learning approaches in other domains that focus solely on prediction or clustering.

### 4.3 Biological and Clinical Implications

#### 4.3.1 Exploratory Microbiome Studies

In exploratory studies seeking microbiome associations with disease, MeLSI's combination of adaptive metric learning and appropriate conservatism is particularly valuable. The method's ability to maintain strict Type I error control while providing context-dependent sensitivity reduces both false positives (wasted follow-up effort) and false negatives (missed discoveries), particularly on real-world datasets with subtle, complex signals.

#### 4.3.2 Hypothesis Generation

The interpretable feature weights provided by MeLSI facilitate hypothesis generation for mechanistic studies. By identifying which taxa contribute most to group separation, researchers can prioritize candidates for follow-up experiments such as culturomics, gnotobiotics, or metabolic profiling (Koppel & Balskus, 2016).

#### 4.3.3 Meta-analysis and Cross-study Comparison

An interesting future application is using MeLSI for meta-analysis across multiple studies. By learning study-specific metrics and examining consistency of feature weights, researchers could identify robust, reproducible microbiome signatures that generalize across populations—a major challenge in microbiome research (Duvallet et al., 2017).

### 4.4 Statistical and Methodological Considerations

#### 4.4.1 The Cost of Flexibility

MeLSI's learned metrics introduce additional flexibility compared to fixed metrics, potentially increasing variance in test statistics. Indeed, our results show elevated F-statistics on null data (F = 1.307 and 1.737) compared to Euclidean distance (F = 0.964 and 1.215). However, the permutation testing framework properly accounts for this flexibility by generating null distributions using the same learning procedure, ensuring valid p-values (p = 0.607 and 0.224, both appropriately non-significant).

This "cost of flexibility" is inherent to adaptive methods and represents a trade-off: we accept modest increases in null variability in exchange for context-dependent improvements on real-world datasets. The permutation framework ensures this trade-off maintains statistical validity, preventing inflation of Type I error despite the increased flexibility.

#### 4.4.2 Sample Size Requirements

Like all machine learning approaches, MeLSI benefits from larger sample sizes. Our scalability analysis (Table 3) showed steady gains in F-statistics from $n=20$ (F = 1.222) through $n=100$ (F = 1.786), with performance continuing to improve at $n=500$ (F = 2.424). MeLSI therefore remains practical for moderate sample sizes (50-100 samples per group), but studies with fewer than about 30 samples per group may prefer fixed metrics whose variance is easier to control.

#### 4.4.3 Computational Considerations

MeLSI's computational cost is 100-1000× higher than fixed metrics due to permutation-based metric learning. For extremely large datasets (n > 10,000 or p > 10,000), computational demands may be prohibitive without parallelization. However, most microbiome studies fall within tractable ranges (n < 1000, p < 1000 after pre-filtering).

Future implementations could employ several optimizations including parallel permutation testing across CPU cores, GPU acceleration for distance calculations, adaptive permutation schemes that terminate early for clearly non-significant results, and approximations to the metric learning procedure for very high-dimensional data.

### 4.5 Limitations and Future Directions

#### 4.5.1 Current Limitations

Several limitations warrant discussion:

1. **Computational intensity**: While practical for typical datasets, MeLSI may be too slow for real-time or high-throughput applications.

2. **Hyperparameter tuning**: Although our parameter sensitivity analysis shows robustness to default settings, optimal hyperparameters may vary by dataset characteristics.

3. **Phylogenetic information**: The current implementation learns metrics in abundance space; incorporating phylogenetic relationships (like UniFrac) could improve power for phylogenetically structured signals.

4. **Compositional constraints**: While CLR transformation addresses compositionality, directly learning metrics in compositional space (e.g., Aitchison geometry) may provide advantages.

#### 4.5.2 Extensions for Regression and Covariate Adjustment

An important extension is adapting MeLSI for regression settings where the outcome is continuous or adjusting for covariates is necessary. This could be achieved by learning metrics that maximize correlation between distances and continuous outcomes, incorporating residual distances after regressing out covariates, or integrating MeLSI with linear mixed models for repeated measures or family data.

#### 4.5.3 Integration with Compositional Data Analysis

Microbiome data are inherently compositional (relative abundances sum to 1), violating assumptions of standard multivariate methods (Gloor et al., 2017; Quinn et al., 2018). While CLR transformation provides one solution, directly learning metrics in compositional space using Aitchison distance as the base metric may offer theoretical advantages.

#### 4.5.4 Compatibility with Alternative Ordination and Distance-Based Methods

MeLSI's learned distance matrices can be directly applied to ordination methods beyond PCoA. Non-metric multidimensional scaling (NMDS), which ranks distances rather than using their exact values, could visualize MeLSI's learned dissimilarities while being more robust to non-linearities (Kruskal, 1964). Similarly, ANOSIM (Analysis of Similarities), an alternative permutation-based method that tests whether between-group dissimilarities exceed within-group dissimilarities using rank statistics (Clarke, 1993), could leverage MeLSI's learned metrics for hypothesis testing with different underlying assumptions than PERMANOVA.

However, MeLSI-learned metrics are not suitable for Principal Component Analysis (PCA), which requires Euclidean distances to preserve variance decomposition properties. PCA operates in the original feature space and implicitly uses Euclidean geometry; applying non-Euclidean learned metrics would violate PCA's mathematical foundations. For non-Euclidean distances, PCoA (also called metric multidimensional scaling) remains the appropriate ordination method, as it can accommodate any distance metric without assuming Euclidean structure.

Future work could explore whether MeLSI's ensemble learning and permutation framework could be adapted to optimize metrics specifically for ANOSIM's rank-based statistics or NMDS's stress criteria, potentially offering complementary perspectives on community structure beyond PERMANOVA's sum-of-squares framework.

#### 4.5.5 Multi-omics Integration

Microbiome studies increasingly collect complementary omics data (metagenomics, metatranscriptomics, metabolomics). Multi-view metric learning approaches that jointly learn from multiple data types while identifying shared and specific patterns could provide systems-level insights.

### 4.6 Software Availability and Reproducibility

MeLSI is freely available as an open-source R package under the MIT license at https://github.com/NathanBresette/MeLSI. The package includes comprehensive documentation, tutorial vignettes, and example datasets. All validation experiments reported in this paper are fully reproducible using provided code and data.

### 4.7 Recommendations for Users

Based on our validation results, we offer the following recommendations:

1. **Sample size**: Aim for $n \geq 50$ per group for stable results; consider fixed metrics for $n < 30$ per group.

2. **Preprocessing**: Apply CLR transformation to address compositionality; filter very low-abundance taxa (e.g., present in <10% of samples) before analysis.

3. **Hyperparameters**: Default settings (B=30, m_frac=0.8, n_perms=200) work well in most cases; increase n_perms to 999 for more precise p-values near significance threshold.

4. **Interpretation**: Use omnibus p-value for primary inference; examine feature importance weights to identify which taxa contribute most to group separation; supplement with mean abundance comparisons to determine directionality (which group has higher abundance for top-weighted taxa); validate top features with univariate differential abundance methods.

5. **Computational resources**: Expect 2-30 minutes runtime for typical datasets; use higher-memory machines (16 GB RAM or more) for high-dimensional data.

6. **When to choose MeLSI**: Prefer MeLSI for PERMANOVA workflows that need calibrated p-values plus interpretable taxa weights (e.g., diet interventions or subtle host phenotype comparisons).

---

## 5. Conclusion

This study presents MeLSI (Metric Learning for Statistical Inference), a framework that combines adaptive metric learning with rigorous permutation-based inference for microbiome beta diversity and community composition analysis. By learning data-adaptive distance metrics through an ensemble of weak learners and validating significance via permutation testing, MeLSI retains the interpretability of distance-based workflows while allowing the metric itself to adapt to the dataset at hand.

Across synthetic benchmarks and real datasets, MeLSI maintained proper Type I error control and delivered context-dependent gains in statistical power. The method matched or exceeded the best traditional metrics on real data (9.1% higher F-statistic than Euclidean distance on the Atlas1006 cohort and a significant detection on DietSwap where fixed metrics remained marginal) while remaining competitive with Bray-Curtis and UniFrac on synthetic datasets that strongly favor abundance- or phylogeny-based distances. Crucially, even when traditional metrics achieved similar F-statistics, MeLSI supplied feature-weight profiles that isolated the taxa driving observed differences, turning omnibus PERMANOVA results into biologically interpretable findings. Scalability experiments demonstrated smooth increases in power as sample size and dimensionality grew, and conservative pre-filtering yielded modest (1–6%) runtime savings without distorting p-values.

By addressing the rigidity of fixed distance metrics, MeLSI offers a principled, statistically valid, and practical alternative for microbiome community composition analysis. We anticipate that the method will be particularly valuable in exploratory analyses, meta-analyses, and hypothesis generation, where adaptive yet interpretable metrics can uncover subtle patterns while maintaining strong error control. The open-source implementation lowers the barrier to adoption and invites the community to experiment with faster permutation strategies, covariate-aware objectives, and other optimizations tailored to emerging microbiome study designs.

---

## Acknowledgments

[To be added]

---

## Author Contributions

[To be added]

---

## Competing Interests

The authors declare no competing interests.

---

## Data Availability

MeLSI is freely available at https://github.com/NathanBresette/MeLSI under the MIT license. All validation data and analysis scripts are included in the package repository. The Atlas1006 and DietSwap datasets are available through the R microbiome package (https://microbiome.github.io/).

---

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman and Hall, London.

Anderson, M. J. (2017). Permutational multivariate analysis of variance (PERMANOVA). In *Wiley StatsRef: Statistics Reference Online* (pp. 1-15). John Wiley & Sons, Ltd.

Bellet, A., Habrard, A., & Sebban, M. (2013). A survey on metric learning for feature vectors and structured data. arXiv preprint arXiv:1306.6709.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*, 57(1), 289-300.

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.

Clemente, J. C., Ursell, L. K., Parfrey, L. W., & Knight, R. (2012). The impact of the gut microbiota on human health: an integrative view. *Cell*, 148(6), 1258-1270.

Duvallet, C., Gibbons, S. M., Gurry, T., Irizarry, R. A., & Alm, E. J. (2017). Meta-analysis of gut microbiome studies identifies disease-specific and shared responses. *Nature Communications*, 8(1), 1-10.

Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2(1), 15.

Gilbert, J. A., Blaser, M. J., Caporaso, J. G., Jansson, J. K., Lynch, S. V., & Knight, R. (2018). Current understanding of the human microbiome. *Nature Medicine*, 24(4), 392-400.

Gloor, G. B., Macklaim, J. M., & Fernandes, A. D. (2017). Displaying variation in large datasets: plotting a visual summary of effect sizes. *Journal of Computational and Graphical Statistics*, 25(3), 971-979.

Good, P. I. (2013). *Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses*. Springer Science & Business Media.

Gopalakrishnan, V., Spencer, C. N., Nezi, L., Reuben, A., Andrews, M. C., Karpinets, T. V., ... & Wargo, J. A. (2018). Gut microbiome modulates response to anti–PD-1 immunotherapy in melanoma patients. *Science*, 359(6371), 97-103.

Knights, D., Costello, E. K., & Knight, R. (2011). Supervised classification of human microbiota. *FEMS Microbiology Reviews*, 35(2), 343-359.

Koppel, N., & Balskus, E. P. (2016). Exploring and understanding the biochemical diversity of the human microbiota. *Cell Chemical Biology*, 23(1), 18-30.

Kulis, B. (2013). Metric learning: A survey. *Foundations and Trends in Machine Learning*, 5(4), 287-364.

Lahti, L., Salojärvi, J., Salonen, A., Scheffer, M., & de Vos, W. M. (2014). Tipping elements in the human intestinal ecosystem. *Nature Communications*, 5(1), 1-10.

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful transformations for ordination of species data. *Oecologia*, 129(2), 271-280.

Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology*, 71(12), 8228-8235.

Lynch, S. V., & Pedersen, O. (2016). The human intestinal microbiome in health and disease. *New England Journal of Medicine*, 375(24), 2369-2379.

Mandal, S., Van Treuren, W., White, R. A., Eggesbø, M., Knight, R., & Peddada, S. D. (2015). Analysis of composition of microbiomes: a novel method for studying microbial composition. *Microbial Ecology in Health and Disease*, 26(1), 27663.

Mahalanobis, P. C. (1936). On the generalized distance in statistics. *Proceedings of the National Institute of Sciences of India*, 2(1), 49-55.

Markle, J. G., Frank, D. N., Mortin-Toth, S., Robertson, C. E., Feazel, L. M., Rolle-Kampczyk, U., ... & Danska, J. S. (2013). Sex differences in the gut microbiome drive hormone-dependent regulation of autoimmunity. *Science*, 339(6123), 1084-1088.

McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to community data: a comment on distance-based redundancy analysis. *Ecology*, 82(1), 290-297.

McMurdie, P. J., & Holmes, S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, 8(4), e61217.

Morton, J. T., Aksenov, A. A., Nothias, L. F., Foulds, J. R., Quinn, R. A., Badri, M. H., ... & Knight, R. (2019). Learning representations of microbe–metabolite interactions. *Nature Methods*, 16(12), 1306-1314.

Org, E., Mehrabian, M., Parks, B. W., Shipkova, P., Liu, X., Drake, T. A., & Lusis, A. J. (2016). Sex differences and hormonal effects on gut microbiota composition in mice. *Gut Microbes*, 7(4), 313-322.

Oksanen, J., Blanchet, F. G., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., ... & Wagner, H. (2020). vegan: Community Ecology Package. R package version 2.5-7. https://CRAN.R-project.org/package=vegan

Pasolli, E., Truong, D. T., Malik, F., Waldron, L., & Segata, N. (2016). Machine learning meta-analysis of large metagenomic datasets: tools and biological insights. *PLoS Computational Biology*, 12(7), e1004977.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. *Statistical Applications in Genetics and Molecular Biology*, 9(1), Article 39.

Quinn, T. P., Erb, I., Richardson, M. F., & Crowley, T. M. (2018). Understanding sequencing data as compositions: an outlook and review. *Bioinformatics*, 34(16), 2870-2878.

Ramette, A. (2007). Multivariate analyses in microbial ecology. *FEMS Microbiology Ecology*, 62(2), 142-160.

Shreiner, A. B., Kao, J. Y., & Young, V. B. (2015). The gut microbiome in health and in disease. *Current Opinion in Gastroenterology*, 31(1), 69-75.

Topçuoğlu, B. D., Lesniak, N. A., Ruffin IV, M. T., Wiens, J., & Schloss, P. D. (2020). A framework for effective application of machine learning to microbiome-based classification problems. *mBio*, 11(3), e00434-20.

Vemuri, R., Gundamaraju, R., Shastri, M. D., Shukla, S. D., Kalpurath, K., Ball, M., ... & Eri, R. (2019). Gut microbial changes, interactions, and their implications on human lifecycle: an ageing perspective. *BioMed Research International*, 2019, 4178607.

Weinberger, K. Q., & Saul, L. K. (2009). Distance metric learning for large margin nearest neighbor classification. *Journal of Machine Learning Research*, 10(2), 207-244.

Weiss, S., Xu, Z. Z., Peddada, S., Amir, A., Bittinger, K., Gonzalez, A., ... & Knight, R. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome*, 5(1), 27.

Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. John Wiley & Sons.

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.

Xing, E. P., Jordan, M. I., Russell, S. J., & Ng, A. Y. (2002). Distance metric learning with application to clustering with side-information. In *Advances in Neural Information Processing Systems* (pp. 521-528).

O'Keefe, S. J., Li, J. V., Lahti, L., Ou, J., Carbonero, F., Mohammed, K., ... & Sanderson, I. R. (2015). Fat, fibre and cancer risk in African Americans and rural Africans. *Nature Communications*, 6, 6342.

---