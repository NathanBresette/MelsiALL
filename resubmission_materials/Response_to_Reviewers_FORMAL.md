# Response to Reviewers
## MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis
**Manuscript ID:** mSystems01710-25

**Authors:** Nathan Bresette, Aaron C. Ericsson, Carter Woods, Ai-Ling Lin

---

## Cover Note

We thank the editor and reviewers for their constructive feedback. We have addressed all concerns, with particular attention to statistical validation. Below we provide point-by-point responses with references to the revised manuscript.

### Summary of Major Revisions

- **Type I error control:** Expanded from a single simulation to 600 total simulations (100 per condition × 3 sample sizes × 2 dataset types), confirming 3-6% rejection rates at α = 0.05
- **Statistical power:** Expanded to 450 simulations across three effect sizes and three sample sizes, with individual method comparisons (Supplementary Table S2)
- **Feature correlation robustness:** Added 200 simulations across four correlation levels (r = 0, 0.3, 0.6, 0.8), demonstrating stable performance (Table 5)
- **Pre-filtering analysis:** Expanded to 150 simulations with empirical power evaluation (Table 6)
- **Parameter sensitivity:** Expanded to 275 simulations with single-learner baseline comparison (Table 4, Supplementary Table S4)
- **Signal taxa recovery:** Added Supplementary Table S1 with Precision, Recall, Mean Rank, and AUC-ROC metrics
- **SKIOME dataset:** Added third real-world validation dataset (skin microbiome, 511 samples, 3 groups) demonstrating multi-group and cross-body-site applicability (Figure 3)
- **DietSwap figures:** Added VIP and PCoA plots (Figure 2)
- **Atlas1006 PCoA:** Removed standalone PCoA figure to avoid overstating visual separation; retained VIP plot (Figure 1)
- **Confidence ellipses:** Updated all PCoA plots from 68% to 95%
- **CLR discussion:** Added explicit discussion of CLR trade-offs and when count-based methods are preferable
- **Notation:** Systematically introduced mathematical notation in the Introduction
- **Table annotations:** Added standardized footnotes with abbreviation definitions to all tables

---

## RESPONSE TO REVIEWER #1

**Reviewer Comment:** This manuscript proposes MeLSI, a supervised, data-driven metric learning approach for constructing distance matrices in microbiome studies, followed by PERMANOVA-based inference for beta-diversity analysis. The manuscript emphasizes (i) rigorous Type I error control, (ii) statistical power, (iii) robustness and overfitting prevention via ensemble learning, and (iv) enhanced biological interpretability. However, while the methodological idea is promising, several key claims are not fully supported by the current experimental design or presented evidence.

**Response:** We appreciate this recognition and have addressed each concern below.

---

### Feature Pre-filtering

**Reviewer Comment:** The pre-filtering strategy relies on t-statistic-like measures for two-group comparisons and ANOVA F-statistics for multi-group settings. Given the count-based, compositional nature of microbiome data, it would be helpful for the authors to discuss the implications of using statistics that assume Gaussian distributions.

**Response:** We clarify that these statistics serve as ranking heuristics, not inferential tests. We use them solely to rank features by between-group signal, without relying on distributional assumptions. All statistical inference uses permutation testing, which is distribution-free.

**Location:** Lines 155-160

---

**Reviewer Comment:** The manuscript does not illustrate multi-group extensions using synthetic data, which would help assess generalizability and statistical properties beyond the two-group setting.

**Response:** We address multi-group generalizability through real-world validation on the SKIOME skin microbiome dataset (511 samples, 3 groups: Atopic Dermatitis, Healthy, Psoriasis). MeLSI detected significant differences (F = 4.895, p = 0.005), with all pairwise comparisons significant after FDR correction. The permutation framework ensures valid inference regardless of group number. Comprehensive multi-group synthetic validation (>1,500 additional simulations) represents future work.

**Location:** Limitations (lines 565-570), Methods (lines 304-311), SKIOME results (lines 501-516)

---

### Double Dipping / Overfitting Concerns

**Reviewer Comment:** The distance metric is learned in a supervised manner and subsequently used for hypothesis testing on the same dataset. This raises a concern regarding "double dipping." While bagging is employed to mitigate overfitting and the authors claim the proposed method prevents overfitting across the manuscript, the potential impact of in-sample estimation and inference remains unclear and should be discussed more explicitly or addressed through additional validation strategies.

**Response:** The metric is relearned on each permutation, preventing double-dipping. Each permutation represents an independent metric learning experiment under the null hypothesis—the metric is never optimized on observed data and then tested on the same data. The p-value properly accounts for the adaptive nature of the method.

**Location:** Lines 218-229

---

**Reviewer Comment:** Statements regarding overfitting prevention across the manuscript would benefit from explicit evaluation or supporting evidence.

**Response:** Two lines of evidence confirm overfitting prevention: (1) Table 4 shows ensemble learning reduces variance 4× compared to a single learner (SD = 0.119-0.128 vs. 0.505), and (2) proper Type I error control across 100 simulations (3-6% rejection rates, Table 1) confirms overfitting does not inflate false positives.

**Location:** Table 4 (line 391), Table 1 (line 343), Parameter Sensitivity (line 387)

---

**Reviewer Comment:** The ensemble learning strategy is motivated by robustness and overfitting prevention, but no direct comparison is provided with existing metric learning approach based on a single-learner.

**Response:** Table 4 now includes a single-learner baseline (B=1), which shows substantially higher variance (SD = 0.505) and higher p-values (mean = 0.421) compared to ensemble approaches (SD = 0.119-0.128 for B≥10). Bootstrap sampling and feature subsampling create diversity among learners, following established ensemble learning practice (20).

**Location:** Table 4 (line 391), Parameter Sensitivity (line 387)

---

### Major Concerns on Experimental Design

#### Type I Error Control

**Reviewer Comment:** Type I error control appears to be evaluated using a single synthetic dataset, with conclusions drawn from a single p-value. Repeated simulations are necessary to support claims of proper error control. Type I error control is a distributional property that must be assessed over repeated realizations of the null hypothesis. Demonstrating that a single test yields p-values greater than 0.05 in one or two null examples does not quantify the probability of false positive findings. Proper evaluation would require repeated simulations under the null, with estimation of the empirical rejection rate at the chosen significance level and, ideally, examination of the null p-value distribution.

**Response:** We have expanded the Type I error analysis to include 100 simulations per condition across three sample sizes (n=50, 100, 200) for both synthetic null data and real shuffled data (600 total simulations: 2 dataset types × 3 sample sizes × 100 simulations). Table 1 now reports empirical rejection rates at α = 0.05 for each sample size and dataset type, confirming Type I error rates near the nominal 5% level (range: 3-6% across all conditions). The Conclusions section has been updated accordingly.

**Location:** Table 1 (line 343), Results (lines 346-350), Conclusions (lines 534-536)

---

#### Statistical Power Analysis

**Reviewer Comment:** Similar concerns apply to the statistical power analysis, which also appears to rely on a single dataset per setting.

**Response:** We have expanded the power analysis to include 50 simulations per condition across three effect sizes (small, medium, large) and three sample sizes (n=50, 100, 200), totaling 450 simulations. Table 2 now reports empirical power (detection rates) and mean F-statistics for each combination. Individual method comparisons in Supplementary Table S2 show MeLSI consistently outperforms Jaccard and Unweighted UniFrac while demonstrating appropriate conservatism for small effects.

**Location:** Table 2 (line 354), method comparisons (lines 360-369), power analysis (lines 351-369)

---

#### Sample Size Exploration

**Reviewer Comment:** Exploring multiple sample sizes in both Type I error and power analyses would strengthen the evaluation. (may consider a few typical sample sizes in microbiome studies)

**Response:** Both Type I error and power analyses now include multiple sample sizes (n=50, 100, 200), covering typical microbiome study sizes. The revised analyses report empirical Type I error rates and statistical power at each sample size, complementing the computational scalability analysis in Table 3. Results confirm proper Type I error control and appropriate power gains with increasing sample size.

**Location:** Table 1 (line 343), Table 2 (line 354), Table 3 (line 375), Results (lines 337-350, 351-369, 370-386)

---

#### Feature Correlation in Simulations

**Reviewer Comment:** Synthetic data generation assumed independent taxa; given the multivariate and correlated structure of microbiome data. Simulations incorporating feature correlation would be important.

**Response:** We have added validation under varying feature correlation levels (Table 5). The analysis evaluated four correlation levels (r=0, 0.3, 0.6, 0.8) using 50 simulations per condition (200 total). MeLSI is robust to correlated features: F-statistics remained stable (±1.7% variation), power was consistent (42-50%), and feature recovery metrics showed minimal variation. MeLSI achieved its best ranking (1/6) at high correlation (r=0.8), suggesting particular suitability for correlated microbiome data. The ensemble approach naturally handles correlated features by aggregating signal across correlated taxa.

**Location:** Table 5 (line 409), Feature Correlation Robustness (lines 401-428)

---

### Interpretability Advantage

**Reviewer Comment:** Emphasizing the interpretability advantage of the proposed method, particularly in identifying taxa that drive group differences, would improve positioning of the manuscript. It would strengthen the manuscript to examine how learned weights recover true signal taxa, possibly across varying effect sizes, sample sizes, and feature dimensions.

**Response:** We have added validation of signal taxa recovery across varying effect sizes and sample sizes. Supplementary Table S1 provides detailed recovery metrics: Precision at k, Recall at k, Mean Rank, and AUC-ROC. MeLSI effectively recovers true signal taxa, with performance improving with effect size and sample size. For large effects, Precision at 5 reached 0.876-1.000 and Mean Rank decreased to 14.4, confirming that true signal taxa are consistently ranked among the top features.

**Location:** Results section (lines 355-359, Table 2 footnote referencing Supplementary Table S1); Supplementary Table S1

---

### Real Data Analyses

**Reviewer Comment:** The observed outperformance of MeLSI on real datasets may partly reflect overfitting or double-dipping, despite the ensemble strategy. This possibility should be discussed more explicitly.

**Response:** MeLSI's outperformance on real datasets does not reflect overfitting because the permutation framework relearns the metric on each permutation. This is directly confirmed by proper Type I error control on real shuffled Atlas1006 data (3-6% rejection rates across 100 simulations, Table 1), demonstrating genuine signal detection rather than overfitting.

**Location:** Table 1 (line 343, "Null Real Shuffled" rows), Null distribution generation (lines 218-229)

---

**Reviewer Comment:** In Figure 2, the manuscript describes a "clear separation" between male and female samples along PCoA1; however, substantial overlap between groups appears to remain. Can authors elaborate this? Is this separation clearer compared to the results of traditional approaches?

**Response:** We acknowledge this concern. We have removed the standalone Atlas1006 PCoA figure (previously Figure 2) to avoid overstating visual separation. The statistical significance (F=5.141, p=0.005) is reported in the text, and the VIP plot (Figure 1) provides the primary interpretability value. We now show PCoA ordination for DietSwap and SKIOME (Figures 2-3), where group separation is more visually apparent, confirming that MeLSI-learned distances are compatible with standard ordination approaches.

**Location:** Figure 1 (lines 474-479), Figures 2-3 (lines 497-500, 511-516)

---

**Reviewer Comment:** Could authors explain the rationale for presenting 68% confidence ellipses, rather than more conventional choice of 95%?

**Response:** All figures now use the conventional 95% confidence ellipses. The manuscript text and figure captions reflect this change, and all figures have been regenerated accordingly. The code has also been updated to generate 95% confidence ellipses for all PCoA plots.

**Location:** Figures 2-3 (DietSwap and SKIOME PCoA plots); Code: `github/R/melsi_robust.R` line 1227

---

**Reviewer Comment:** The significant PERMANOVA result (F=5.141, p=0.005) may be influenced by the large sample size (N>1000). I was wondering if the seemingly subtle separation shown in the figure is considered biologically meaningful besides statistical significance.

**Response:** We acknowledge that the large sample size (n=1,114) contributes to statistical significance. The revised text (lines 458-461) states that "MeLSI's improvement over the best fixed metric suggests that learned metrics can capture biologically relevant patterns in subtle, high-dimensional comparisons, consistent with previously documented sex-associated microbiome differences (29, 30)." The learned feature weights (Figure 1) provide actionable biological insight regardless of sample size.

**Location:** Lines 458-461, Figure 1 (lines 474-479)

---

**Reviewer Comment:** There is an inconsistency in Figure 2, where PCoA1 variance is listed as 18.4% in the legend but 21.5% on the x-axis.

**Response:** This issue has been resolved by removing the standalone Atlas1006 PCoA figure (previously Figure 2) in response to the reviewer's concern about visual separation. All remaining PCoA figures (Figures 2-3) have been verified to ensure consistency between x-axis labels (dynamically generated from PCoA calculations) and figure captions.

**Location:** Figures 2-3 (DietSwap and SKIOME PCoA plots)

---

**Reviewer Comment:** Transition to the DietSwap data analysis is abrupt, and corresponding results including VIP feature and PCoA plots are not shown.

**Response:** We have improved the transition by adding introductory text (lines 451-454): "To evaluate MeLSI's utility in real-world applications, we analyzed three published microbiome datasets: Atlas1006 (sex-associated differences), DietSwap (dietary intervention), and SKIOME (multi-group skin microbiome validation)." VIP and PCoA plots for DietSwap have been added (Figure 2).

**Location:** Lines 451-454, Figure 2 (lines 497-500)

---

### Ensemble Size Analysis

**Reviewer Comment:** In the ensemble size analysis, the authors refer to "variance" in performance, but Table 4 does not report variance estimates. If results are based on a single synthetic dataset, clarification of what is meant by variance is needed.

**Response:** We have expanded the parameter sensitivity analysis to include 25 replications per parameter value (275 total simulations), enabling robust variance estimation. Table 4 now reports mean F-statistics, p-values, and computation times, with a single-learner baseline (B=1) showing substantially higher variance (SD = 0.505) compared to ensemble approaches (SD = 0.119-0.128). Standard deviations are provided in Supplementary Table S4. F-statistics remained stable across ensemble sizes (B=10-100) compared to the single-learner baseline.

**Location:** Table 4 (line 391), Parameter Sensitivity (line 387)

---

### Pre-filtering Analysis

**Reviewer Comment:** The sample sizes used in the synthetic pre-filtering analyses are not clearly stated.

**Response:** We have clarified that all pre-filtering analyses use n=100 samples per condition. The text describing the experimental conditions for Table 6 now explicitly states "n=100 samples per condition" along with the dimensionalities for each effect size.

**Location:** Pre-filtering analysis section (lines 429-434)

---

**Reviewer Comment:** Feature dimensionality differs across effect size settings (p=500 for small, 200 for medium, and 100 for large). Could authors provide the details on this design choice to help interpretation?

**Response:** We have clarified in the pre-filtering analysis section that the varying dimensionality (p=500 for small, p=200 for medium, p=100 for large) is now explicitly stated alongside sample sizes. This design tests pre-filtering benefits across different dimensionalities, with higher-dimensional datasets (p=500) showing greater time savings (39.8%) compared to lower-dimensional datasets (p=100, 16.5% time savings).

**Location:** Pre-filtering analysis section (lines 429-434)

---

**Reviewer Comment:** As with other simulation sections, results based on a single synthetic dataset cannot be interpreted as statistical power.

**Response:** We have expanded the pre-filtering analysis to include 50 simulations per scenario (150 total). Table 6 now reports empirical power, mean F-statistics, and computational time savings.

Pre-filtering is most impactful in high-dimensional, sparse-signal settings. In the high-dimensional setting (p=500, small effect: 1.5× fold change in 5 signal taxa), pre-filtering removed 150 noise features, increasing power from 4% to 100%. In the lower-dimensional setting (p=100, large effect: 3.0× fold change in 20 signal taxa), power increased from 14% to 84%. The higher power in the small-effect scenario reflects the experimental design: pre-filtering provides maximum benefit when the signal-to-noise ratio is lowest (higher dimensionality), as removing 30% of non-informative features in the p=500 environment recovers signal masked by noise.

**Location:** Table 6 (line 435), Pre-filtering analysis (lines 429-450)

---

### Directionality and Log2 Fold-Change

**Reviewer Comment:** The manuscript states that MeLSI provides directionality and log2 fold-change information for each taxon. Additional details on how this information is derived would be helpful.

**Response:** We have clarified the derivation of directionality and effect size in the Results section. **Directionality:** Is determined by identifying which group has the higher mean abundance on the CLR-transformed data, ensuring consistency with the metric learning process. **Effect Size:** We report the log2 fold change computed from CLR-transformed group means: $\log_2(\mu_{\text{CLR,1}} / \mu_{\text{CLR,2}})$. This provides a standardized measure of the magnitude of difference between groups for each taxon.

**Location:** Lines 480-486

---

### Computational Time Justification

**Reviewer Comment:** The claim that increased computational time is justified by improved statistical power is not fully supported, especially given inadequate validation study designs and given that Table 2 shows performance comparable to the best traditional methods.

**Response:** The revised Limitations section (lines 559-564) justifies MeLSI's computational time (2-30 minutes for typical datasets) based on: (1) interpretability gains through learned feature weights, (2) pre-filtering providing 16-40% time savings while improving power by 36-37% (Table 6), and (3) the modest time investment relative to overall study timelines. For large-scale screening studies, traditional methods may be more appropriate. The justification emphasizes interpretability as the primary benefit, consistent with Table 2 showing comparable power to traditional methods.

**Location:** Limitations (lines 559-564), Table 6 (line 435)

---

### Prevalence Threshold

**Reviewer Comment:** The mention of a 10% prevalence threshold appears only in the computational efficiency section. Clarification is needed on whether prevalence filtering is an expected preprocessing step and how it plays a role along with the pre-filtering process.

**Response:** We have added clarification in the Methods section (Real data sources) explaining that "Prevalence filtering (retaining features present in ≥10% of samples) is an optional preprocessing step distinct from MeLSI's variance-based pre-filtering; when applied, prevalence filtering removes rare taxa before analysis, while MeLSI's pre-filtering focuses on variance-based feature selection after preprocessing." This clarifies the distinction between the two filtering steps and their roles.

**Location:** Methods section, Real data sources (lines 281-284)

---

### Tables

**Reviewer Comment:** Tables would benefit from clearer annotations and definitions of abbreviations.

**Response:** All tables (Tables 1-6) now include concise footnotes defining abbreviations: n (sample size), p (number of taxa/features), F (PERMANOVA F-statistic), Power (empirical statistical power), Time (computation time in seconds), Rank (MeLSI's rank among 6 methods), and other table-specific terms. Footnotes are standardized for consistency across tables.

**Location:** Table 1: Line 344; Table 2: Line 355; Table 3: Line 376; Table 4: Line 392; Table 5: Line 410; Table 6: Line 436

---

### Notation

**Reviewer Comment:** Introduce notation more systematically in the "Metric learning: an emerging paradigm" section.

**Response:** The revised text (lines 88-96) now formally defines: (1) the feature abundance matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$ with $n$ samples and $p$ taxa, (2) group labels $\mathbf{y} = (y_1, \ldots, y_n)$, (3) the positive semi-definite metric matrix $\mathbf{M} \in \mathbb{R}^{p \times p}$, (4) the Mahalanobis distance formula $d_M(\mathbf{x}_i, \mathbf{x}_j) = \sqrt{(\mathbf{x}_i - \mathbf{x}_j)^T \mathbf{M} (\mathbf{x}_i - \mathbf{x}_j)}$, and (5) the special case of diagonal $\mathbf{M}$ reducing to weighted Euclidean distance with feature-specific weights $M_{jj}$.

**Location:** Introduction, "Metric learning: an emerging paradigm" subsection (lines 81-96)

---

### CLR Transformation

**Reviewer Comment:** The role and impact of the CLR transformation in conjunction with the proposed method should be discussed more explicitly. CLR transformation for MeLSI is mentioned only in the software availability section as part of the recommended usage.

**Response:** We have added explicit discussion of the CLR transformation in the Methods section, explaining that: (1) CLR converts relative abundances to log-ratios suitable for Euclidean distance, (2) CLR treats abundance ratios more equitably than count-based metrics, which can be dominated by highly abundant taxa, (3) CLR may attenuate large fold-change signals, as evidenced by traditional count-based methods achieving higher F-statistics on synthetic data with large effects (3× fold change), and (4) CLR is most appropriate when signals are distributed across multiple taxa and interpretability is prioritized. The Results section now states when CLR-based versus count-based methods are preferable.

**Location:** Methods (lines 269-284), Results (lines 360-369)

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** While we can appreciate the novelty of this work, several critical technical concerns need to be addressed for improvement or clarification.

**Response:** We thank the reviewer for the detailed technical feedback. We address each concern below.

---

### 1. Computational Cost vs. Benefit Trade-off

**Reviewer Comment:** MeLSI looks slower (minutes to hours) than traditional PERMANOVA (seconds). The paper argues this is acceptable given the gain in interpretability, but for large-scale or screening studies, this may be prohibitive. No clear power-time trade-off analysis is provided to justify when MeLSI is worth the extra computation.

**Response:** We have added an explicit power-time trade-off analysis (Computational Performance, lines 517-529). Pre-filtering increases power by 36-37% while reducing computation time by 16-40% (Table 6). For typical studies (n=50-200, p=100-500), MeLSI completes in 2-30 minutes (Table 3). The trade-off is most favorable when: (1) sample sizes are moderate, (2) interpretability is prioritized, and (3) pre-filtering is applied. For very large studies (n>500) or rapid screening, traditional methods may be preferable.

**Location:** Computational Performance (lines 517-529), Tables 2 (line 354), 3 (line 375), and 6 (line 435)

---

### 2. Limited Real-Data Performance Gains

**Reviewer Comment:** On Atlas1006, MeLSI's F-statistic was only 9.1% higher than Euclidean distance (5.141 vs. 4.711)-a modest improvement. On DietSwap, MeLSI reached significance (p=0.015) while Bray-Curtis was marginal (p=0.058), but the effect size difference is small. The authors acknowledge that on synthetic data with large effects, traditional metrics (Bray-Curtis, UniFrac) often outperformed MeLSI.

**Response:** The revised Conclusions (lines 545-557) explicitly position when MeLSI versus traditional methods are preferable. MeLSI is recommended when: (1) effect sizes are moderate (2× fold change), (2) interpretability through feature weights is needed, (3) traditional methods yield marginal results (p-values near 0.05), and (4) signals are distributed across multiple taxa. Traditional methods are preferable for: (1) large, obvious effects (3× fold change), (2) large-scale screening studies where speed is critical, and (3) when only omnibus testing is needed without feature-level interpretation.

**Location:** Conclusions (lines 545-557), Table 2 (line 354)

---

### 3. CLR Transformation May Limit Sensitivity

**Reviewer Comment:** MeLSI uses CLR-transformed data for Euclidean distance, which may attenuate large fold-change signals compared to count-based metrics (Bray-Curtis, UniFrac). This could explain why MeLSI underperforms on synthetic data with strong effects.

**Response:** We have added discussion of CLR trade-offs in both the Methods (lines 269-284) and Results (lines 364-369) sections. CLR may attenuate large fold-change signals, as evidenced by traditional count-based methods achieving higher F-statistics with large effects (3× fold change). The CLR-based approach is most appropriate when: (1) signals are distributed across multiple taxa, (2) interpretability is prioritized, and (3) effect sizes are moderate. For large, obvious effects, count-based methods (Bray-Curtis, UniFrac) may be preferable.

**Location:** Methods (lines 269-284), Results (lines 364-369)

---

### 4. Lack of Covariate Adjustment

**Reviewer Comment:** MeLSI currently only supports simple group comparisons. No ability to adjust for confounders (age, BMI, antibiotics) or model continuous outcomes. This limits its utility in observational studies where confounding is common.

**Response:** We acknowledge this limitation (Limitations section, lines 571-577). The current implementation focuses on simple group comparisons, which represent the primary use case for microbiome beta diversity analysis.

Covariate adjustment is a high-priority future extension. The underlying PERMANOVA framework (via vegan's `adonis2`) already supports covariates through formula syntax (e.g., `adonis2(dist ~ group + age + BMI)`), providing a natural pathway. The technical approach would involve learning the metric on residuals after regressing out covariates, or incorporating covariates directly into the objective function. This extension would require separate validation specific to that design type.

**Location:** Limitations (lines 571-577)

---

### 5. Sparse and Compositional Data Challenges

**Reviewer Comment:** While pre-filtering is used, the method does not explicitly model compositionality or zero-inflation beyond CLR. Aitchison geometry or Dirichlet-multinomial models might be more appropriate for compositional data.

**Response:** MeLSI addresses compositionality through CLR transformation. The Aitchison distance is defined as the Euclidean distance in CLR space: $d_A(x,y) = \sqrt{\sum(\text{clr}(x) - \text{clr}(y))^2}$. Our method computes Mahalanobis distance on CLR-transformed data, which is a generalized Aitchison distance: $d_M(x,y) = \sqrt{(\text{clr}(x) - \text{clr}(y))^T \mathbf{M}^{-1} (\text{clr}(x) - \text{clr}(y))}$. When $\mathbf{M} = \mathbf{I}$, this reduces exactly to Aitchison distance; when $\mathbf{M} \neq \mathbf{I}$ (learned from data), this is a weighted Aitchison distance that adaptively weights dimensions based on their contribution to group separation. MeLSI thus uses Aitchison geometry as the foundation but extends it with learned feature-specific weights. Zero-inflation is handled through pseudocounts (adding 1 before log transformation). This approach maintains statistical validity through permutation testing, as demonstrated by proper Type I error control (Table 1). Explicit zero-inflation models represent a potential future enhancement.

**Location:** Methods (lines 269-284), Limitations (lines 575-577)

---

### 6. Validation on Limited Real Datasets

**Reviewer Comment:** Only two real datasets (Atlas1006, DietSwap) were used for benchmarking, both from gut microbiome studies. Performance in other body sites (oral, skin) or disease cohorts is unknown.

**Response:** We have expanded validation to three published datasets: (1) **Atlas1006** (gut microbiome, 1,114 samples, sex-associated differences), (2) **DietSwap** (gut microbiome, 74 baseline samples, dietary intervention), and (3) **SKIOME** (skin microbiome, 511 samples, three-group comparison: Atopic Dermatitis, Healthy, Psoriasis).

On SKIOME, MeLSI detected significant differences (F = 4.895, p = 0.005), comparable to Euclidean distance (F = 4.897) but lower than count-based methods (Bray-Curtis: F = 16.275, Jaccard: F = 11.058). All pairwise comparisons remained significant after FDR correction. While count-based methods achieved higher F-statistics, MeLSI provides unique interpretability through learned feature weights identifying which taxa drive group separation. Validation in additional body sites (oral, vaginal) would further strengthen generalizability, though the current three datasets span different body sites (gut, skin), study designs, and sample sizes (n=74 to n=1,114).

**Location:** Real data validation (lines 451-516), SKIOME subsection (lines 501-516), Figure 3 (lines 511-516)

---

### 7. Parameter Sensitivity and Defaults

**Reviewer Comment:** Sensitivity analysis shows robustness, but defaults (e.g., B=30, m_frac=0.8) are not rigorously justified. For some datasets, different hyperparameters might yield better performance.

**Response:** The revised Parameter Sensitivity section explicitly justifies defaults based on Table 4: (1) B=30 provides F-statistics (mean = 1.530, SD = 0.123) comparable to larger ensembles (B=50-100) with reasonable computation time (576.8s), (2) m_frac=0.8 balances performance with learner diversity (lower values show slightly higher F but reduced diversity; higher values show slightly lower F), and (3) robustness across wide parameter ranges (B=10-100, m_frac=0.5-1.0) indicates defaults provide good performance, though users may optimize for specific datasets.

**Location:** Parameter Sensitivity (lines 387-400), Table 4 (line 391)

---

### 8. Comparison to Other ML Methods

**Reviewer Comment:** No comparison to other interpretable ML models (e.g., Random Forest feature importance, logistic regression with regularization) that also provide taxa rankings. MeLSI is presented as unique, but similar insights might be obtained from simpler models. If authors do so, the ground truth should be well defined for benchmarking.

**Response:** MeLSI addresses a fundamentally different question than RF or penalized regression: inference on community-level differences rather than sample classification. RF identifies features maximizing classification accuracy, whereas MeLSI learns a distance metric optimizing group separation in beta-diversity space — enabling PERMANOVA-based hypothesis testing with interpretable feature weights. Additionally, standard ML importance scores lack a formal framework for $p$-value generation in community composition testing; MeLSI provides this through its permutation-based inference framework.

The Conclusions now clarify that "unlike prediction-focused machine learning (e.g., Random Forest, neural networks), MeLSI is an inference-focused approach: every learned metric undergoes rigorous permutation testing to ensure that p-values remain valid despite the adaptive nature of the method." The appropriate comparisons for MeLSI are other beta diversity methods used with PERMANOVA (Bray-Curtis, Euclidean, Jaccard, UniFrac), which we comprehensively evaluate.

**Location:** Conclusions, Summary subsection (lines 552-557)

---

### 9. Scalability to Very High Dimensions

**Reviewer Comment:** The method scales as O(p²) in the number of features. For shotgun metagenomics (thousands of species or genes), this could be computationally prohibitive even with pre-filtering.

**Response:** The revised Scalability section (lines 370-386) acknowledges that O(p²) scaling becomes prohibitive for p>1000, with Table 3 showing computation time increases from 244.8s at p=50 to 8633.0s at p=1000. Pre-filtering (retaining 70% of features) substantially mitigates this. For shotgun metagenomics, we recommend: (1) applying pre-filtering, (2) feature aggregation (e.g., species-level rather than gene-level), or (3) traditional methods if interpretability is not prioritized. The current implementation is most suitable for typical 16S rRNA datasets (p<1000) and metagenomic datasets with moderate dimensionality after preprocessing.

**Location:** Scalability (lines 370-386), Table 3 (line 375)

---

### 10. No Longitudinal/Paired Design Support

**Reviewer Comment:** MeLSI does not currently support repeated measures or paired samples, which are common in microbiome intervention studies.

**Response:** PERMANOVA (via vegan's `adonis2`) already supports paired and longitudinal designs through the `strata` argument, which restricts permutations within blocks. MeLSI can be extended using the same mechanism by implementing block permutations and passing the strata argument to `adonis2`. The current implementation uses unrestricted permutations for independent group comparisons, which represent the majority of microbiome studies. Paired design support would require separate validation and is listed as an immediate future extension.

**Location:** Limitations and future work section (lines 558-584, specifically lines 573-575)

---

## Closing Remarks

We believe these revisions substantially strengthen the manuscript through expanded simulation-based validation (over 1,675 total simulations), additional real-world validation (SKIOME dataset), and clarified methodological positioning. We remain available for any further discussion the editor or reviewers may find helpful.
