# Response to Reviewers
## MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis
**Manuscript ID:** mSystems01710-25

**Authors:** Nathan Bresette, Aaron Ericsson, Carter Woods, Ai-Ling Lin

---

## Cover Note

We thank the editor and reviewers for their thorough and constructive feedback. We have carefully addressed all concerns raised, with particular attention to the critical statistical validation issues. Below we provide point-by-point responses to each comment.

---

## RESPONSE TO REVIEWER #1

### Major Concerns on Experimental Design

#### Type I Error Control
**Reviewer concern:** Type I error control appears to be evaluated using a single synthetic dataset, with conclusions drawn from a single p-value. Repeated simulations are necessary to support claims of proper error control. Type I error control is a distributional property that must be assessed over repeated realizations of the null hypothesis.

**Response:**
[TBD - Address with expanded Type I error analysis]

#### Statistical Power Analysis
**Reviewer concern:** Similar concerns apply to the statistical power analysis, which also appears to rely on a single dataset per setting.

**Response:**
[TBD - Address with expanded power analysis across multiple simulations]

#### Sample Size Exploration
**Reviewer concern:** Exploring multiple sample sizes in both Type I error and power analyses would strengthen the evaluation. (may consider a few typical sample sizes in microbiome studies)

**Response:**
[TBD - Address with systematic sample size variation]

#### Feature Correlation in Simulations
**Reviewer concern:** Synthetic data generation assumed independent taxa; given the multivariate and correlated structure of microbiome data. Simulations incorporating feature correlation would be important.

**Response:**
[TBD - Address with correlated feature simulations or discussion]

### Feature Pre-filtering

**Reviewer concern:** The pre-filtering strategy relies on t-statistic-like measures for two-group comparisons and ANOVA F-statistics for multi-group settings. Given the count-based, compositional nature of microbiome data, it would be helpful for the authors to discuss the implications of using statistics that assume Gaussian distributions.

**Response:**
[TBD - Add discussion of Gaussian assumption implications]

### Multi-Group Extensions

**Reviewer concern:** The manuscript does not illustrate multi-group extensions using synthetic data, which would help assess generalizability and statistical properties beyond the two-group setting.

**Response:**
[TBD - Address with multi-group synthetic validation or acknowledge as limitation]

### Double Dipping / Overfitting Concerns

**Reviewer concern:** The distance metric is learned in a supervised manner and subsequently used for hypothesis testing on the same dataset. This raises a concern regarding "double dipping." While bagging is employed to mitigate overfitting and the authors claim the proposed method prevents overfitting across the manuscript, the potential impact of in-sample estimation and inference remains unclear and should be discussed more explicitly or addressed through additional validation strategies.

**Response:**
[TBD - Clarify that metric is relearned on each permutation, address overfitting explicitly]

**Reviewer concern:** Statements regarding overfitting prevention across the manuscript would benefit from explicit evaluation or supporting evidence.

**Response:**
[TBD - Provide explicit evaluation of overfitting prevention]

### Ensemble Learning Strategy

**Reviewer concern:** The ensemble learning strategy is motivated by robustness and overfitting prevention, but no direct comparison is provided with existing metric learning approach based on a single-learner.

**Response:**
[TBD - Add comparison or justify ensemble approach]

### Interpretability Advantage

**Reviewer concern:** Emphasizing the interpretability advantage of the proposed method, particularly in identifying taxa that drive group differences, would improve positioning of the manuscript. It would strengthen the manuscript to examine how learned weights recover true signal taxa, possibly across varying effect sizes, sample sizes, and feature dimensions.

**Response:**
[TBD - Add validation of interpretability/recovery of true signal taxa]

### Real Data Analyses

#### Overfitting on Real Data
**Reviewer concern:** The observed outperformance of MeLSI on real datasets may partly reflect overfitting or double-dipping, despite the ensemble strategy. This possibility should be discussed more explicitly.

**Response:**
[TBD - Explicitly discuss this possibility]

#### Figure 2 Concerns
**Reviewer concern:** In Figure 2, the manuscript describes a "clear separation" between male and female samples along PCoA1; however, substantial overlap between groups appears to remain. Can authors elaborate this? Is this separation clearer compared to the results of traditional approaches?

**Response:**
[TBD - Clarify separation claims, compare to traditional methods]

**Reviewer concern:** Could authors explain the rationale for presenting 68% confidence ellipses, rather than more conventional choice of 95%?

**Response:**
[TBD - Explain rationale for 68% ellipses]

**Reviewer concern:** The significant PERMANOVA result (F=5.141, p=0.005) may be influenced by the large sample size (N>1000). I was wondering if the seemingly subtle separation shown in the figure is considered biologically meaningful besides statistical significance.

**Response:**
[TBD - Discuss biological vs statistical significance]

**Reviewer concern:** There is an inconsistency in Figure 2, where PCoA1 variance is listed as 18.4% in the legend but 21.5% on the x-axis.

**Response:**
[TBD - Fix Figure 2 inconsistency]

#### DietSwap Results
**Reviewer concern:** Transition to the DietSwap data analysis is abrupt, and corresponding results including VIP feature and PCoA plots are not shown.

**Response:**
[TBD - Add DietSwap VIP and PCoA plots, improve transition]

### Ensemble Size Analysis

**Reviewer concern:** In the ensemble size analysis, the authors refer to "variance" in performance, but Table 4 does not report variance estimates. If results are based on a single synthetic dataset, clarification of what is meant by variance is needed.

**Response:**
[TBD - Clarify variance reporting in Table 4]

### Pre-filtering Analysis

**Reviewer concern:** The sample sizes used in the synthetic pre-filtering analyses are not clearly stated.

**Response:**
[TBD - Clearly state sample sizes in pre-filtering analysis]

**Reviewer concern:** Feature dimensionality differs across effect size settings (p=500 for small, 200 for medium, and 100 for large). Could authors provide the details on this design choice to help interpretation?

**Response:**
[TBD - Explain design choice for varying dimensionality]

**Reviewer concern:** As with other simulation sections, results based on a single synthetic dataset cannot be interpreted as statistical power.

**Response:**
[TBD - Address with multiple simulations]

### Directionality and Log2 Fold-Change

**Reviewer concern:** The manuscript states that MeLSI provides directionality and log2 fold-change information for each taxon. Additional details on how this information is derived would be helpful.

**Response:**
[TBD - Add detailed explanation of directionality/log2FC calculation]

### Computational Time Justification

**Reviewer concern:** The claim that increased computational time is justified by improved statistical power is not fully supported, especially given inadequate validation study designs and given that Table 2 shows performance comparable to the best traditional methods.

**Response:**
[TBD - Revise justification for computational time, acknowledge limitations]

### Prevalence Threshold

**Reviewer concern:** The mention of a 10% prevalence threshold appears only in the computational efficiency section. Clarification is needed on whether prevalence filtering is an expected preprocessing step and how it plays a role along with the pre-filtering process.

**Response:**
[TBD - Clarify prevalence filtering role and relationship to pre-filtering]

### Tables

**Reviewer concern:** Tables would benefit from clearer annotations and definitions of abbreviations.

**Response:**
[TBD - Add clear annotations and abbreviation definitions to all tables]

### Notation

**Reviewer concern:** Introduce notation more systematically in the "Metric learning: an emerging paradigm" section.

**Response:**
[TBD - Systematically introduce notation]

### CLR Transformation

**Reviewer concern:** The role and impact of the CLR transformation in conjunction with the proposed method should be discussed more explicitly. CLR transformation for MeLSI is mentioned only in the software availability section as part of the recommended usage.

**Response:**
[TBD - Add explicit discussion of CLR transformation role and impact]

---

## RESPONSE TO REVIEWER #3

### Computational Cost vs. Benefit Trade-off

**Reviewer concern:** MeLSI looks slower (minutes to hours) than traditional PERMANOVA (seconds). The paper argues this is acceptable given the gain in interpretability, but for large-scale or screening studies, this may be prohibitive. No clear power-time trade-off analysis is provided to justify when MeLSI is worth the extra computation.

**Response:**
[TBD - Add power-time trade-off analysis or better justification]

### Limited Real-Data Performance Gains

**Reviewer concern:** On Atlas1006, MeLSI's F-statistic was only 9.1% higher than Euclidean distance (5.141 vs. 4.711)-a modest improvement. On DietSwap, MeLSI reached significance (p=0.015) while Bray-Curtis was marginal (p=0.058), but the effect size difference is small. The authors acknowledge that on synthetic data with large effects, traditional metrics (Bray-Curtis, UniFrac) often outperformed MeLSI.

**Response:**
[TBD - Better position MeLSI's value proposition, acknowledge when traditional methods are preferable]

### CLR Transformation May Limit Sensitivity

**Reviewer concern:** MeLSI uses CLR-transformed data for Euclidean distance, which may attenuate large fold-change signals compared to count-based metrics (Bray-Curtis, UniFrac). This could explain why MeLSI underperforms on synthetic data with strong effects.

**Response:**
[TBD - Discuss CLR trade-offs explicitly, explain when CLR approach is appropriate vs. count-based]

### Lack of Covariate Adjustment

**Reviewer concern:** MeLSI currently only supports simple group comparisons. No ability to adjust for confounders (age, BMI, antibiotics) or model continuous outcomes. This limits its utility in observational studies where confounding is common.

**Response:**
[TBD - Acknowledge limitation, discuss as future work]

### Sparse and Compositional Data Challenges

**Reviewer concern:** While pre-filtering is used, the method does not explicitly model compositionality or zero-inflation beyond CLR. Aitchison geometry or Dirichlet-multinomial models might be more appropriate for compositional data.

**Response:**
[TBD - Discuss current approach and acknowledge potential alternatives as future directions]

### Validation on Limited Real Datasets

**Reviewer concern:** Only two real datasets (Atlas1006, DietSwap) were used for benchmarking, both from gut microbiome studies. Performance in other body sites (oral, skin) or disease cohorts is unknown.

**Response:**
[TBD - Acknowledge limitation, justify current validation scope]

### Parameter Sensitivity and Defaults

**Reviewer concern:** Sensitivity analysis shows robustness, but defaults (e.g., B=30, m_frac=0.8) are not rigorously justified. For some datasets, different hyperparameters might yield better performance.

**Response:**
[TBD - Better justify defaults or provide guidance on hyperparameter selection]

### Comparison to Other ML Methods

**Reviewer concern:** No comparison to other interpretable ML models (e.g., Random Forest feature importance, logistic regression with regularization) that also provide taxa rankings. MeLSI is presented as unique, but similar insights might be obtained from simpler models. If authors do so, the ground truth should be well defined for benchmarking.

**Response:**
[TBD - Add comparison to other ML methods or justify why not included]

### Scalability to Very High Dimensions

**Reviewer concern:** The method scales as O(p²) in the number of features. For shotgun metagenomics (thousands of species or genes), this could be computationally prohibitive even with pre-filtering.

**Response:**
[TBD - Discuss scalability limits, acknowledge constraints]

### No Longitudinal/Paired Design Support

**Reviewer concern:** MeLSI does not currently support repeated measures or paired samples, which are common in microbiome intervention studies.

**Response:**
[TBD - Acknowledge limitation, discuss as future work]

---

## SUMMARY OF CHANGES

### Major Additions:
1. [TBD - Expanded Type I error validation with 100+ simulations]
2. [TBD - Expanded power analysis with multiple simulations per condition]
3. [TBD - Sample size exploration in validation]
4. [TBD - DietSwap VIP and PCoA figures]
5. [TBD - Interpretation/recovery validation]

### Major Revisions:
1. [TBD - Clarify double-dipping/overfitting prevention]
2. [TBD - Fix Figure 2 inconsistencies and improve descriptions]
3. [TBD - Explicit CLR transformation discussion]
4. [TBD - Better justification for computational costs]
5. [TBD - Improved limitations discussion]

### Minor Revisions:
1. [TBD - Table annotations and abbreviations]
2. [TBD - Systematic notation introduction]
3. [TBD - Pre-filtering assumptions discussion]
4. [TBD - Directionality calculation details]

---

**Manuscript pages revised:** [TBD]
**New figures added:** [TBD]
**New tables added:** [TBD]
**New supplementary materials:** [TBD]
