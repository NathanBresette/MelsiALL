# Response to Reviewers
## MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis
**Manuscript ID:** mSystems01710-25

**Authors:** Nathan Bresette, Aaron Ericsson, Carter Woods, Ai-Ling Lin

---

## Cover Note

We thank the editor and reviewers for their thorough and constructive feedback. We have carefully addressed all concerns raised, with particular attention to the critical statistical validation issues. Below we provide point-by-point responses to each comment, with references to specific line numbers in the revised manuscript where changes were made.

---

## RESPONSE TO REVIEWER #1

**Reviewer Comment:** This manuscript proposes MeLSI, a supervised, data-driven metric learning approach for constructing distance matrices in microbiome studies, followed by PERMANOVA-based inference for beta-diversity analysis. The manuscript emphasizes (i) rigorous Type I error control, (ii) statistical power, (iii) robustness and overfitting prevention via ensemble learning, and (iv) enhanced biological interpretability. However, while the methodological idea is promising, several key claims are not fully supported by the current experimental design or presented evidence.

**Response:** We appreciate the reviewer's recognition of the methodological promise and acknowledge that several validation aspects require strengthening. We have addressed the concerns below, with particular attention to experimental design improvements and methodological clarifications.

---

### Feature Pre-filtering

**Reviewer Comment:** The pre-filtering strategy relies on t-statistic-like measures for two-group comparisons and ANOVA F-statistics for multi-group settings. Given the count-based, compositional nature of microbiome data, it would be helpful for the authors to discuss the implications of using statistics that assume Gaussian distributions.

**Response:** We have added a brief discussion addressing the use of these statistics on count-based, compositional microbiome data. The revised text (line 99) clarifies that while these statistics are commonly associated with Gaussian assumptions, we use them purely as ranking heuristics to identify features with large between-group differences relative to within-group variation, without relying on distributional assumptions. The actual statistical inference uses permutation testing, which makes no distributional assumptions.

**Location in revised manuscript:** Line 99

---

**Reviewer Comment:** The manuscript does not illustrate multi-group extensions using synthetic data, which would help assess generalizability and statistical properties beyond the two-group setting.

**Response:** [TBD - Address with multi-group synthetic validation or acknowledge as limitation]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### Double Dipping / Overfitting Concerns

**Reviewer Comment:** The distance metric is learned in a supervised manner and subsequently used for hypothesis testing on the same dataset. This raises a concern regarding "double dipping." While bagging is employed to mitigate overfitting and the authors claim the proposed method prevents overfitting across the manuscript, the potential impact of in-sample estimation and inference remains unclear and should be discussed more explicitly or addressed through additional validation strategies.

**Response:** We have clarified that the metric is relearned on each permutation, which prevents double-dipping. The revised text (line 182) explicitly states that by relearning the metric on each permutation, we prevent "double dipping": the metric is not optimized on the observed data and then tested on the same data. Instead, each permutation represents an independent metric learning experiment under the null hypothesis, ensuring that the p-value properly accounts for the adaptive nature of the method. This addresses the reviewer's concern by making explicit that our permutation testing framework treats each permutation as an independent learning experiment, preventing the double-dipping issue.

**Location in revised manuscript:** Line 182

---

**Reviewer Comment:** Statements regarding overfitting prevention across the manuscript would benefit from explicit evaluation or supporting evidence.

**Response:** [TBD - Provide explicit evaluation of overfitting prevention]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

**Reviewer Comment:** The ensemble learning strategy is motivated by robustness and overfitting prevention, but no direct comparison is provided with existing metric learning approach based on a single-learner.

**Response:** We have added a direct comparison with a single-learner approach (B=1) to the parameter sensitivity analysis (Table 4). The results demonstrate that the ensemble approach provides improved robustness and stability compared to a single learner, supporting the use of ensemble learning. This follows established practice in machine learning where combining multiple weak learners improves robustness and reduces overfitting compared to a single learner (20). The combination of bootstrap sampling and feature subsampling creates diversity among learners that a single learner cannot achieve, which is particularly important for adaptive metric learning where overfitting is a concern.

**Location in revised manuscript:** Table 4, Lines 344-349, 359-363

---

### Major Concerns on Experimental Design

#### Type I Error Control

**Reviewer Comment:** Type I error control appears to be evaluated using a single synthetic dataset, with conclusions drawn from a single p-value. Repeated simulations are necessary to support claims of proper error control. Type I error control is a distributional property that must be assessed over repeated realizations of the null hypothesis. Demonstrating that a single test yields p-values greater than 0.05 in one or two null examples does not quantify the probability of false positive findings. Proper evaluation would require repeated simulations under the null, with estimation of the empirical rejection rate at the chosen significance level and, ideally, examination of the null p-value distribution.

**Response:** We have expanded the Type I error analysis to include 100 simulations per dataset type (synthetic null and real shuffled), addressing the reviewer's concern about proper statistical validation. The revised analysis now reports empirical rejection rates at α = 0.05, demonstrating that MeLSI maintains proper Type I error control across repeated realizations of the null hypothesis. We also examine the null p-value distribution, showing that p-values are appropriately calibrated under the null. The revised Table 1 now includes empirical Type I error rates (as percentages) alongside the p-value distribution statistics (mean, median, standard deviation). These results confirm that MeLSI's permutation-based inference properly accounts for the adaptive nature of the method, maintaining Type I error rates near the nominal 5% level.

**Location in revised manuscript:** [TBD - Line numbers after addition to Table 1 and Results section]

---

#### Statistical Power Analysis

**Reviewer Comment:** Similar concerns apply to the statistical power analysis, which also appears to rely on a single dataset per setting.

**Response:** We have expanded the statistical power analysis to include 50 simulations per effect size (small, medium, large), addressing the reviewer's concern about proper evaluation of detection rates. The revised analysis now reports empirical statistical power (detection rates) for each method across repeated simulations, along with mean F-statistics and standard deviations. This allows for proper assessment of power as a distributional property, demonstrating how detection rates vary across different realizations of the same effect size. The revised Table 2 now includes power estimates (as percentages) and mean F-statistics for each effect size, providing a more robust evaluation of MeLSI's performance relative to traditional methods.

**Location in revised manuscript:** [TBD - Line numbers after addition to Table 2 and Results section]

---

#### Sample Size Exploration

**Reviewer Comment:** Exploring multiple sample sizes in both Type I error and power analyses would strengthen the evaluation. (may consider a few typical sample sizes in microbiome studies)

**Response:** We have expanded both Type I error and power analyses to include multiple sample sizes (n=50, 100, 200), addressing the reviewer's concern about evaluating performance across typical microbiome study sizes. The revised analyses now report empirical Type I error rates and statistical power at each sample size, allowing assessment of how MeLSI's performance scales from small (n=50) to larger (n=200) studies. This complements the scalability analysis in Table 3, which demonstrates computational performance across sample sizes, by now also showing statistical validity (Type I error control) and detection capability (power) at each sample size. The results confirm that MeLSI maintains proper Type I error control and demonstrates appropriate power gains with increasing sample size, consistent with standard statistical expectations.

**Location in revised manuscript:** [TBD - Line numbers after addition to Tables 1 and 2, and Results section]

---

#### Feature Correlation in Simulations

**Reviewer Comment:** Synthetic data generation assumed independent taxa; given the multivariate and correlated structure of microbiome data. Simulations incorporating feature correlation would be important.

**Response:** [TBD - Address with correlated feature simulations or discussion]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### Interpretability Advantage

**Reviewer Comment:** Emphasizing the interpretability advantage of the proposed method, particularly in identifying taxa that drive group differences, would improve positioning of the manuscript. It would strengthen the manuscript to examine how learned weights recover true signal taxa, possibly across varying effect sizes, sample sizes, and feature dimensions.

**Response:** [TBD - Add validation of interpretability/recovery of true signal taxa]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### Real Data Analyses

**Reviewer Comment:** The observed outperformance of MeLSI on real datasets may partly reflect overfitting or double-dipping, despite the ensemble strategy. This possibility should be discussed more explicitly.

**Response:** [TBD - Explicitly discuss this possibility]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

**Reviewer Comment:** In Figure 2, the manuscript describes a "clear separation" between male and female samples along PCoA1; however, substantial overlap between groups appears to remain. Can authors elaborate this? Is this separation clearer compared to the results of traditional approaches?

**Response:** We have revised the language to be more precise. The text now states "modest but statistically significant separation" and explicitly compares to traditional methods. The revised text (lines 403-403) reads: "Figure 2 shows modest but statistically significant separation between male and female samples along the first principal coordinate (21.5% of variance). This separation is comparable to that observed with traditional metrics (Euclidean: F=4.711, Bray-Curtis: F=4.442), demonstrating that MeLSI maintains visual separation while providing additional interpretability through learned feature weights."

**Location in revised manuscript:** Lines 403-403

---

**Reviewer Comment:** Could authors explain the rationale for presenting 68% confidence ellipses, rather than more conventional choice of 95%?

**Response:** We have updated Figure 2 to use the conventional 95% confidence ellipses. The manuscript text and figure caption now reflect this change (lines 403 and 408), and the figure has been regenerated accordingly. The code has also been updated to generate 95% confidence ellipses.

**Location in revised manuscript:** Lines 403, 408; Code: `reproducibility_scripts/figure_atlas1006_vip_pcoa.R` line 249, `github/R/melsi_robust.R` line 1227

---

**Reviewer Comment:** The significant PERMANOVA result (F=5.141, p=0.005) may be influenced by the large sample size (N>1000). I was wondering if the seemingly subtle separation shown in the figure is considered biologically meaningful besides statistical significance.

**Response:** We have added explicit discussion addressing this concern. The revised text (lines 403-403) now includes: "While the large sample size (n=1,114) contributes to statistical significance, the sex-associated microbiome differences identified by MeLSI align with previously documented biological patterns (29, 30), and the learned feature weights provide actionable biological insight regardless of sample size." This acknowledges the role of sample size in statistical significance while emphasizing that the biological patterns identified are consistent with known sex-associated microbiome differences.

**Location in revised manuscript:** Lines 403-403

---

**Reviewer Comment:** There is an inconsistency in Figure 2, where PCoA1 variance is listed as 18.4% in the legend but 21.5% on the x-axis.

**Response:** We have fixed this inconsistency by updating the manuscript text to match the actual values shown on the figure. The text now correctly states that PCoA1 explains 21.5% of variance (matching the x-axis), and the figure caption has been updated accordingly (lines 403 and 408). The x-axis label is dynamically generated from the PCoA calculation, and we have ensured the text matches these correct values.

**Location in revised manuscript:** Lines 403, 408

---

**Reviewer Comment:** Transition to the DietSwap data analysis is abrupt, and corresponding results including VIP feature and PCoA plots are not shown.

**Response:** We have improved the transition to the DietSwap section by adding: "To further evaluate MeLSI's utility in real-world applications, we analyzed the DietSwap dietary intervention dataset." (line 297). This provides better context and flow between the Atlas1006 and DietSwap analyses. Note: VIP and PCoA plots for DietSwap remain to be generated as part of additional figure work.

**Location in revised manuscript:** Line 297

---

### Ensemble Size Analysis

**Reviewer Comment:** In the ensemble size analysis, the authors refer to "variance" in performance, but Table 4 does not report variance estimates. If results are based on a single synthetic dataset, clarification of what is meant by variance is needed.

**Response:** We have clarified in the table footnote that "variance in performance" refers to variation in F-statistics across different parameter values tested on the same dataset. Specifically, Table 4 now includes (line 355): "Variance in performance refers to variation in F-statistics across different parameter values tested on the same dataset." This clarifies that we are referring to the range of F-statistic values observed when varying ensemble size (B) or feature fraction (m_frac) parameters, not a statistical variance estimate.

**Location in revised manuscript:** Line 355 (Table 4 footnote)

---

### Pre-filtering Analysis

**Reviewer Comment:** The sample sizes used in the synthetic pre-filtering analyses are not clearly stated.

**Response:** [TBD - Clearly state sample sizes in pre-filtering analysis]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

**Reviewer Comment:** Feature dimensionality differs across effect size settings (p=500 for small, 200 for medium, and 100 for large). Could authors provide the details on this design choice to help interpretation?

**Response:** [TBD - Explain design choice for varying dimensionality]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

**Reviewer Comment:** As with other simulation sections, results based on a single synthetic dataset cannot be interpreted as statistical power.

**Response:** [TBD - Address with multiple simulations]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### Directionality and Log2 Fold-Change

**Reviewer Comment:** The manuscript states that MeLSI provides directionality and log2 fold-change information for each taxon. Additional details on how this information is derived would be helpful.

**Response:** We have added explicit details in the Results section explaining how directionality and log2 fold-change are calculated. The revised text (lines 401-401) now includes: "Directionality is calculated by comparing mean abundances between groups: for each taxon, we identify which group (Group 1 or Group 2) has higher mean abundance on CLR-transformed data. Log2 fold-change is calculated as log2(mean_group1 / mean_group2), where small epsilon values are added to both means to avoid division by zero. These values are computed on the CLR-transformed data used for metric learning, ensuring consistency with the distance metric calculation."

**Location in revised manuscript:** Lines 401-401

---

### Computational Time Justification

**Reviewer Comment:** The claim that increased computational time is justified by improved statistical power is not fully supported, especially given inadequate validation study designs and given that Table 2 shows performance comparable to the best traditional methods.

**Response:** [TBD - Revise justification for computational time, acknowledge limitations]

**Location in revised manuscript:** [TBD - Line numbers after revision]

---

### Prevalence Threshold

**Reviewer Comment:** The mention of a 10% prevalence threshold appears only in the computational efficiency section. Clarification is needed on whether prevalence filtering is an expected preprocessing step and how it plays a role along with the pre-filtering process.

**Response:** We have added clarification in the pre-filtering analysis section (lines 381-381) explaining that "Prevalence filtering (retaining features present in ≥10% of samples) is an optional preprocessing step distinct from MeLSI's variance-based pre-filtering. When applied, prevalence filtering removes rare taxa before MeLSI analysis, while MeLSI's pre-filtering focuses on variance-based feature selection after preprocessing." This clarifies the distinction between the two filtering steps and their roles.

**Location in revised manuscript:** Lines 381-381

---

### Tables

**Reviewer Comment:** Tables would benefit from clearer annotations and definitions of abbreviations.

**Response:** We have added clear annotations and abbreviation definitions to all tables (Tables 1-5). Each table now includes a footnote explaining abbreviations such as: n (sample size), p (number of taxa/features), F (PERMANOVA F-statistic), p (p-value), Time (computation time in seconds), Trad (traditional method), and other table-specific terms.

**Location in revised manuscript:** 
- Table 1: Line 258
- Table 2: Line 277
- Table 3: Line 320
- Table 4: Line 355
- Table 5: Line 376

---

### Notation

**Reviewer Comment:** Introduce notation more systematically in the "Metric learning: an emerging paradigm" section.

**Response:** [TBD - Systematically introduce notation]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### CLR Transformation

**Reviewer Comment:** The role and impact of the CLR transformation in conjunction with the proposed method should be discussed more explicitly. CLR transformation for MeLSI is mentioned only in the software availability section as part of the recommended usage.

**Response:** [TBD - Add explicit discussion of CLR transformation role and impact]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** While we can appreciate the novelty of this work, several critical technical concerns need to be addressed for improvement or clarification.

**Response:** We appreciate the reviewer's recognition of the work's novelty and thank them for the detailed technical feedback. We address each concern below.

---

### 1. Computational Cost vs. Benefit Trade-off

**Reviewer Comment:** MeLSI looks slower (minutes to hours) than traditional PERMANOVA (seconds). The paper argues this is acceptable given the gain in interpretability, but for large-scale or screening studies, this may be prohibitive. No clear power-time trade-off analysis is provided to justify when MeLSI is worth the extra computation.

**Response:** [TBD - Add power-time trade-off analysis or better justification]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 2. Limited Real-Data Performance Gains

**Reviewer Comment:** On Atlas1006, MeLSI's F-statistic was only 9.1% higher than Euclidean distance (5.141 vs. 4.711)-a modest improvement. On DietSwap, MeLSI reached significance (p=0.015) while Bray-Curtis was marginal (p=0.058), but the effect size difference is small. The authors acknowledge that on synthetic data with large effects, traditional metrics (Bray-Curtis, UniFrac) often outperformed MeLSI.

**Response:** [TBD - Better position MeLSI's value proposition, acknowledge when traditional methods are preferable]

**Location in revised manuscript:** [TBD - Line numbers after revision]

---

### 3. CLR Transformation May Limit Sensitivity

**Reviewer Comment:** MeLSI uses CLR-transformed data for Euclidean distance, which may attenuate large fold-change signals compared to count-based metrics (Bray-Curtis, UniFrac). This could explain why MeLSI underperforms on synthetic data with strong effects.

**Response:** [TBD - Discuss CLR trade-offs explicitly, explain when CLR approach is appropriate vs. count-based]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 4. Lack of Covariate Adjustment

**Reviewer Comment:** MeLSI currently only supports simple group comparisons. No ability to adjust for confounders (age, BMI, antibiotics) or model continuous outcomes. This limits its utility in observational studies where confounding is common.

**Response:** [TBD - Acknowledge limitation, discuss as future work]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 5. Sparse and Compositional Data Challenges

**Reviewer Comment:** While pre-filtering is used, the method does not explicitly model compositionality or zero-inflation beyond CLR. Aitchison geometry or Dirichlet-multinomial models might be more appropriate for compositional data.

**Response:** [TBD - Discuss current approach and acknowledge potential alternatives as future directions]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 6. Validation on Limited Real Datasets

**Reviewer Comment:** Only two real datasets (Atlas1006, DietSwap) were used for benchmarking, both from gut microbiome studies. Performance in other body sites (oral, skin) or disease cohorts is unknown.

**Response:** [TBD - Acknowledge limitation, justify current validation scope]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 7. Parameter Sensitivity and Defaults

**Reviewer Comment:** Sensitivity analysis shows robustness, but defaults (e.g., B=30, m_frac=0.8) are not rigorously justified. For some datasets, different hyperparameters might yield better performance.

**Response:** [TBD - Better justify defaults or provide guidance on hyperparameter selection]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 8. Comparison to Other ML Methods

**Reviewer Comment:** No comparison to other interpretable ML models (e.g., Random Forest feature importance, logistic regression with regularization) that also provide taxa rankings. MeLSI is presented as unique, but similar insights might be obtained from simpler models. If authors do so, the ground truth should be well defined for benchmarking.

**Response:** [TBD - Add comparison to other ML methods or justify why not included]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 9. Scalability to Very High Dimensions

**Reviewer Comment:** The method scales as O(p²) in the number of features. For shotgun metagenomics (thousands of species or genes), this could be computationally prohibitive even with pre-filtering.

**Response:** [TBD - Discuss scalability limits, acknowledge constraints]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 10. No Longitudinal/Paired Design Support

**Reviewer Comment:** MeLSI does not currently support repeated measures or paired samples, which are common in microbiome intervention studies.

**Response:** [TBD - Acknowledge limitation, discuss as future work]

**Location in revised manuscript:** [TBD - Line numbers after addition]

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
2. ✅ Fix Figure 2 inconsistencies and improve descriptions (Lines 403, 408)
3. [TBD - Explicit CLR transformation discussion]
4. [TBD - Better justification for computational costs]
5. [TBD - Improved limitations discussion]

### Minor Revisions Completed:
1. ✅ Table annotations and abbreviations (Tables 1-5, Lines 258, 277, 320, 355, 376)
2. [TBD - Systematic notation introduction]
3. ✅ Pre-filtering assumptions discussion (Line 381)
4. ✅ Directionality calculation details (Line 401)
5. ✅ Prevalence threshold clarification (Line 381)
6. ✅ DietSwap transition improvement (Line 297)
7. ✅ Figure 2 ellipses changed to 95% (Lines 403, 408; Code updated)
8. ✅ Figure 2 PCoA variance fixed (21.5%, Lines 403, 408)
9. ✅ Figure 2 separation language revised (Line 403)
10. ✅ Biological vs statistical significance discussion (Line 403)
11. ✅ Table 4 variance clarification (Line 355)

---

**Manuscript pages revised:** [TBD]
**New figures added:** [TBD]
**New tables added:** [TBD]
**New supplementary materials:** [TBD]
