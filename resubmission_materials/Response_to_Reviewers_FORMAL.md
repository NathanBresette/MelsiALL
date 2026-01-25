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

**Response:** We have added explicit evaluation of overfitting prevention using two complementary lines of evidence. First, the comparison with a single-learner baseline (B=1) in Table 4 demonstrates that ensemble learning substantially reduces variance: the single-learner approach shows 4× higher variance (SD = 0.505) compared to ensemble approaches (SD = 0.119-0.128), providing direct quantitative evidence of overfitting prevention. Second, proper Type I error control across 100 simulations (3-6% rejection rates, Table 1) confirms that overfitting does not inflate false positive rates, as inflated Type I error would be expected if overfitting occurred. The permutation testing framework, which relearns the metric on each permutation, ensures that the null distribution properly accounts for the adaptive nature of the method, preventing overfitting from affecting statistical validity.

**Location in revised manuscript:** Parameter Sensitivity section after Table 4 (line 386), Table 4 (line 360), and Table 1 (line 254)

---

**Reviewer Comment:** The ensemble learning strategy is motivated by robustness and overfitting prevention, but no direct comparison is provided with existing metric learning approach based on a single-learner.

**Response:** We have added a direct comparison with a single-learner approach (B=1) to the parameter sensitivity analysis (Table 4, line 352). The results demonstrate that the ensemble approach provides improved robustness and stability compared to a single learner: B=1 shows substantially higher variance (SD = 0.505) and higher p-values (mean = 0.421, SD = 0.29) compared to ensemble approaches (SD = 0.119-0.128 for B≥10), supporting the use of ensemble learning. This follows established practice in machine learning where combining multiple weak learners improves robustness and reduces overfitting compared to a single learner (20). The combination of bootstrap sampling and feature subsampling creates diversity among learners that a single learner cannot achieve, which is particularly important for adaptive metric learning where overfitting is a concern.

**Location in revised manuscript:** Table 4 (line 352), Ensemble size section (lines 363-367)

---

### Major Concerns on Experimental Design

#### Type I Error Control

**Reviewer Comment:** Type I error control appears to be evaluated using a single synthetic dataset, with conclusions drawn from a single p-value. Repeated simulations are necessary to support claims of proper error control. Type I error control is a distributional property that must be assessed over repeated realizations of the null hypothesis. Demonstrating that a single test yields p-values greater than 0.05 in one or two null examples does not quantify the probability of false positive findings. Proper evaluation would require repeated simulations under the null, with estimation of the empirical rejection rate at the chosen significance level and, ideally, examination of the null p-value distribution.

**Response:** We have expanded the Type I error analysis to include 100 simulations per condition across three sample sizes (n=50, 100, 200) for both synthetic null data and real shuffled data (600 total simulations: 2 dataset types × 3 sample sizes × 100 simulations), addressing the reviewer's concern about proper statistical validation. The revised analysis now reports empirical rejection rates at α = 0.05, demonstrating that MeLSI maintains proper Type I error control across repeated realizations of the null hypothesis. The revised Table 1 (line 254) now includes empirical Type I error rates (as percentages) for each sample size and dataset type. These results confirm that MeLSI's permutation-based inference properly accounts for the adaptive nature of the method, maintaining Type I error rates near the nominal 5% level (range: 3-6% across all conditions). We have also updated the Conclusions section (line 436) to reflect this rigorous validation, replacing the previous mention of specific p-values from single simulations with a statement about empirical rejection rates across 100 simulations per condition.

**Location in revised manuscript:** Table 1 (line 254), Results section (lines 265-267), and Conclusions section (line 436)

---

#### Statistical Power Analysis

**Reviewer Comment:** Similar concerns apply to the statistical power analysis, which also appears to rely on a single dataset per setting.

**Response:** We have expanded the statistical power analysis to include 50 simulations per condition across three effect sizes (small, medium, large) and three sample sizes (n=50, 100, 200) (450 total simulations: 3 effect sizes × 3 sample sizes × 50 simulations), addressing the reviewer's concern about proper evaluation of detection rates. The revised analysis now reports empirical statistical power (detection rates) for each method across repeated simulations, along with mean F-statistics. This allows for proper assessment of power as a distributional property, demonstrating how detection rates vary across different realizations of the same effect size. The revised Table 2 (line 275) now includes power estimates (as percentages) and mean F-statistics for each effect size and sample size combination, providing a more robust evaluation of MeLSI's performance relative to traditional methods. We also added a supplementary section (lines 285-287) comparing MeLSI to each traditional method individually, showing that MeLSI consistently outperforms Jaccard and Unweighted UniFrac while demonstrating appropriate conservatism for small effects.

**Location in revised manuscript:** Table 2 (line 275), Individual method comparisons section (lines 285-287), and Synthetic power analysis section (lines 289-305)

---

#### Sample Size Exploration

**Reviewer Comment:** Exploring multiple sample sizes in both Type I error and power analyses would strengthen the evaluation. (may consider a few typical sample sizes in microbiome studies)

**Response:** We have expanded both Type I error and power analyses to include multiple sample sizes (n=50, 100, 200), addressing the reviewer's concern about evaluating performance across typical microbiome study sizes. The revised analyses now report empirical Type I error rates and statistical power at each sample size, allowing assessment of how MeLSI's performance scales from small (n=50) to larger (n=200) studies. This complements the scalability analysis in Table 3 (line 317), which demonstrates computational performance across sample sizes, by now also showing statistical validity (Type I error control) and detection capability (power) at each sample size. The results confirm that MeLSI maintains proper Type I error control and demonstrates appropriate power gains with increasing sample size, consistent with standard statistical expectations.

**Location in revised manuscript:** Table 1 (line 254), Table 2 (line 275), Table 3 (line 323), and Results sections (lines 265-267, 289-305, 328-338)

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

**Response:** We have added explicit discussion addressing this concern in the Real Data: Atlas1006 section. The revised text (lines 310-312) explains that MeLSI's outperformance on real datasets does not reflect overfitting because: (1) the permutation testing framework relearns the metric on each permutation, ensuring that the null distribution properly accounts for the adaptive nature of the method, and (2) proper Type I error control on real shuffled data (3-6% rejection rates across 100 simulations, Table 1) confirms that overfitting does not occur, as inflated Type I error would be expected if overfitting inflated false positive rates on real data. The permutation framework treats each permutation as an independent metric learning experiment under the null hypothesis, preventing overfitting from affecting statistical validity.

**Location in revised manuscript:** Real Data: Atlas1006 section (line 313) and Table 1 (line 254)

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

**Response:** We have significantly expanded the parameter sensitivity analysis to include 25 replications per parameter value (11 parameter values × 25 replications = 275 total simulations), addressing the reviewer's concern about variance estimation. The revised Table 4 (line 352) now reports mean F-statistics, p-values, and computation times with standard deviations (SD) across the 25 replications. We also added a direct comparison with a single-learner approach (B=1), which shows substantially higher variance (SD = 0.505) compared to ensemble approaches (SD = 0.119-0.128), supporting the use of ensemble learning. The table footnote (line 361) now clarifies that values are shown as "mean (SD) across 25 replications per parameter value."

**Location in revised manuscript:** Table 4 (line 352), Ensemble size section (lines 363-367)

---

### Pre-filtering Analysis

**Reviewer Comment:** The sample sizes used in the synthetic pre-filtering analyses are not clearly stated.

**Response:** We have clarified that all pre-filtering analyses use n=100 samples, as stated in Table 5 (line 386). The table now explicitly shows the sample size (n=100) for all three effect size conditions.

**Location in revised manuscript:** Table 5 (line 398)

---

**Reviewer Comment:** Feature dimensionality differs across effect size settings (p=500 for small, 200 for medium, and 100 for large). Could authors provide the details on this design choice to help interpretation?

**Response:** We have clarified in the Results section (line 375) that the varying dimensionality (p=500 for small, p=200 for medium, p=100 for large) reflects the design choice to test pre-filtering benefits across different dimensionalities. This allows assessment of how pre-filtering performance scales with dimensionality, with higher-dimensional datasets (p=500) showing greater time savings (39.8%) compared to lower-dimensional datasets (p=100, 16.5% time savings).

**Location in revised manuscript:** Pre-filtering analysis section (line 375), Table 5 (line 386)

---

**Reviewer Comment:** As with other simulation sections, results based on a single synthetic dataset cannot be interpreted as statistical power.

**Response:** We have significantly expanded the pre-filtering analysis to include 50 simulations per effect size (3 effect sizes × 50 simulations = 150 total), addressing the reviewer's concern about proper statistical evaluation. The revised analysis now reports empirical statistical power (detection rates) for both pre-filtered and non-filtered approaches across repeated simulations, along with mean F-statistics and standard deviations. The revised Table 5 (line 398) now includes power estimates, mean F-statistics, and mean time reduction percentages, providing a rigorous evaluation of the pre-filtering strategy's impact on both statistical power and computational efficiency. The results demonstrate substantial benefits: pre-filtering increases power from 4% to 100% for small effects, from 14% to 94% for medium effects, and from 14% to 84% for large effects, while providing 16-40% time savings.

**Location in revised manuscript:** Table 5 (line 398) and Pre-filtering analysis section (lines 388-395)

---

### Directionality and Log2 Fold-Change

**Reviewer Comment:** The manuscript states that MeLSI provides directionality and log2 fold-change information for each taxon. Additional details on how this information is derived would be helpful.

**Response:** We have added explicit details in the Results section explaining how directionality and log2 fold-change are calculated. The revised text (lines 401-401) now includes: "Directionality is calculated by comparing mean abundances between groups: for each taxon, we identify which group (Group 1 or Group 2) has higher mean abundance on CLR-transformed data. Log2 fold-change is calculated as log2(mean_group1 / mean_group2), where small epsilon values are added to both means to avoid division by zero. These values are computed on the CLR-transformed data used for metric learning, ensuring consistency with the distance metric calculation."

**Location in revised manuscript:** Lines 401-401

---

### Computational Time Justification

**Reviewer Comment:** The claim that increased computational time is justified by improved statistical power is not fully supported, especially given inadequate validation study designs and given that Table 2 shows performance comparable to the best traditional methods.

**Response:** We have revised the computational time justification in the Limitations section (lines 456-456) to provide a more comprehensive and evidence-based rationale. The revised text justifies MeLSI's computational time (2-30 minutes for typical datasets) based on: (1) substantial interpretability gains through learned feature weights that identify biologically relevant taxa, (2) pre-filtering benefits that provide 16-40% time savings while improving power by 36-37% (Table 5), and (3) the modest time investment relative to the overall study timeline (weeks to months for sample collection and sequencing). We acknowledge that for large-scale screening studies with thousands of samples, traditional methods may be more appropriate. The justification now emphasizes interpretability as the primary benefit rather than power alone, which aligns with Table 2 showing comparable power to traditional methods while providing unique interpretability advantages.

**Location in revised manuscript:** Limitations section (line 460) and Table 5 (line 398)

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

**Response:** We have systematically introduced notation in the "Metric learning: an emerging paradigm" section (lines 41-43). The revised text now formally defines: (1) the feature abundance matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$ with $n$ samples and $p$ taxa, (2) group labels $\mathbf{y} = (y_1, \ldots, y_n)$, (3) the positive semi-definite metric matrix $\mathbf{M} \in \mathbb{R}^{p \times p}$, (4) the Mahalanobis distance formula $d_M(\mathbf{x}_i, \mathbf{x}_j) = \sqrt{(\mathbf{x}_i - \mathbf{x}_j)^T \mathbf{M} (\mathbf{x}_i - \mathbf{x}_j)}$, and (5) the special case of diagonal $\mathbf{M}$ reducing to weighted Euclidean distance with feature-specific weights $M_{jj}$. This systematic introduction provides a clear mathematical foundation before the notation is used throughout the Methods section.

**Location in revised manuscript:** Introduction section, "Metric learning: an emerging paradigm" subsection (line 43)

---

### CLR Transformation

**Reviewer Comment:** The role and impact of the CLR transformation in conjunction with the proposed method should be discussed more explicitly. CLR transformation for MeLSI is mentioned only in the software availability section as part of the recommended usage.

**Response:** We have added explicit discussion of the CLR transformation role and impact in the Methods section (lines 234-234). The revised text explains that: (1) CLR transformation converts relative abundances to log-ratios, making the data suitable for Euclidean distance while preserving relative relationships between taxa, (2) CLR treats abundance ratios more equitably than count-based metrics, which can be dominated by highly abundant taxa, (3) CLR may attenuate large fold-change signals compared to count-based metrics, as evidenced by results showing traditional count-based methods achieve higher F-statistics on synthetic data with large effects (3× fold change), and (4) CLR is particularly appropriate when signals are distributed across multiple taxa rather than concentrated in highly abundant taxa, and when interpretability through feature weights is prioritized. We also added discussion in the Results section (lines 305-307) explicitly stating when the CLR-based approach is most appropriate versus when count-based methods may be preferable.

**Location in revised manuscript:** Methods section (lines 234-234) and Results section (lines 305-307)

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** While we can appreciate the novelty of this work, several critical technical concerns need to be addressed for improvement or clarification.

**Response:** We appreciate the reviewer's recognition of the work's novelty and thank them for the detailed technical feedback. We address each concern below.

---

### 1. Computational Cost vs. Benefit Trade-off

**Reviewer Comment:** MeLSI looks slower (minutes to hours) than traditional PERMANOVA (seconds). The paper argues this is acceptable given the gain in interpretability, but for large-scale or screening studies, this may be prohibitive. No clear power-time trade-off analysis is provided to justify when MeLSI is worth the extra computation.

**Response:** We have added an explicit power-time trade-off analysis in the Computational Performance section (lines 444-448). The analysis demonstrates that MeLSI provides substantial value: pre-filtering increases statistical power by 36-37% while reducing computation time by 16-40% (Table 5). For typical microbiome studies (n=50-200, p=100-500), MeLSI completes in 2-30 minutes (Table 3), representing a modest time investment that yields both improved power (particularly for medium effect sizes, Table 2) and interpretability through feature weights. The power-time trade-off is most favorable when: (1) sample sizes are moderate (n=50-200), (2) interpretability is prioritized, and (3) pre-filtering is applied. For very large studies (n>500) or when only rapid screening is needed, traditional methods may be preferable. This analysis provides clear guidance on when MeLSI's computational investment is justified.

**Location in revised manuscript:** Computational Performance section, Power-time trade-off analysis subsection (line 446), Tables 2 (line 275), 3 (line 323), and 5 (line 398)

---

### 2. Limited Real-Data Performance Gains

**Reviewer Comment:** On Atlas1006, MeLSI's F-statistic was only 9.1% higher than Euclidean distance (5.141 vs. 4.711)-a modest improvement. On DietSwap, MeLSI reached significance (p=0.015) while Bray-Curtis was marginal (p=0.058), but the effect size difference is small. The authors acknowledge that on synthetic data with large effects, traditional metrics (Bray-Curtis, UniFrac) often outperformed MeLSI.

**Response:** We have revised the Conclusions section (lines 456-456) to explicitly position MeLSI's value proposition and acknowledge when traditional methods are preferable. The revised text clearly states that MeLSI is recommended when: (1) effect sizes are moderate (2× fold change) rather than very large, (2) interpretability through feature weights is needed to identify biologically relevant taxa, (3) traditional methods yield marginal results (p-values near 0.05), and (4) signals are distributed across multiple taxa rather than concentrated in highly abundant taxa. Traditional methods (Bray-Curtis, UniFrac) are preferable for: (1) large, obvious effects (3× fold change) where any method succeeds, (2) large-scale screening studies where speed is critical, and (3) when only omnibus significance testing is needed without feature-level interpretation. This positioning acknowledges that MeLSI's primary value is interpretability rather than marginal power gains, while clearly stating when each approach is most appropriate.

**Location in revised manuscript:** Conclusions section (line 456) and Table 2 (line 275)

---

### 3. CLR Transformation May Limit Sensitivity

**Reviewer Comment:** MeLSI uses CLR-transformed data for Euclidean distance, which may attenuate large fold-change signals compared to count-based metrics (Bray-Curtis, UniFrac). This could explain why MeLSI underperforms on synthetic data with strong effects.

**Response:** We have added explicit discussion of CLR trade-offs in both the Methods section (line 234) and Results section (line 305). The revised text acknowledges that CLR transformation may attenuate large fold-change signals compared to count-based metrics, as evidenced by our results showing that traditional count-based methods achieve higher F-statistics on synthetic data with large effects (3× fold change). We explicitly state when the CLR-based approach is most appropriate: (1) when signals are distributed across multiple taxa rather than concentrated in highly abundant taxa, (2) when interpretability through feature weights is prioritized, and (3) when effect sizes are moderate rather than very large. For large, obvious effects (3× fold change), count-based methods (Bray-Curtis, UniFrac) may be preferable due to their sensitivity to abundance dominance. This positioning acknowledges the trade-off while clarifying when each approach is most appropriate.

**Location in revised manuscript:** Methods section (line 234) and Results section (line 305)

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

**Response:** We have expanded the Parameter Sensitivity section (lines 388-390) to explicitly justify the default parameters based on Table 4 results. The revised text explains that: (1) B=30 provides F-statistics (mean = 1.530, SD = 0.123) comparable to larger ensembles (B=50-100) while maintaining reasonable computation time (mean = 576.8s), (2) m_frac=0.8 balances performance (mean F = 1.530) with diversity among weak learners, as lower values (m_frac=0.5) show slightly higher F-statistics but reduced diversity, while higher values (m_frac=0.9-1.0) show slightly lower F-statistics, and (3) the robustness demonstrated across wide parameter ranges (B=10-100, m_frac=0.5-1.0) indicates that default parameters provide good performance across diverse datasets, though users may optimize for specific datasets if needed. This provides clear justification for the defaults while acknowledging that optimization may be beneficial for specific datasets.

**Location in revised manuscript:** Parameter Sensitivity section (line 388) and Table 4 (line 360)

---

### 8. Comparison to Other ML Methods

**Reviewer Comment:** No comparison to other interpretable ML models (e.g., Random Forest feature importance, logistic regression with regularization) that also provide taxa rankings. MeLSI is presented as unique, but similar insights might be obtained from simpler models. If authors do so, the ground truth should be well defined for benchmarking.

**Response:** [TBD - Add comparison to other ML methods or justify why not included]

**Location in revised manuscript:** [TBD - Line numbers after addition]

---

### 9. Scalability to Very High Dimensions

**Reviewer Comment:** The method scales as O(p²) in the number of features. For shotgun metagenomics (thousands of species or genes), this could be computationally prohibitive even with pre-filtering.

**Response:** We have expanded the Scalability section (lines 352-354) to explicitly discuss scalability limits and acknowledge constraints. The revised text explains that MeLSI's O(p²) scaling becomes computationally prohibitive for very high-dimensional datasets (p>1000), with Table 3 demonstrating that computation time increases from 227.9s at p=50 to 8405.6s at p=1000. However, pre-filtering (retaining 70% of features) substantially mitigates this scaling, reducing effective dimensionality. For shotgun metagenomics with thousands of features, we recommend: (1) applying pre-filtering to reduce dimensionality, (2) considering feature aggregation (e.g., species-level rather than gene-level), or (3) using traditional methods if interpretability is not prioritized. The current implementation is most suitable for typical 16S rRNA datasets (p<1000) and metagenomic datasets with moderate dimensionality after preprocessing. This acknowledges the scalability constraint while providing practical guidance for high-dimensional applications.

**Location in revised manuscript:** Scalability section, Dimensionality scaling subsection (line 354) and Table 3 (line 323)

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
12. ✅ Conclusions section updated to reflect rigorous validation (Line 436)
13. ✅ Parameter sensitivity mention enhanced with B=1 comparison (Line 438)

---

**Manuscript pages revised:** [TBD]
**New figures added:** [TBD]
**New tables added:** [TBD]
**New supplementary materials:** [TBD]
