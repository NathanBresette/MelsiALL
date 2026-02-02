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

**Response:** We have added discussion clarifying that these statistics serve as ranking heuristics, not inferential tests. The revised text (line 99) explains that while t-statistics and F-statistics are commonly associated with Gaussian assumptions, we use them solely to rank features by between-group differences relative to within-group variation, without relying on distributional assumptions. All statistical inference uses permutation testing, which makes no distributional assumptions and is appropriate for count-based, compositional data.

**Location in revised manuscript:** Line 99

---

**Reviewer Comment:** The manuscript does not illustrate multi-group extensions using synthetic data, which would help assess generalizability and statistical properties beyond the two-group setting.

**Response:** We acknowledge this limitation and have added explicit discussion in the Limitations section. The current synthetic validation focuses on two-group comparisons, which represent the primary use case for microbiome beta diversity analysis. Comprehensive multi-group synthetic validation would require duplicating all validation tables (over 1,500 additional simulations), representing substantial computational resources.

We address multi-group generalizability through real-world validation on the SKIOME skin microbiome dataset (PRJNA554499), which includes three groups (Atopic Dermatitis, Healthy, Psoriasis) with 511 samples and 1,856 taxa. MeLSI's omnibus test detected significant differences (F = 4.895, p = 0.005), with all pairwise comparisons remaining significant after FDR correction (p = 0.005 for all pairs). This demonstrates MeLSI's utility beyond two-group comparisons and provides evidence of generalizability to diverse body sites (skin vs. gut microbiome) and disease cohorts. The statistical framework (permutation testing, Type I error control) is identical for two-group and multi-group analyses (Methods section, Multi-group extensions subsection). The permutation approach that relearns the metric on each permutation ensures valid statistical inference regardless of the number of groups, as the null distribution properly accounts for the adaptive nature of the method.

For the current resubmission, we focus on rigorous two-group validation, which allows comprehensive statistical evaluation across multiple effect sizes, sample sizes, and dimensionalities. Multi-group synthetic validation represents an important future enhancement that would complement the real-world SKIOME validation. A brief description of multi-group capabilities is included in the Methods section, with acknowledgment that synthetic validation is not included.

**Location in revised manuscript:** Limitations and future work section (line 445), Methods section (Multi-group extensions subsection, line 229), and Real data validation section (SKIOME dataset subsection, lines 416-427)

---

### Double Dipping / Overfitting Concerns

**Reviewer Comment:** The distance metric is learned in a supervised manner and subsequently used for hypothesis testing on the same dataset. This raises a concern regarding "double dipping." While bagging is employed to mitigate overfitting and the authors claim the proposed method prevents overfitting across the manuscript, the potential impact of in-sample estimation and inference remains unclear and should be discussed more explicitly or addressed through additional validation strategies.

**Response:** We have clarified that the metric is relearned on each permutation, which prevents double-dipping. The revised text (line 182) explicitly states that each permutation represents an independent metric learning experiment under the null hypothesis: the metric is not optimized on the observed data and then tested on the same data. Instead, the permutation framework relearns the metric on each permuted dataset, ensuring that the p-value properly accounts for the adaptive nature of the method. This prevents double-dipping by treating each permutation as an independent learning experiment.

**Location in revised manuscript:** Line 182

---

**Reviewer Comment:** Statements regarding overfitting prevention across the manuscript would benefit from explicit evaluation or supporting evidence.

**Response:** We have added explicit evaluation of overfitting prevention using two complementary lines of evidence. First, Table 4 demonstrates that ensemble learning substantially reduces variance compared to a single-learner baseline (B=1): the single-learner approach shows 4× higher variance (SD = 0.505) compared to ensemble approaches (SD = 0.119-0.128), providing direct quantitative evidence of overfitting prevention. Second, proper Type I error control across 100 simulations (3-6% rejection rates, Table 1) confirms that overfitting does not inflate false positive rates, as inflated Type I error would be expected if overfitting occurred. The permutation testing framework relearns the metric on each permutation, ensuring the null distribution properly accounts for the adaptive nature of the method.

**Location in revised manuscript:** Parameter Sensitivity section after Table 4 (line 338), Table 4 (line 318), and Table 1 (line 254)

---

**Reviewer Comment:** The ensemble learning strategy is motivated by robustness and overfitting prevention, but no direct comparison is provided with existing metric learning approach based on a single-learner.

**Response:** We have added a direct comparison with a single-learner approach (B=1) to the parameter sensitivity analysis (Table 4, line 318). The results demonstrate that the ensemble approach provides improved robustness and stability compared to a single learner: B=1 shows substantially higher variance (SD = 0.505) and higher p-values (mean = 0.421, SD = 0.29) compared to ensemble approaches (SD = 0.119-0.128 for B≥10), supporting the use of ensemble learning. This follows established practice in machine learning where combining multiple weak learners improves robustness and reduces overfitting compared to a single learner (20). The combination of bootstrap sampling and feature subsampling creates diversity among learners that a single learner cannot achieve, which is particularly important for adaptive metric learning where overfitting is a concern.

**Location in revised manuscript:** Table 4 (line 318), Parameter Sensitivity section text (line 338)

---

### Major Concerns on Experimental Design

#### Type I Error Control

**Reviewer Comment:** Type I error control appears to be evaluated using a single synthetic dataset, with conclusions drawn from a single p-value. Repeated simulations are necessary to support claims of proper error control. Type I error control is a distributional property that must be assessed over repeated realizations of the null hypothesis. Demonstrating that a single test yields p-values greater than 0.05 in one or two null examples does not quantify the probability of false positive findings. Proper evaluation would require repeated simulations under the null, with estimation of the empirical rejection rate at the chosen significance level and, ideally, examination of the null p-value distribution.

**Response:** We have expanded the Type I error analysis to include 100 simulations per condition across three sample sizes (n=50, 100, 200) for both synthetic null data and real shuffled data (600 total simulations: 2 dataset types × 3 sample sizes × 100 simulations). The revised analysis reports empirical rejection rates at α = 0.05, demonstrating that MeLSI maintains proper Type I error control across repeated realizations of the null hypothesis. Table 1 (line 254) now includes empirical Type I error rates (as percentages) for each sample size and dataset type. Results confirm that MeLSI's permutation-based inference properly accounts for the adaptive nature of the method, maintaining Type I error rates near the nominal 5% level (range: 3-6% across all conditions). We have updated the Conclusions section, Summary subsection (line 437) to reflect this rigorous validation.

**Location in revised manuscript:** Table 1 (line 254), Results section (lines 265-267), and Conclusions section, Summary subsection (line 437)

---

#### Statistical Power Analysis

**Reviewer Comment:** Similar concerns apply to the statistical power analysis, which also appears to rely on a single dataset per setting.

**Response:** We have expanded the statistical power analysis to include 50 simulations per condition across three effect sizes (small, medium, large) and three sample sizes (n=50, 100, 200) (450 total simulations: 3 effect sizes × 3 sample sizes × 50 simulations). The revised analysis reports empirical statistical power (detection rates) for each method across repeated simulations, along with mean F-statistics, allowing proper assessment of power as a distributional property. Table 2 (line 275) now includes power estimates (as percentages) and mean F-statistics for each effect size and sample size combination. We also added a supplementary section (lines 285-287) comparing MeLSI to each traditional method individually, showing that MeLSI consistently outperforms Jaccard and Unweighted UniFrac while demonstrating appropriate conservatism for small effects.

**Location in revised manuscript:** Table 2 (line 275), Individual method comparisons section (lines 285-287), and Synthetic power analysis section (lines 289-305)

---

#### Sample Size Exploration

**Reviewer Comment:** Exploring multiple sample sizes in both Type I error and power analyses would strengthen the evaluation. (may consider a few typical sample sizes in microbiome studies)

**Response:** We have expanded both Type I error and power analyses to include multiple sample sizes (n=50, 100, 200), evaluating performance across typical microbiome study sizes. The revised analyses now report empirical Type I error rates and statistical power at each sample size, allowing assessment of how MeLSI's performance scales from small (n=50) to larger (n=200) studies. This complements the scalability analysis in Table 3 (line 295), which demonstrates computational performance across sample sizes, by now also showing statistical validity (Type I error control) and detection capability (power) at each sample size. The results confirm that MeLSI maintains proper Type I error control and demonstrates appropriate power gains with increasing sample size, consistent with standard statistical expectations.

**Location in revised manuscript:** Table 1 (line 254), Table 2 (line 275), Table 3 (line 295), and Results sections (lines 265-267, 289-305, 328-338)

---

#### Feature Correlation in Simulations

**Reviewer Comment:** Synthetic data generation assumed independent taxa; given the multivariate and correlated structure of microbiome data. Simulations incorporating feature correlation would be important.

**Response:** We have added comprehensive validation of MeLSI's performance under varying levels of feature correlation (Table 5, line 344). The analysis evaluated four correlation levels (r=0, 0.3, 0.6, 0.8) using 50 simulations per condition (200 total simulations). Results demonstrate that MeLSI is robust to correlated features: F-statistics remained stable across correlation levels (±1.7% variation), statistical power was consistent (42-50%), and feature recovery metrics (Precision at 10, AUC-ROC) showed minimal variation. MeLSI achieved its best ranking (1/6) at high correlation (r=0.8), suggesting the method is particularly well-suited for correlated microbiome data. The ensemble learning approach naturally handles correlated features by aggregating signal across correlated taxa, which is appropriate for microbiome data where correlated taxa often represent biologically related groups.

**Location in revised manuscript:** Table 5 (line 344), Feature Correlation Robustness section (lines 340-356)

---

### Interpretability Advantage

**Reviewer Comment:** Emphasizing the interpretability advantage of the proposed method, particularly in identifying taxa that drive group differences, would improve positioning of the manuscript. It would strengthen the manuscript to examine how learned weights recover true signal taxa, possibly across varying effect sizes, sample sizes, and feature dimensions.

**Response:** We have added comprehensive validation of interpretability and recovery of true signal taxa across varying effect sizes and sample sizes (lines 295-304). The analysis evaluates how well MeLSI's learned feature weights identify true signal taxa using four metrics: Precision at k, Recall at k, Mean Rank, and AUC-ROC. Results demonstrate that MeLSI effectively recovers true signal taxa, with performance improving substantially with effect size and sample size. For large effects, Precision at 5 reached 0.876-1.000 and Mean Rank decreased to 14.4, confirming that true signal taxa are consistently ranked among the top features. These results validate MeLSI's interpretability advantage: the learned feature weights reliably identify biologically relevant taxa that drive group differences, with recovery performance scaling appropriately with signal strength and sample size.

**Location in revised manuscript:** Lines 295-304 (Recovery of true signal taxa subsection)

---

### Real Data Analyses

**Reviewer Comment:** The observed outperformance of MeLSI on real datasets may partly reflect overfitting or double-dipping, despite the ensemble strategy. This possibility should be discussed more explicitly.

**Response:** We have added explicit discussion in the Real Data: Atlas1006 section. The revised text explains that MeLSI's outperformance on real datasets does not reflect overfitting because the permutation testing framework relearns the metric on each permutation, ensuring the null distribution properly accounts for the adaptive nature of the method. This is confirmed by proper Type I error control on real shuffled data (3-6% rejection rates across 100 simulations, Table 1), demonstrating that MeLSI's outperformance reflects genuine signal detection rather than overfitting. The permutation framework treats each permutation as an independent metric learning experiment under the null hypothesis, which serves as the inherent guardrail against overfitting.

**Location in revised manuscript:** Real Data: Atlas1006 section and Table 1 (line 254)

---

**Reviewer Comment:** In Figure 2, the manuscript describes a "clear separation" between male and female samples along PCoA1; however, substantial overlap between groups appears to remain. Can authors elaborate this? Is this separation clearer compared to the results of traditional approaches?

**Response:** We acknowledge the reviewer's valid concern about the visual separation in the Atlas1006 PCoA plot. Given the substantial overlap and the reviewer's question about whether the separation is meaningful, we have removed the standalone Atlas1006 PCoA figure (previously Figure 2) to avoid overstating visual separation. The statistical significance (F=5.141, p=0.005) is already reported in the text, and the VIP plot (Figure 1) provides the primary interpretability value by identifying key taxa. PCoA ordination is a standard visualization tool for traditional distance-based methods (Bray-Curtis, Euclidean, etc.) and works equally well with MeLSI's learned distance matrices. We now show PCoA ordination for DietSwap and SKIOME datasets (Figures 2-3), where group separation is more visually apparent, demonstrating MeLSI's utility for ordination when separation is clearer and reaffirming that MeLSI-learned distances are compatible with standard ordination approaches used throughout the microbiome field.

**Location in revised manuscript:** Figure 1 (Atlas1006 VIP only), Figures 2-3 (DietSwap and SKIOME with combined VIP+PCoA)

---

**Reviewer Comment:** Could authors explain the rationale for presenting 68% confidence ellipses, rather than more conventional choice of 95%?

**Response:** All figures now use the conventional 95% confidence ellipses. The manuscript text and figure captions reflect this change, and all figures have been regenerated accordingly. The code has also been updated to generate 95% confidence ellipses for all PCoA plots.

**Location in revised manuscript:** Figures 2-3 (DietSwap and SKIOME PCoA plots); Code: `github/R/melsi_robust.R` line 1227

---

**Reviewer Comment:** The significant PERMANOVA result (F=5.141, p=0.005) may be influenced by the large sample size (N>1000). I was wondering if the seemingly subtle separation shown in the figure is considered biologically meaningful besides statistical significance.

**Response:** We acknowledge that the large sample size (n=1,114) contributes to statistical significance. The revised text (line 379) explicitly states that "MeLSI's improvement over the best fixed metric suggests that learned metrics can capture biologically relevant patterns in subtle, high-dimensional comparisons, consistent with previously documented sex-associated microbiome differences (29, 30)." This emphasizes that the biological patterns identified are consistent with known sex-associated microbiome differences, and the learned feature weights (Figure 1) provide actionable biological insight regardless of sample size.

**Location in revised manuscript:** Line 379, Figure 1

---

**Reviewer Comment:** There is an inconsistency in Figure 2, where PCoA1 variance is listed as 18.4% in the legend but 21.5% on the x-axis.

**Response:** This issue has been resolved by removing the standalone Atlas1006 PCoA figure (previously Figure 2) in response to the reviewer's concern about visual separation. All remaining PCoA figures (Figures 2-3) have been verified to ensure consistency between x-axis labels (dynamically generated from PCoA calculations) and figure captions.

**Location in revised manuscript:** Figures 2-3 (DietSwap and SKIOME PCoA plots)

---

**Reviewer Comment:** Transition to the DietSwap data analysis is abrupt, and corresponding results including VIP feature and PCoA plots are not shown.

**Response:** We have improved the transition to the DietSwap section by adding: "To evaluate MeLSI's utility in real-world applications, we analyzed three published microbiome datasets: Atlas1006 (sex-associated differences), DietSwap (dietary intervention), and SKIOME (multi-group skin microbiome validation)." (line 375). This provides better context and flow between the Atlas1006 and DietSwap analyses. VIP and PCoA plots for DietSwap have been added to the manuscript (Figure 2, lines 399-407).

**Location in revised manuscript:** Line 375, Figure 2 (lines 399-407)

---

### Ensemble Size Analysis

**Reviewer Comment:** In the ensemble size analysis, the authors refer to "variance" in performance, but Table 4 does not report variance estimates. If results are based on a single synthetic dataset, clarification of what is meant by variance is needed.

**Response:** We have significantly expanded the parameter sensitivity analysis to include 25 replications per parameter value (11 parameter values × 25 replications = 275 total simulations), enabling robust variance estimation. The revised Table 4 (line 318) now reports mean F-statistics, p-values, and computation times across the 25 replications. We also added a direct comparison with a single-learner approach (B=1), which shows substantially higher variance (SD = 0.505) compared to ensemble approaches (SD = 0.119-0.128), supporting the use of ensemble learning. Standard deviations for all parameters are provided in Supplementary Table S4, which demonstrates that ensemble learning reduces variance and prevents overfitting. The analysis shows that F-statistics remained stable across ensemble sizes (B=10-100) compared to the single-learner baseline (B=1). The table footnote (line 336) clarifies that values are shown as mean across 25 replications per parameter value, with SD values in Supplementary Table S4.

**Location in revised manuscript:** Table 4 (line 318), Parameter Sensitivity section text (line 338)

---

### Pre-filtering Analysis

**Reviewer Comment:** The sample sizes used in the synthetic pre-filtering analyses are not clearly stated.

**Response:** We have clarified that all pre-filtering analyses use n=100 samples, as stated in Table 6 (line 361). The table now explicitly shows the sample size (n=100) for all three effect size conditions.

**Location in revised manuscript:** Table 6 (line 361)

---

**Reviewer Comment:** Feature dimensionality differs across effect size settings (p=500 for small, 200 for medium, and 100 for large). Could authors provide the details on this design choice to help interpretation?

**Response:** We have clarified in the Results section (line 375) that the varying dimensionality (p=500 for small, p=200 for medium, p=100 for large) reflects the design choice to test pre-filtering benefits across different dimensionalities. This allows assessment of how pre-filtering performance scales with dimensionality, with higher-dimensional datasets (p=500) showing greater time savings (39.8%) compared to lower-dimensional datasets (p=100, 16.5% time savings).

**Location in revised manuscript:** Pre-filtering analysis section (line 357), Table 6 (line 361)

---

**Reviewer Comment:** As with other simulation sections, results based on a single synthetic dataset cannot be interpreted as statistical power.

**Response:** We have significantly expanded the pre-filtering analysis to include 50 simulations per scenario (150 total), enabling robust statistical power evaluation. The revised Table 6 (line 361) now reports empirical power, mean F-statistics (with standard deviations), and computational time savings.

Our results demonstrate that pre-filtering is most impactful in high-dimensional, sparse-signal settings—the exact conditions where standard microbiome analyses often fail. In the high-dimensional setting (p=500 features, small effect: 1.5× fold change in 5 signal taxa), pre-filtering removed 150 noise features, increasing power from 4% to 100%. In the lower-dimensional setting (p=100 features, large effect: 3.0× fold change in 20 signal taxa), power increased from 14% to 84%. The higher power observed in the small-effect scenario relative to the large-effect scenario is a function of the experimental design: pre-filtering provides maximum benefit when the initial signal-to-noise ratio is lowest (i.e., higher dimensionality). By removing 30% of non-informative features in the p=500 environment, the model successfully recovered signal that was previously masked by noise. This highlights the utility of our approach for typical high-dimensional microbiome datasets.

**Location in revised manuscript:** Table 6 (line 361) and Pre-filtering analysis section (lines 357-373)

---

### Directionality and Log2 Fold-Change

**Reviewer Comment:** The manuscript states that MeLSI provides directionality and log2 fold-change information for each taxon. Additional details on how this information is derived would be helpful.

**Response:** We have clarified the derivation of directionality and effect size in the Results section (lines 395-396). **Directionality:** Is determined by identifying which group has the higher mean abundance on the CLR-transformed data, ensuring consistency with the metric learning process. **Effect Size:** We report the difference in CLR-transformed means between groups. Because CLR data is already in log-ratio space, the difference between means ($\mu_{\text{CLR,1}} - \mu_{\text{CLR,2}}$) is the standard way to represent the log-fold change for compositional data. This ensures that the reported effect sizes are directly derived from the same feature space used to calculate the distances in MeLSI, providing a mathematically consistent representation of the taxa's contribution to group separation.

**Location in revised manuscript:** Lines 395-396

---

### Computational Time Justification

**Reviewer Comment:** The claim that increased computational time is justified by improved statistical power is not fully supported, especially given inadequate validation study designs and given that Table 2 shows performance comparable to the best traditional methods.

**Response:** We have revised the computational time justification in the Limitations section (line 443) to provide a more comprehensive and evidence-based rationale. The revised text justifies MeLSI's computational time (2-30 minutes for typical datasets) based on: (1) substantial interpretability gains through learned feature weights that identify biologically relevant taxa, (2) pre-filtering benefits that provide 16-40% time savings while improving power by 36-37% (Table 6), and (3) the modest time investment relative to the overall study timeline (weeks to months for sample collection and sequencing). We acknowledge that for large-scale screening studies with thousands of samples, traditional methods may be more appropriate. The justification emphasizes interpretability as the primary benefit rather than power alone, which aligns with Table 2 showing comparable power to traditional methods while providing unique interpretability advantages.

**Location in revised manuscript:** Limitations and future work section (line 443) and Table 6 (line 361)

---

### Prevalence Threshold

**Reviewer Comment:** The mention of a 10% prevalence threshold appears only in the computational efficiency section. Clarification is needed on whether prevalence filtering is an expected preprocessing step and how it plays a role along with the pre-filtering process.

**Response:** We have added clarification in the pre-filtering analysis section (line 381) explaining that "Prevalence filtering (retaining features present in ≥10% of samples) is an optional preprocessing step distinct from MeLSI's variance-based pre-filtering. When applied, prevalence filtering removes rare taxa before MeLSI analysis, while MeLSI's pre-filtering focuses on variance-based feature selection after preprocessing." This clarifies the distinction between the two filtering steps and their roles.

**Location in revised manuscript:** Line 381

---

### Tables

**Reviewer Comment:** Tables would benefit from clearer annotations and definitions of abbreviations.

**Response:** We have added clear annotations and abbreviation definitions to all tables (Tables 1-6). Each table now includes a concise footnote explaining abbreviations such as: n (sample size), p (number of taxa/features), F (PERMANOVA F-statistic), Power (empirical statistical power), Time (computation time in seconds), Rank (MeLSI's rank among 6 methods), and other table-specific terms. Footnotes are standardized for consistency, with repeated abbreviations (e.g., Rank) using the same concise format across tables.

**Location in revised manuscript:** Table 1: Line 262; Table 2: Line 284; Table 3: Line 310; Table 4: Line 336; Table 5: Line 353; Table 6: Line 369

---

### Notation

**Reviewer Comment:** Introduce notation more systematically in the "Metric learning: an emerging paradigm" section.

**Response:** We have systematically introduced notation in the "Metric learning: an emerging paradigm" section (lines 41-43). The revised text now formally defines: (1) the feature abundance matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$ with $n$ samples and $p$ taxa, (2) group labels $\mathbf{y} = (y_1, \ldots, y_n)$, (3) the positive semi-definite metric matrix $\mathbf{M} \in \mathbb{R}^{p \times p}$, (4) the Mahalanobis distance formula $d_M(\mathbf{x}_i, \mathbf{x}_j) = \sqrt{(\mathbf{x}_i - \mathbf{x}_j)^T \mathbf{M} (\mathbf{x}_i - \mathbf{x}_j)}$, and (5) the special case of diagonal $\mathbf{M}$ reducing to weighted Euclidean distance with feature-specific weights $M_{jj}$. This systematic introduction provides a clear mathematical foundation before the notation is used throughout the Methods section.

**Location in revised manuscript:** Introduction section, "Metric learning: an emerging paradigm" subsection (line 43)

---

### CLR Transformation

**Reviewer Comment:** The role and impact of the CLR transformation in conjunction with the proposed method should be discussed more explicitly. CLR transformation for MeLSI is mentioned only in the software availability section as part of the recommended usage.

**Response:** We have added explicit discussion of the CLR transformation role and impact in the Methods section. The revised text explains that: (1) CLR transformation converts relative abundances to log-ratios, making the data suitable for Euclidean distance while preserving relative relationships between taxa, (2) CLR treats abundance ratios more equitably than count-based metrics, which can be dominated by highly abundant taxa, (3) CLR may attenuate large fold-change signals compared to count-based metrics, as evidenced by results showing traditional count-based methods achieve higher F-statistics on synthetic data with large effects (3× fold change), and (4) CLR is particularly appropriate when signals are distributed across multiple taxa rather than concentrated in highly abundant taxa, and when interpretability through feature weights is prioritized. We also added explicit discussion in the Results section stating when the CLR-based approach is most appropriate versus when count-based methods may be preferable.

**Location in revised manuscript:** Methods section (Real data preprocessing) and Results section (Statistical power analysis)

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** While we can appreciate the novelty of this work, several critical technical concerns need to be addressed for improvement or clarification.

**Response:** We appreciate the reviewer's recognition of the work's novelty and thank them for the detailed technical feedback. We address each concern below.

---

### 1. Computational Cost vs. Benefit Trade-off

**Reviewer Comment:** MeLSI looks slower (minutes to hours) than traditional PERMANOVA (seconds). The paper argues this is acceptable given the gain in interpretability, but for large-scale or screening studies, this may be prohibitive. No clear power-time trade-off analysis is provided to justify when MeLSI is worth the extra computation.

**Response:** We have added an explicit power-time trade-off analysis in the Computational Performance section (line 431). The analysis demonstrates that MeLSI provides substantial value: pre-filtering increases statistical power by 36-37% while reducing computation time by 16-40% (Table 6). For typical microbiome studies (n=50-200, p=100-500), MeLSI completes in 2-30 minutes (Table 3), representing a modest time investment that yields both improved power (particularly for medium effect sizes, Table 2) and interpretability through feature weights. The power-time trade-off is most favorable when: (1) sample sizes are moderate (n=50-200), (2) interpretability is prioritized, and (3) pre-filtering is applied. For very large studies (n>500) or when only rapid screening is needed, traditional methods may be preferable.

**Location in revised manuscript:** Computational Performance section (line 431), Tables 2 (line 275), 3 (line 295), and 6 (line 361)

---

### 2. Limited Real-Data Performance Gains

**Reviewer Comment:** On Atlas1006, MeLSI's F-statistic was only 9.1% higher than Euclidean distance (5.141 vs. 4.711)-a modest improvement. On DietSwap, MeLSI reached significance (p=0.015) while Bray-Curtis was marginal (p=0.058), but the effect size difference is small. The authors acknowledge that on synthetic data with large effects, traditional metrics (Bray-Curtis, UniFrac) often outperformed MeLSI.

**Response:** We have revised the Conclusions section, Summary subsection (line 439) to explicitly position MeLSI's value proposition and acknowledge when traditional methods are preferable. The revised text clearly states that MeLSI is recommended when: (1) effect sizes are moderate (2× fold change) rather than very large, (2) interpretability through feature weights is needed to identify biologically relevant taxa, (3) traditional methods yield marginal results (p-values near 0.05), and (4) signals are distributed across multiple taxa rather than concentrated in highly abundant taxa. Traditional methods (Bray-Curtis, UniFrac) are preferable for: (1) large, obvious effects (3× fold change) where any method succeeds, (2) large-scale screening studies where speed is critical, and (3) when only omnibus significance testing is needed without feature-level interpretation. This positioning acknowledges that MeLSI's primary value is interpretability rather than marginal power gains, while clearly stating when each approach is most appropriate.

**Location in revised manuscript:** Conclusions section, Summary subsection (line 439) and Table 2 (line 275)

---

### 3. CLR Transformation May Limit Sensitivity

**Reviewer Comment:** MeLSI uses CLR-transformed data for Euclidean distance, which may attenuate large fold-change signals compared to count-based metrics (Bray-Curtis, UniFrac). This could explain why MeLSI underperforms on synthetic data with strong effects.

**Response:** We have added explicit discussion of CLR trade-offs in both the Methods section (line 234) and Results section (line 305). The revised text acknowledges that CLR transformation may attenuate large fold-change signals compared to count-based metrics, as evidenced by our results showing that traditional count-based methods achieve higher F-statistics on synthetic data with large effects (3× fold change). We explicitly state when the CLR-based approach is most appropriate: (1) when signals are distributed across multiple taxa rather than concentrated in highly abundant taxa, (2) when interpretability through feature weights is prioritized, and (3) when effect sizes are moderate rather than very large. For large, obvious effects (3× fold change), count-based methods (Bray-Curtis, UniFrac) may be preferable due to their sensitivity to abundance dominance. This positioning acknowledges the trade-off while clarifying when each approach is most appropriate.

**Location in revised manuscript:** Methods section (line 234) and Results section (line 305)

---

### 4. Lack of Covariate Adjustment

**Reviewer Comment:** MeLSI currently only supports simple group comparisons. No ability to adjust for confounders (age, BMI, antibiotics) or model continuous outcomes. This limits its utility in observational studies where confounding is common.

**Response:** We acknowledge this limitation and have added explicit discussion in the Limitations and future work section (line 447). The current implementation focuses on simple group comparisons, which represent the primary use case for microbiome beta diversity analysis. However, we recognize that observational studies often require adjustment for confounders (age, BMI, medication use, etc.) or modeling of continuous outcomes.

Covariate adjustment is a high-priority future extension. The underlying PERMANOVA framework (via vegan's `adonis2`) already supports covariates through formula syntax (e.g., `adonis2(dist ~ group + age + BMI)`), which provides a natural pathway for MeLSI extension. The technical approach would involve learning the metric on residuals after regressing out covariates, or incorporating covariates directly into the metric learning objective function. This extension would enable MeLSI to integrate with epidemiological frameworks where confounding is common.

For the current resubmission, we focus on independent group comparisons, which represent the majority of microbiome studies and allow for rigorous validation of the core method. The validation datasets (Atlas1006, DietSwap, SKIOME) demonstrate MeLSI's utility for this primary use case. Covariate adjustment support would require separate validation specific to that design type, which represents important follow-up work.

**Location in revised manuscript:** Limitations and future work section (line 445)

---

### 5. Sparse and Compositional Data Challenges

**Reviewer Comment:** While pre-filtering is used, the method does not explicitly model compositionality or zero-inflation beyond CLR. Aitchison geometry or Dirichlet-multinomial models might be more appropriate for compositional data.

**Response:** We address compositionality through CLR transformation, which converts compositional data to log-ratios. The Aitchison distance between two compositions $x$ and $y$ is defined as the Euclidean distance in CLR space: $d_A(x,y) = \sqrt{\sum(\text{clr}(x) - \text{clr}(y))^2}$. Our method computes Mahalanobis distance on CLR-transformed data, which is mathematically equivalent to a generalized Aitchison distance: $d_M(x,y) = \sqrt{(\text{clr}(x) - \text{clr}(y))^T \mathbf{M}^{-1} (\text{clr}(x) - \text{clr}(y))}$, where $\mathbf{M}$ is a learned positive-definite metric matrix. When $\mathbf{M} = \mathbf{I}$ (identity matrix), this reduces exactly to Aitchison distance. When $\mathbf{M} \neq \mathbf{I}$ (learned from data), this is a weighted Aitchison distance that adaptively weights dimensions based on their contribution to group separation. Thus, MeLSI uses Aitchison geometry as the foundation but extends it with an adaptive metric that learns feature-specific weights rather than treating all taxa equally—maintaining compositional properties while allowing the distance metric to adapt to dataset-specific signal structure. Zero-inflation is handled through pseudocounts (adding 1 before log transformation) rather than explicit Dirichlet-multinomial modeling. This approach works well in practice, as demonstrated by proper Type I error control (Table 1) and performance on real data (Atlas1006, DietSwap, SKIOME). We acknowledge that explicit zero-inflation models represent a potential future enhancement, though the current CLR + pseudocount approach maintains statistical validity through permutation testing.

**Location in revised manuscript:** Methods section (line 234) and Limitations and future work section (line 447)

---

### 6. Validation on Limited Real Datasets

**Reviewer Comment:** Only two real datasets (Atlas1006, DietSwap) were used for benchmarking, both from gut microbiome studies. Performance in other body sites (oral, skin) or disease cohorts is unknown.

**Response:** We have expanded real-world validation to include three published datasets representing diverse body sites and study designs: (1) **Atlas1006** (gut microbiome, 1,114 samples, sex-associated differences), (2) **DietSwap** (gut microbiome, 43 samples, dietary intervention), and (3) **SKIOME** (skin microbiome, 511 samples, three-group comparison: Atopic Dermatitis, Healthy, Psoriasis). The addition of SKIOME demonstrates MeLSI's utility across diverse body sites, validating performance beyond gut microbiome studies.

On the SKIOME dataset, MeLSI's omnibus test detected significant differences (F = 4.895, p = 0.005), comparable to Euclidean distance (F = 4.897, p = 0.001) but lower than count-based methods (Bray-Curtis: F = 16.275, Jaccard: F = 11.058, both p = 0.001). All pairwise comparisons remained significant after FDR correction (p = 0.005 for all pairs). While count-based methods achieved higher F-statistics on this dataset, MeLSI provides unique interpretability through learned feature weights that identify which taxa drive group separation—a capability that fixed metrics cannot supply. The SKIOME validation demonstrates MeLSI's generalizability across body sites (skin vs. gut) and validates multi-group capability beyond two-group comparisons.

We acknowledge that validation in additional body sites (oral, vaginal) or disease cohorts would further strengthen generalizability. However, the current validation across three datasets representing different body sites (gut, skin), study designs (observational, intervention, multi-group), and sample sizes (n=43 to n=1,114) provides evidence of MeLSI's broad applicability. The statistical framework (permutation testing, Type I error control) is consistent across all datasets, ensuring that the method's statistical properties are maintained regardless of body site or study design.

**Location in revised manuscript:** Real data validation section (lines 373-427), including SKIOME dataset subsection (lines 416-427) with Figure 3 (lines 420-427)

---

### 7. Parameter Sensitivity and Defaults

**Reviewer Comment:** Sensitivity analysis shows robustness, but defaults (e.g., B=30, m_frac=0.8) are not rigorously justified. For some datasets, different hyperparameters might yield better performance.

**Response:** We have expanded the Parameter Sensitivity section (lines 314-338) to explicitly justify the default parameters based on Table 4 results. The revised text explains that: (1) B=30 provides F-statistics (mean = 1.530, SD = 0.123) comparable to larger ensembles (B=50-100) while maintaining reasonable computation time (mean = 576.8s), (2) m_frac=0.8 balances performance (mean F = 1.530) with diversity among weak learners, as lower values (m_frac=0.5) show slightly higher F-statistics but reduced diversity, while higher values (m_frac=0.9-1.0) show slightly lower F-statistics, and (3) the robustness demonstrated across wide parameter ranges (B=10-100, m_frac=0.5-1.0) indicates that default parameters provide good performance across diverse datasets, though users may optimize for specific datasets if needed. This provides clear justification for the defaults while acknowledging that optimization may be beneficial for specific datasets.

**Location in revised manuscript:** Parameter Sensitivity section (line 314) and Table 4 (line 318)

---

### 8. Comparison to Other ML Methods

**Reviewer Comment:** No comparison to other interpretable ML models (e.g., Random Forest feature importance, logistic regression with regularization) that also provide taxa rankings. MeLSI is presented as unique, but similar insights might be obtained from simpler models. If authors do so, the ground truth should be well defined for benchmarking.

**Response:** While interpretable ML models like Random Forest (RF) and penalized regression provide feature rankings, MeLSI serves a distinct and complementary role in microbiome research:

**Metric-Based Inference:** Methods like RF identify features that maximize classification accuracy. In contrast, MeLSI learns a distance metric that optimizes group separation in beta-diversity space. This allows researchers to perform traditional statistical hypothesis testing (PERMANOVA) while simultaneously gaining the interpretability usually reserved for classification models.

**Global vs. Local Structure:** RF rankings are optimized for predicting group labels, whereas MeLSI feature weights describe the specific taxa driving the Aitchison geometry of the entire community. MeLSI captures global community-level differences (beta diversity), while RF focuses on local, per-feature contributions to classification.

**Type I Error Control:** Standard ML importance scores lack a formal framework for $p$-value generation in community composition testing. MeLSI bridges this gap by providing interpretable weights within a permutation-based inference framework that maintains proper Type I error control.

We have updated the Introduction (line 47) to clarify that while RF and logistic regression are excellent for prediction, MeLSI is specifically tailored for metric-based statistical inference, where the goal is to define and test community-level differences. The appropriate comparisons for MeLSI are other beta diversity methods used with PERMANOVA (Bray-Curtis, Euclidean, Jaccard, UniFrac), which we comprehensively evaluate.

**Location in revised manuscript:** Introduction (line 47) - distinction between prediction and inference

---

### 9. Scalability to Very High Dimensions

**Reviewer Comment:** The method scales as O(p²) in the number of features. For shotgun metagenomics (thousands of species or genes), this could be computationally prohibitive even with pre-filtering.

**Response:** We have expanded the Scalability section (line 312) to explicitly discuss scalability limits and acknowledge constraints. The revised text explains that MeLSI's O(p²) scaling becomes computationally prohibitive for very high-dimensional datasets (p>1000), with Table 3 demonstrating that computation time increases from 227.9s at p=50 to 8405.6s at p=1000. However, pre-filtering (retaining 70% of features) substantially mitigates this scaling, reducing effective dimensionality. For shotgun metagenomics with thousands of features, we recommend: (1) applying pre-filtering to reduce dimensionality, (2) considering feature aggregation (e.g., species-level rather than gene-level), or (3) using traditional methods if interpretability is not prioritized. The current implementation is most suitable for typical 16S rRNA datasets (p<1000) and metagenomic datasets with moderate dimensionality after preprocessing. This acknowledges the scalability constraint while providing practical guidance for high-dimensional applications.

**Location in revised manuscript:** Scalability section, Dimensionality scaling subsection (line 312) and Table 3 (line 295)

---

### 10. No Longitudinal/Paired Design Support

**Reviewer Comment:** MeLSI does not currently support repeated measures or paired samples, which are common in microbiome intervention studies.

**Response:** PERMANOVA (via vegan's adonis2) supports paired and longitudinal designs through the `strata` argument, which restricts permutations within blocks (e.g., within pairs or within subjects). MeLSI uses the same underlying PERMANOVA framework and can be extended to support paired designs using the same strata mechanism. The current implementation uses unrestricted permutations for independent group comparisons, which represent the majority of microbiome studies. Adding paired design support would require implementing block permutations (permuting labels within pairs/blocks rather than across all samples) and passing the strata argument to adonis2—a straightforward implementation enhancement that leverages PERMANOVA's existing capabilities. The current validations (Type I error control, power analysis) remain valid for independent group comparisons, as they were specifically designed to test that use case. Paired design support would require separate validation specific to that design type. For the resubmission, we focus on independent group comparisons, which represent the primary use case for microbiome beta diversity analysis, with paired design support as a future enhancement.

**Location in revised manuscript:** Methods section (line 206) - note on permutation strategies

---
