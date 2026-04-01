# Response to Reviewers
## MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis
**Manuscript ID:** mSystems00130-26

---

## RESPONSE TO EDITOR

**Editor Comment:** There seems to be some major issues with Table 2 and S2. I think the authors should significantly justify why all methods reach the same power, etc. And the authors should address all remaining concerns.

**Response:** We identified and corrected two issues. First, Table S2 contained transcription errors in the traditional method F-statistics at medium and large effect sizes (e.g., Bray-Curtis F was listed as 2.200 at Large/n=200; the correct simulation value is 8.236). All values have been verified against simulation output CSVs and corrected. Second, Table 2 has been restructured to show empirical power for all six methods side-by-side. With correct data, the power convergence pattern is explained by two factors: (1) at large effect sizes with sufficient sample sizes, all abundance-sensitive methods converge to 100% power, which is expected statistical behavior; (2) Jaccard and Unweighted UniFrac show near-zero power because binary metrics cannot detect fold-change effects when taxa are already present in both groups. The F-statistic ranking column has been removed, as F-statistics from different distance metrics are not directly comparable.

**Location:** Table 2 (line 358), Table S2 (Supplementary Tables), Results: Table 2 discussion (lines 364-382)

---

## RESPONSE TO REVIEWER #1

### Comment 1 — Table 2 / Table S2 Ranking Inconsistency

**Reviewer Comment:** The MeLSI rank reported in Table 2 does not appear to align with the results shown in Table S2. For example, under the large effect size with n = 200 scenario, MeLSI has the largest mean F-statistic (3.129), which would correspond to rank 1/6, rather than 3/6 as reported in Table 2. The authors should verify the ranking calculations and ensure consistency between the tables.

**Response:** The root cause was transcription errors in Table S2. At Large/n=200, the correct values are Bray-Curtis F=8.236 and Weighted UniFrac F=8.659, both substantially higher than MeLSI's F=3.129. The previous Table S2 incorrectly listed these as F=2.200, making MeLSI appear to have the highest F-statistic. With corrected values, MeLSI's rank of 3/6 is consistent. We have also taken the reviewer's broader point (Comment 3) and removed F-statistic rankings entirely, restructuring Table 2 to report empirical power for all six methods directly.

**Location:** Table 2 (line 358), Table S2 (Supplementary Tables)

---

### Comment 2 — Text Claims vs. Table S2 Results (Large Effect F-Statistics)

**Reviewer Comment:** In lines 362-365, 530-531, and 534-535, the manuscript states that traditional count-based methods achieve higher F-statistics than MeLSI in synthetic datasets with large effects (3× fold change). However, the simulation results reported in Table S2 do not appear to fully support this statement. In the large-effect settings, MeLSI seems to demonstrate comparable performance/power relative to other methods. The authors should clarify this discrepancy and revise the text if necessary.

**Response:** This discrepancy arose from the Table S2 transcription errors described in Comment 1. With corrected data, count-based methods do achieve substantially higher F-statistics than MeLSI at medium and large effects (e.g., BC F=8.236 vs MeLSI F=3.129 at Large/n=200). We have revised the text to focus on empirical power rather than F-statistics, since cross-metric F-statistic comparisons are not statistically rigorous (see Comment 3). The revised text states that count-based methods show higher empirical power at medium effects and that power converges across abundance-sensitive methods at large effects.

**Location:** Results: Table 2 discussion (lines 364-382), Conclusions: Summary (lines 562-563)

---

### Comment 3 — Ranking Methods by F-Statistic vs. Empirical Power

**Reviewer Comment:** The manuscript reports empirical power as the primary evaluation metric across several validation studies. However, the methods are compared and ranked based on the mean PERMANOVA F-statistic rather than empirical power. Because the PERMANOVA F-statistic is not a direct measure of power and may be influenced by differences in dispersion or preprocessing steps, ranking methods using the F-statistic may not accurately reflect their detection ability. Furthermore, because different distance quantify dissimilarities in fundamentally different ways, the resulting distance structures, and therefore the corresponding pseudo-F statistics, may not be directly comparable across metrics. It would therefore be preferable to rank methods based on empirical power, which more directly reflects the probability of detecting true differences. Alternatively, if the authors choose to retain the F-statistic as a ranking criterion, the manuscript should clearly justify this choice and discuss its limitations.

**Response:** The reviewer raises a valid point. F-statistics from different distance metrics are not directly comparable because each metric defines a distinct geometric space. We made three changes: (1) Table 2 now reports empirical power for all six methods directly. (2) F-statistic ranking columns have been removed from Tables 2, 3, and 5. (3) The manuscript now states that cross-metric F-statistic comparisons should be interpreted cautiously. F-statistics are retained in Supplementary Tables for within-method assessment across conditions.

**Location:** Table 2 (line 358), Results: Table 2 discussion (lines 364-382), Table S2 (Supplementary Tables)

---

### Comment 4 — Interpretability / Signal Taxa Recovery Results Need More Discussion

**Reviewer Comment:** The manuscript briefly states that "the learned feature weights reliably identify true signal taxa," (lines 362-363) and the corresponding table report several recovery metrics. Because interpretability through identifying signal taxa is presented as a key advantage of the proposed method, I encourage the authors to elaborate more on these results in the manuscript. In particular, the manuscript should discuss how recovery performance varies across the simulated conditions (small, medium, and large effect sizes; varying sample sizes) and what these results imply about the reliability of the learned feature weights in realistic microbiome settings, where effect sizes are often modest. Providing a more detailed interpretation of these results would strengthen the claim that the method yields interpretable feature weights and can help identify biologically meaningful taxa. It may be helpful to provide context or benchmarks for interpreting these metrics (e.g., whether the reported precision/recall values are considered strong).

**Response:** We have added a new "Signal taxa recovery" subsection. At small effects (1.5x fold change), recovery is modest (Precision@5 = 0.104-0.148, AUC-ROC = 0.641-0.673). At medium effects, performance improves with sample size (Precision@5 from 0.356 at n=50 to 0.660 at n=200, AUC-ROC from 0.733 to 0.842). At large effects, recovery is strong (Precision@5 = 0.876-1.000, AUC-ROC = 0.858-0.960). For context, random assignment yields Precision@5 of 0.025-0.100 and AUC-ROC of 0.500, so even modest small-effect recovery represents meaningful detection above chance. MeLSI's learned weights are most reliable at moderate to large effects, which are the conditions where identifying specific signal taxa is most biologically actionable.

**Location:** Signal taxa recovery subsection (lines 383-397)

---

### Comment 5 — Correlation Simulation Setup Lacks Detail

**Reviewer Comment:** The manuscript now evaluates robustness to feature correlation by varying correlation levels. However, the simulation setup does not provide sufficient detail to understand how correlation was introduced. For example, it is unclear whether correlations were imposed among all taxa, only among signal taxa, or within specific blocks of taxa. Providing more details about the data-generating mechanism and correlation structure would improve the reproducibility and interpretation of these simulations.

**Response:** We have added detailed description of the correlation structure. The 200 taxa were divided into 10 blocks of 20 taxa each, with uniform pairwise correlation $r$ imposed among taxa within each block and independence maintained between blocks. This block structure mimics realistic microbiome correlation patterns where functionally related or co-occurring taxa exhibit positive correlation while taxonomically distant groups remain independent. Within each block, correlated multivariate normal noise was generated using Cholesky decomposition and added to log-transformed abundances, scaled to preserve the original signal structure. Signal taxa were randomly distributed across blocks, reflecting the realistic scenario where correlated taxa may or may not include differentially abundant species.

**Location:** Feature correlation robustness section (line 426), Table 5 (line 444)

---

### Comment 6 — "9.1% Statistical Improvement" Framing Is Misleading

**Reviewer Comment:** In lines 455-459, the manuscript reports that MeLSI achieves a PERMANOVA F-statistic of 5.141 compared with 4.711 for Euclidean distance and interprets this as a "9.1% statistical improvement." Because the PERMANOVA F-statistic is not an effect size/a direct measure of power and depends on the underlying distance structure, expressing this difference as a percentage improvement may be misleading. It would be more appropriate to simply report the observed F-statistics and corresponding permutation p-values. More generally, statements throughout the manuscript that interpret mean F-statistics as evidence of improved power (e.g., lines 509-510) should be reconsidered or revised.

**Response:** This framing was imprecise and we have corrected it. All percentage-based F-statistic comparisons have been removed. Results are now reported as F-statistics and p-values directly: MeLSI (F = 4.841, p = 0.005), Euclidean (F = 4.711, p = 0.001), Bray-Curtis (F = 4.442, p = 0.001), Jaccard (F = 1.791, p = 0.144). We note that the MeLSI F-statistic differs slightly from the reviewer's cited value of 5.141 because Figures 1-3 were regenerated with a fixed random seed (set.seed(42)) for reproducibility; both values yield p = 0.005 and the same biological conclusion. The manuscript now includes an explicit caveat that cross-metric F-statistic comparisons are not directly comparable. The Conclusions have been revised to focus on empirical power and interpretability rather than F-statistic differences.

**Location:** Atlas1006 dataset (lines 487-510), DietSwap dataset (lines 510-521), Conclusions: Summary (lines 562-563)

---

## RESPONSE TO REVIEWER #2

### Comment 1 — All Methods Reaching the Same Power; F-Statistic Ranking Not Rigorous

**Reviewer Comment:** In the simulation studies, the authors compared the power across different methods. However, based on Table 2 and S2, all traditional methods reached the same power as the proposed MeLSI under all scenarios, which should not be a coincident. The authors should explore and interpret the reasons, and revisit the simulation procedures. Then, they made a rank based on F statistics. Considering the same power for all methods, using F statistics for ranking is not statistical rigorous to me. Please clarify.

**Response:** Table 2 has been restructured to show empirical power for all six methods, which reveals that methods do not all reach the same power. Three patterns emerge: (1) Bray-Curtis and Weighted UniFrac achieve higher power than MeLSI at medium effects (e.g., BC 100% vs MeLSI 50% at Medium/n=100), reflecting their direct sensitivity to abundance fold changes in raw count data. (2) Jaccard and Unweighted UniFrac show near-zero power because binary metrics cannot detect fold-change effects: multiplying the abundance of an already-present taxon does not change its presence/absence profile. (3) Power converges to 100% for all abundance-sensitive methods at large effects with sufficient sample size, which is expected statistical behavior. F-statistic ranking has been removed, as the reviewer correctly notes it is not statistically rigorous across different distance metrics.

**Location:** Table 2 (line 358), Table 2 discussion (lines 364-382), Table S2 (Supplementary Tables)

---

### Comment 2 — "Importance" vs. "Weight" Terminology Is Confusing

**Reviewer Comment:** In real data application, the authors provided more details on feature importance (line 464-482). However, it is hard to differentiate this "importance" term with the importance score introduced in the pre-filtering section (line 146-152). Then the authors used a term "weight" to indicate the importance in Figure 1-3, which make me more confused. It is suggested to clarify each terminology. Also, the content of line 464-482 should be defined and discussed earlier in the Methods section.

**Response:** Two distinct quantities are now explicitly defined in the Methods. (1) The **pre-filtering score** ($I_j = |\mu_{1j} - \mu_{2j}| / \sqrt{\sigma_{1j}^2 + \sigma_{2j}^2}$) is used solely for feature selection and ranks features by between-group signal relative to within-group variation. (2) The **learned metric weight** ($M_{jj}$) is the diagonal element of the metric matrix optimized during gradient-based metric learning and represents each feature's contribution to the final distance. The Methods now explicitly states that $I_j$ and $M_{jj}$ are distinct quantities. The term "importance" has been replaced with "learned metric weight" throughout, figure axes now read "Learned Metric Weight," and the explanation of metric weight computation has been moved to the Problem Formulation section of the Methods.

**Location:** Methods: Conservative pre-filtering (lines 149-150), Problem formulation (lines 141-143), Figure legends (lines 757, 765, 771), Results: Learned metric weights subsections (lines 495, 518, 541)

---

### Comment 3 — PCoA Plots for Traditional Methods Missing from Figures 2 & 3

**Reviewer Comment:** For Figure 2&3, the PCoA plot using other traditional distance metrics should be added and compared. The authors should discuss the different performance across methods in real data.

**Response:** Figures 2 and 3 now include PCoA ordination plots for Euclidean distance and Bray-Curtis dissimilarity alongside the MeLSI PCoA, enabling direct visual comparison. Each figure has four panels: (A) VIP plot with directionality, (B) MeLSI PCoA, (C) Euclidean PCoA, (D) Bray-Curtis PCoA.

For DietSwap, the MeLSI PCoA shows better visual separation between diet groups compared to both Euclidean and Bray-Curtis, consistent with MeLSI being the only method achieving significance (p=0.015 vs p≥0.066).

For SKIOME, the most informative comparison is MeLSI (F=4.972) versus Euclidean (F=4.897), since both operate in CLR-transformed space. Their near-identical performance is a meaningful biological finding: it indicates the signal is broadly distributed across many taxa rather than concentrated in a subset, and MeLSI correctly learns this, assigning relatively uniform weights rather than artificially inflating a few taxa. This reflects statistical integrity, not a limitation. Bray-Curtis (F=16.275) yields a larger pseudo-F but operates in a fundamentally different geometric space and is not directly comparable. Regardless of detection equivalence, MeLSI provides learned feature weights identifying which taxa drive group separation (e.g., Staphylococcus, consistent with its known role in atopic dermatitis), which no fixed-distance method can supply.

**Location:** DietSwap dataset (lines 510-521), SKIOME dataset (lines 522-541), Figure 2 legend (line 765), Figure 3 legend (line 771)

---

### Comment 4 — Figure 1 Left Panel Is Redundant; Add PCoA Plot

**Reviewer Comment:** For Figure 1, the two panels are redundant. It it not necessary to present left panel. And the similar PCoA plot as Figure 2&3 should be looked at.

**Response:** The redundant left panel (feature weights without directionality coloring) has been removed. Figure 1 now has four panels: (A) VIP bar plot with directionality coloring (the previous right panel), (B) PCoA using MeLSI learned distance, (C) PCoA using Euclidean distance (CLR), (D) PCoA using Bray-Curtis dissimilarity. This enables direct visual comparison of ordination patterns across methods for the Atlas1006 dataset, consistent with the format of Figures 2 and 3.

**Location:** Figure 1 legend (line 757)

---

### Comment 5 — Insufficient Discussion of How Weighted Distance Improves Performance

**Reviewer Comment:** As stated by the authors, one of the major strength of MeLSI is learning the features' importance and using a kind of weighted distance metric. However, it is not very clear to me how does this "weighted distance" improve the statistical performance and biological interpretation of a beta-diversity analysis merely from their simulation and real data application. I suggest the authors include more in depth discussion.

**Response:** We have added a new paragraph in the Conclusions. Standard PERMANOVA with fixed metrics treats all taxa equally, which is problematic in high-dimensional data: uninformative taxa add noise to the distance calculation and dilute signal from the few taxa that actually differ. MeLSI learns diagonal metric weights $M_{jj}$ that upweight informative taxa and downweight noise, concentrating statistical power on the most relevant features. This is why MeLSI detected significance on DietSwap (p=0.015) while Euclidean, Bray-Curtis, and Jaccard remained marginal (p>=0.066). The learned weights also provide a ranked list of taxa driving group separation, whereas fixed distance methods yield only an omnibus p-value. This unification of hypothesis testing and feature attribution within a single permutation framework is MeLSI's core contribution.

**Location:** Conclusions (lines 576-586)

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** No further comments.

**Response:** We thank Reviewer #3 for their assessment.

---

## Closing Remarks

We thank the editor and reviewers for their continued engagement with this work. We note that since the previous submission, MeLSI has been accepted into Bioconductor, providing independent validation of the software's quality and reproducibility standards.
