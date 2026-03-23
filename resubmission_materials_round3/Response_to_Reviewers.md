# Response to Reviewers
## MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis
**Manuscript ID:** mSystems00130-26

**Authors:** Nathan Bresette, Aaron C. Ericsson, Carter Woods, Ai-Ling Lin

---

## Cover Note

We thank the editor and reviewers for their continued engagement with this manuscript. We have carefully addressed all concerns raised in this round of review. Below we provide point-by-point responses with references to the revised manuscript.

### Summary of Major Revisions

- **Table 2 restructured:** Now reports empirical power (%) for all six methods side-by-side, replacing the previous MeLSI-only format with F-statistic ranking. This directly addresses the editor's and both reviewers' concerns about transparency and ranking methodology.
- **Table S2 corrected:** All traditional method F-statistics and power values have been verified against simulation output files and corrected. The previous version contained transcription errors in the traditional method columns at medium and large effect sizes.
- **F-statistic ranking removed:** Cross-metric F-statistic ranking columns have been removed from Tables 2, 3, and 5. Table 2 now reports empirical power for all six methods directly, letting the data speak for itself. The manuscript acknowledges that F-statistics from different distance metrics are not directly comparable because each defines a distinct geometric space.
- **"9.1% improvement" language removed:** All percentage-based F-statistic comparisons have been replaced with direct reporting of F-statistics and p-values throughout the manuscript and conclusions.
- **Terminology clarified:** Pre-filtering score ($I_j$) and learned metric weight ($M_{jj}$) are now explicitly distinguished in the Methods section, with consistent terminology throughout.
- **Signal taxa recovery discussion expanded:** New subsection provides detailed interpretation of recovery metrics across all simulated conditions with context benchmarks.
- **Correlation simulation details added:** Block correlation structure (10 blocks of 20 taxa, uniform within-block correlation, between-block independence) now fully described.
- **Figures 1-3 redesigned:** Each figure now includes four panels: (A) VIP plot with directionality, (B) MeLSI PCoA, (C) Euclidean PCoA, (D) Bray-Curtis PCoA, enabling direct visual comparison across methods. Figure 1 redundant left panel removed.
- **Weighted distance mechanism discussed:** New paragraph in Conclusions explains how learned weights improve signal-to-noise ratio in the distance calculation.
- **Real data F-statistics added:** Atlas1006 and DietSwap sections now report F-statistics and p-values for all six methods, not just MeLSI.

---

## RESPONSE TO EDITOR

**Editor Comment:** There seems to be some major issues with Table 2 and S2. I think the authors should significantly justify why all methods reach the same power, etc. And the authors should address all remaining concerns.

**Response:** We identified and corrected two issues. First, Table S2 contained transcription errors in the traditional method F-statistics at medium and large effect sizes (e.g., Bray-Curtis F was listed as 2.200 at Large/n=200, but the actual simulation value is 8.236). All values have been verified against the simulation output CSVs and corrected. Second, Table 2 has been completely restructured to show empirical power for all six methods side-by-side, replacing the previous MeLSI-only format that used F-statistic rankings. With correct data, the power convergence pattern is explained by two factors: (1) at large effect sizes with sufficient sample sizes, all abundance-sensitive methods (MeLSI, Euclidean, Bray-Curtis, Weighted UniFrac) converge to 100% power — this is expected statistical behavior, not an artifact; (2) Jaccard and Unweighted UniFrac show near-zero power because these binary (presence/absence) metrics cannot detect fold-change effects in abundance when taxa are already present in both groups. The F-statistic ranking column has been removed entirely, as reviewers correctly noted that F-statistics from different distance metrics are not directly comparable.

**Location:** Table 2 (restructured), Table S2 (corrected), Results section

---

## RESPONSE TO REVIEWER #1

### Comment 1 — Table 2 / Table S2 Ranking Inconsistency

**Reviewer Comment:** The MeLSI rank reported in Table 2 does not appear to align with the results shown in Table S2. For example, under the large effect size with n = 200 scenario, MeLSI has the largest mean F-statistic (3.129), which would correspond to rank 1/6, rather than 3/6 as reported in Table 2. The authors should verify the ranking calculations and ensure consistency between the tables.

**Response:** We thank the reviewer for identifying this inconsistency. The root cause was transcription errors in Table S2: the traditional method F-statistics at medium and large effect sizes were incorrect. For example, at Large/n=200, the correct values from our simulation output are Bray-Curtis F=8.236 and Weighted UniFrac F=8.659, both substantially higher than MeLSI's F=3.129. The previous Table S2 incorrectly listed these as F=2.200, which made MeLSI appear to have the highest F-statistic. With the corrected values, MeLSI's rank of 3/6 by F-statistic is consistent. However, we have taken the reviewer's broader point (Comment 3) and removed F-statistic rankings entirely, restructuring Table 2 to report empirical power for all six methods directly.

**Location:** Table 2 (restructured to show all methods' power), Table S2 (corrected with verified simulation data)

---

### Comment 2 — Text Claims vs. Table S2 Results (Large Effect F-Statistics)

**Reviewer Comment:** In lines 362-365, 530-531, and 534-535, the manuscript states that traditional count-based methods achieve higher F-statistics than MeLSI in synthetic datasets with large effects (3× fold change). However, the simulation results reported in Table S2 do not appear to fully support this statement. In the large-effect settings, MeLSI seems to demonstrate comparable performance/power relative to other methods. The authors should clarify this discrepancy and revise the text if necessary.

**Response:** This discrepancy arose because Table S2 contained incorrect traditional method F-statistics. With the corrected data, traditional count-based methods (Bray-Curtis, Weighted UniFrac) do indeed achieve substantially higher F-statistics than MeLSI at medium and large effects (e.g., BC F=8.236 vs MeLSI F=3.129 at Large/n=200), confirming the manuscript's original claim. However, we have revised the relevant text to focus on empirical power rather than F-statistics, since F-statistics from different distance metrics are not directly comparable (see Comment 3). The revised text states that count-based methods show higher empirical power at medium effects (Table 2), and that power converges across abundance-sensitive methods at large effects.

**Location:** Methods (CLR discussion), Results (Table 2 discussion), Conclusions (Summary)

---

### Comment 3 — Ranking Methods by F-Statistic vs. Empirical Power

**Reviewer Comment:** The manuscript reports empirical power as the primary evaluation metric across several validation studies. However, the methods are compared and ranked based on the mean PERMANOVA F-statistic rather than empirical power. Because the PERMANOVA F-statistic is not a direct measure of power and may be influenced by differences in dispersion or preprocessing steps, ranking methods using the F-statistic may not accurately reflect their detection ability. Furthermore, because different distance quantify dissimilarities in fundamentally different ways, the resulting distance structures, and therefore the corresponding pseudo-F statistics, may not be directly comparable across metrics. It would therefore be preferable to rank methods based on empirical power, which more directly reflects the probability of detecting true differences. Alternatively, if the authors choose to retain the F-statistic as a ranking criterion, the manuscript should clearly justify this choice and discuss its limitations.

**Response:** We agree with this assessment. F-statistics from different distance metrics are not directly comparable because each metric defines a distinct geometric space with different scale properties. We have made three changes: (1) Table 2 now reports empirical power for all six methods, allowing direct comparison of detection rates without cross-metric F-statistic ranking. (2) The F-statistic ranking column has been removed from Tables 2, 3, and 5 entirely. (3) The manuscript now includes an explicit statement that "PERMANOVA F-statistics computed from different distance metrics are not directly comparable because each metric defines a distinct geometric space, so cross-metric F-statistic comparisons should be interpreted cautiously." F-statistics are retained in Supplementary Tables for within-method assessment across conditions, with a note about their limitations for cross-metric comparison.

**Location:** Tables 2, 3, and 5 (rank columns removed), Results section (new caveat statement), Table S2 (added limitation note)

---

### Comment 4 — Interpretability / Signal Taxa Recovery Results Need More Discussion

**Reviewer Comment:** The manuscript briefly states that "the learned feature weights reliably identify true signal taxa," (lines 362-363) and the corresponding table report several recovery metrics. Because interpretability through identifying signal taxa is presented as a key advantage of the proposed method, I encourage the authors to elaborate more on these results in the manuscript. In particular, the manuscript should discuss how recovery performance varies across the simulated conditions (small, medium, and large effect sizes; varying sample sizes) and what these results imply about the reliability of the learned feature weights in realistic microbiome settings, where effect sizes are often modest. Providing a more detailed interpretation of these results would strengthen the claim that the method yields interpretable feature weights and can help identify biologically meaningful taxa. It may be helpful to provide context or benchmarks for interpreting these metrics (e.g., whether the reported precision/recall values are considered strong).

**Response:** We have added a new "Signal taxa recovery" subsection with detailed interpretation. At small effects (1.5× fold change), recovery is modest (Precision@5 = 0.104-0.148, AUC-ROC = 0.641-0.673), reflecting the inherent difficulty of identifying specific signal taxa when effects are subtle. At medium effects, performance improves substantially with sample size (Precision@5 from 0.356 at n=50 to 0.660 at n=200, AUC-ROC from 0.733 to 0.842). At large effects, recovery is strong (Precision@5 = 0.876-1.000, Mean Rank = 14.4-26.2 out of 200 taxa, AUC-ROC = 0.858-0.960). For context, random assignment would yield Precision@5 of 0.025-0.100 and AUC-ROC of 0.500, so even the modest small-effect recovery represents meaningful signal detection above chance. We conclude that MeLSI's learned weights are most reliable as a feature ranking tool when effects are moderate to large — precisely the conditions where identifying the specific taxa driving group separation is most biologically actionable.

**Location:** New "Signal taxa recovery" subsection in Results, between Table 2 discussion and Scalability analysis

---

### Comment 5 — Correlation Simulation Setup Lacks Detail

**Reviewer Comment:** The manuscript now evaluates robustness to feature correlation by varying correlation levels. However, the simulation setup does not provide sufficient detail to understand how correlation was introduced. For example, it is unclear whether correlations were imposed among all taxa, only among signal taxa, or within specific blocks of taxa. Providing more details about the data-generating mechanism and correlation structure would improve the reproducibility and interpretation of these simulations.

**Response:** We have added detailed description of the correlation structure. The 200 taxa were divided into 10 blocks of 20 taxa each, with uniform pairwise correlation $r$ imposed among taxa within each block and independence maintained between blocks. This block structure mimics realistic microbiome correlation patterns where functionally related or co-occurring taxa exhibit positive correlation while taxonomically distant groups remain independent. Within each block, correlated multivariate normal noise was generated using Cholesky decomposition and added to log-transformed abundances, scaled to preserve the original signal structure. Signal taxa were randomly distributed across blocks, reflecting the realistic scenario where correlated taxa may or may not include differentially abundant species.

**Location:** Feature correlation robustness section, Table 5 description

---

### Comment 6 — "9.1% Statistical Improvement" Framing Is Misleading

**Reviewer Comment:** In lines 455-459, the manuscript reports that MeLSI achieves a PERMANOVA F-statistic of 5.141 compared with 4.711 for Euclidean distance and interprets this as a "9.1% statistical improvement." Because the PERMANOVA F-statistic is not an effect size/a direct measure of power and depends on the underlying distance structure, expressing this difference as a percentage improvement may be misleading. It would be more appropriate to simply report the observed F-statistics and corresponding permutation p-values. More generally, statements throughout the manuscript that interpret mean F-statistics as evidence of improved power (e.g., lines 509-510) should be reconsidered or revised.

**Response:** We agree. All percentage-based F-statistic comparisons have been removed. The Atlas1006 section now simply reports F-statistics and p-values for all six methods without computing percentage differences: "MeLSI achieved F = 5.141 (p = 0.005) compared with Euclidean distance (F = 4.711, p = 0.001), Bray-Curtis (F = 3.789, p = 0.001)..." We have also added the caveat that "F-statistics from different distance metrics are not directly comparable because each defines a distinct geometric space." The Conclusions section has been revised to remove all percentage improvement claims, focusing instead on empirical power from Table 2 and the interpretability advantage.

**Location:** Atlas1006 results section, DietSwap results section, Computational performance section, Conclusions

---

## RESPONSE TO REVIEWER #2

### Comment 1 — All Methods Reaching the Same Power; F-Statistic Ranking Not Rigorous

**Reviewer Comment:** In the simulation studies, the authors compared the power across different methods. However, based on Table 2 and S2, all traditional methods reached the same power as the proposed MeLSI under all scenarios, which should not be a coincident. The authors should explore and interpret the reasons, and revisit the simulation procedures. Then, they made a rank based on F statistics. Considering the same power for all methods, using F statistics for ranking is not statistical rigorous to me. Please clarify.

**Response:** We have restructured Table 2 to show empirical power for all six methods, which reveals that the methods do NOT all reach the same power — the pattern is more nuanced. Three key findings emerge: (1) Bray-Curtis and Weighted UniFrac consistently achieve higher power than MeLSI at medium effects (e.g., BC 100% vs MeLSI 50% at Medium/n=100), reflecting their direct sensitivity to abundance fold changes in raw count data. (2) Jaccard and Unweighted UniFrac show near-zero power across all conditions because these binary (presence/absence) metrics cannot detect fold-change effects — multiplying the abundance of an already-present taxon does not change its presence/absence profile. (3) Power converges to 100% for all abundance-sensitive methods at large effects with sufficient sample size (n≥100), which is expected statistical behavior: sufficiently strong signals are detectable by any reasonable method. We have removed the F-statistic ranking entirely, as the reviewer correctly notes that ranking by F-statistics across different distance metrics is not statistically rigorous.

**Location:** Table 2 (restructured with all methods' power), Table 2 discussion text, Table S2 (corrected)

---

### Comment 2 — "Importance" vs. "Weight" Terminology Is Confusing

**Reviewer Comment:** In real data application, the authors provided more details on feature importance (line 464-482). However, it is hard to differentiate this "importance" term with the importance score introduced in the pre-filtering section (line 146-152). Then the authors used a term "weight" to indicate the importance in Figure 1-3, which make me more confused. It is suggested to clarify each terminology. Also, the content of line 464-482 should be defined and discussed earlier in the Methods section.

**Response:** We have clarified the terminology throughout the manuscript. Two distinct quantities are now explicitly defined: (1) The **pre-filtering score** ($I_j = |\mu_{1j} - \mu_{2j}| / \sqrt{\sigma_{1j}^2 + \sigma_{2j}^2}$), used solely for feature selection during the pre-filtering step. This is a ranking heuristic that identifies features with high between-group signal relative to within-group variation. (2) The **learned metric weight** ($M_{jj}$), which is the diagonal element of the metric matrix optimized during gradient-based metric learning. This represents each feature's contribution to the final distance calculation and is the basis for feature importance interpretation. The Methods section now states: "This pre-filtering score $I_j$ is used solely for feature selection and is distinct from the learned metric weights $M_{jj}$, which are optimized during metric learning and represent each feature's contribution to the final distance metric." The term "importance" has been replaced with "learned metric weight" where it referred to $M_{jj}$, and all figure axes now read "Learned Metric Weight." The explanation of how metric weights are computed and interpreted has been moved to the Methods (Problem Formulation) section, with application-specific discussion retained in the Results.

**Location:** Methods: Conservative pre-filtering section (terminology clarification added), Problem Formulation (wording updated), Figure 1-3 captions (updated axis labels), Feature importance subsections (terminology updated)

---

### Comment 3 — PCoA Plots for Traditional Methods Missing from Figures 2 & 3

**Reviewer Comment:** For Figure 2&3, the PCoA plot using other traditional distance metrics should be added and compared. The authors should discuss the different performance across methods in real data.

**Response:** Figures 2 and 3 now include PCoA ordination plots for Euclidean distance and Bray-Curtis dissimilarity alongside the MeLSI PCoA, enabling direct visual comparison. Each figure has four panels: (A) VIP plot with directionality, (B) MeLSI PCoA, (C) Euclidean PCoA, (D) Bray-Curtis PCoA.

For DietSwap, the MeLSI PCoA shows better visual separation between diet groups compared to both Euclidean and Bray-Curtis, consistent with MeLSI being the only method achieving significance (p=0.015 vs p≥0.058).

For SKIOME, Bray-Curtis yields a substantially larger pseudo-F statistic (F=16.275 vs MeLSI F=4.895), but as Reviewer 1 correctly noted, F-statistics from different distance metrics are not directly comparable because each defines a distinct geometric space. The most informative comparison is MeLSI (F=4.895) versus Euclidean (F=4.897), since both operate in CLR-transformed space: MeLSI performed essentially identically to the unweighted Euclidean baseline, indicating that the signal in the SKIOME dataset is broadly distributed across many taxa rather than concentrated in a small subset that MeLSI's weighting could amplify. Critically, while Bray-Curtis provides only an omnibus p-value, MeLSI additionally provides learned feature weights identifying which taxa drive group separation (e.g., Staphylococcus as the top-weighted taxon, consistent with its known role in atopic dermatitis), enabling biological interpretation within a single unified framework.

**Location:** Figure 2 (4-panel layout), Figure 3 (4-panel layout), DietSwap and SKIOME results text

---

### Comment 4 — Figure 1 Left Panel Is Redundant; Add PCoA Plot

**Reviewer Comment:** For Figure 1, the two panels are redundant. It it not necessary to present left panel. And the similar PCoA plot as Figure 2&3 should be looked at.

**Response:** The redundant left panel (feature weights without directionality coloring) has been removed. Figure 1 now has four panels: (A) VIP bar plot with directionality coloring (the previous right panel), (B) PCoA using MeLSI learned distance, (C) PCoA using Euclidean distance (CLR), (D) PCoA using Bray-Curtis dissimilarity. This enables direct visual comparison of ordination patterns across methods for the Atlas1006 dataset, consistent with the format of Figures 2 and 3.

**Location:** Figure 1 (redesigned 4-panel layout), Figure 1 caption

---

### Comment 5 — Insufficient Discussion of How Weighted Distance Improves Performance

**Reviewer Comment:** As stated by the authors, one of the major strength of MeLSI is learning the features' importance and using a kind of weighted distance metric. However, it is not very clear to me how does this "weighted distance" improve the statistical performance and biological interpretation of a beta-diversity analysis merely from their simulation and real data application. I suggest the authors include more in depth discussion.

**Response:** We have added a new paragraph in the Conclusions explaining the mechanism. The weighted distance works by assigning higher metric weights ($M_{jj}$) to taxa with strong between-group differences and lower weights to noise taxa. This effectively increases the signal-to-noise ratio in the distance calculation, concentrating the PERMANOVA analysis on the most informative features. The analogy is to how a domain expert might focus on key taxa — but learned directly from data rather than requiring prior knowledge. This weighting is particularly beneficial when a small subset of taxa drives group differences while many uninformative taxa add noise, as is typical in microbiome studies. The DietSwap results provide a concrete example: MeLSI achieved significance (p=0.015) while all fixed metrics remained marginal (p≥0.058), demonstrating that learned weights can elevate weak but real signals above the detection threshold. Additionally, the learned weights directly provide a ranked list of taxa driving separation, whereas fixed metrics provide only an omnibus p-value without feature-level attribution.

**Location:** Conclusions (new paragraph on weighted distance mechanism)

---

## RESPONSE TO REVIEWER #3

**Reviewer Comment:** No further comments.

**Response:** We thank Reviewer #3 for their assessment.

---

## Closing Remarks

We believe these revisions substantially strengthen the manuscript through: corrected and transparent simulation reporting (Table 2 showing all methods, Table S2 with verified data), removal of misleading cross-metric F-statistic comparisons, expanded discussion of signal taxa recovery and the weighted distance mechanism, detailed correlation simulation methodology, and redesigned multi-panel figures enabling direct visual comparison across methods. The total validation includes over 1,675 simulations across Type I error, power, scalability, parameter sensitivity, feature correlation, and pre-filtering analyses, plus three real-world datasets spanning two body sites (gut and skin). We remain available for any further discussion the editor or reviewers may find helpful.
