# MeLSI Resubmission Checklist

## Folder Structure Created ✓
- [x] `resubmission_materials/` folder created
- [x] `Response_to_Reviewers.md` template created

## Documents to Prepare (Per Editor Instructions)

### Required Files:
- [ ] **Response to Reviewers** - Point-by-point responses (NOT in cover letter)
- [ ] **Marked-Up Manuscript** - Compare copy (without figures) showing changes
- [ ] **Clean DOC/DOCX** - Revised manuscript in Word format
- [ ] **Figure Files** - Each as separate, editable, high-resolution file (TIFF or EPS preferred)
  - [ ] Figure 1
  - [ ] Figure 2 (needs fixing - variance inconsistency)
  - [ ] Figure 3
  - [ ] Figure 4
  - [ ] Figure 5
  - [ ] DietSwap VIP plot (NEW)
  - [ ] DietSwap PCoA plot (NEW)
- [ ] **Supplemental Material** - With legends separate (if applicable)

## Critical Issues to Address (Must Fix for Resubmission)

### Type I Error Control
- [ ] Run 100+ repeated simulations under null hypothesis
- [ ] Calculate empirical rejection rate at α=0.05
- [ ] Examine null p-value distribution
- [ ] Update Table 1 with proper statistical validation
- [ ] Add discussion of results

### Statistical Power Analysis
- [ ] Run multiple simulations per effect size condition
- [ ] Calculate power (proportion of significant results)
- [ ] Create power curves if possible
- [ ] Update Table 2 with power estimates
- [ ] Explore multiple sample sizes

### Sample Size Exploration
- [ ] Add systematic sample size variation to Type I error analysis
- [ ] Add systematic sample size variation to power analysis
- [ ] Report results across sample sizes (e.g., n=20, 50, 100, 200, 500)

### Figure 2 Issues
- [x] Fix PCoA1 variance inconsistency (18.4% vs 21.5%) - Updated to correct value: 21.5% for PCoA1, 10.9% for PCoA2
- [x] Revise "clear separation" language (compare to traditional methods)
- [x] Explain rationale for 68% confidence ellipses (Changed to 95% instead)
- [x] Add discussion of biological vs statistical significance

### DietSwap Results
- [ ] Create VIP plot for DietSwap data
- [ ] Create PCoA plot for DietSwap data
- [x] Improve transition in text
- [ ] Add complete presentation of DietSwap results

## Important Issues (Should Address)

### Double-Dipping/Overfitting
- [ ] Clarify that metric is relearned on each permutation
- [ ] Add explicit discussion of why this prevents double-dipping
- [ ] Discuss overfitting prevention explicitly
- [ ] Add evaluation/supporting evidence for overfitting prevention

### CLR Transformation
- [ ] Add explicit discussion of CLR role and impact in main methods/results
- [ ] Discuss trade-offs (sensitivity vs. compositionality handling)
- [ ] Explain when CLR approach is appropriate vs. count-based

### Pre-filtering Assumptions
- [ ] Discuss implications of using Gaussian-assumed statistics (t-test, ANOVA)
- [ ] Explain design choice for varying dimensionality across effect sizes
- [x] Clarify prevalence filtering role vs. pre-filtering

### Interpretability Validation
- [ ] Validate recovery of true signal taxa in synthetic data
- [ ] Test across varying effect sizes, sample sizes, dimensions
- [ ] Strengthen positioning around interpretability advantage

### Computational Justification
- [ ] Better justify computational time cost
- [ ] Add power-time trade-off analysis (if possible)
- [ ] Acknowledge limitations more clearly

## Nice-to-Have (Can Address or Acknowledge)

### Multi-Group Extensions
- [ ] Add multi-group synthetic validation OR
- [ ] Acknowledge as limitation/future work

### Feature Correlation
- [ ] Add simulations with correlated features OR
- [ ] Acknowledge as limitation/future work

### Additional Real Datasets
- [ ] Add more datasets OR
- [ ] Acknowledge limitation clearly

### Ensemble vs Single-Learner
- [ ] Add comparison OR
- [ ] Better justify ensemble approach in discussion

### Comparison to Other ML Methods
- [ ] Add comparison to Random Forest, regularized regression OR
- [ ] Justify why not included

## Minor Revisions

### Text Improvements
- [x] Add clearer table annotations and abbreviations
- [ ] Systematically introduce notation in "Metric learning" section
- [x] Clarify directionality/log2FC calculation details
- [ ] Better justify parameter defaults (B=30, m_frac=0.8)
- [x] Clarify "variance" in ensemble size analysis (Table 4)
- [ ] Clearly state sample sizes in pre-filtering analysis

### Limitations Discussion
- [ ] Add explicit discussion of covariate adjustment limitation
- [ ] Acknowledge longitudinal/paired design limitation
- [ ] Discuss scalability limits (O(p²))
- [ ] Better position when MeLSI is vs. isn't recommended

## Status Tracking

**Phase 1: Critical Statistical Validation** (MUST DO)
- [ ] Type I error expansion
- [ ] Power analysis expansion
- [ ] Sample size exploration

**Phase 2: Figure/Results Fixes** (MUST DO)
- [x] Figure 2 fixes (ellipses changed to 95%, language revised, biological significance added)
- [ ] DietSwap figures (transition improved, figures still needed)

**Phase 3: Method Clarifications** (SHOULD DO)
- [ ] Double-dipping clarification
- [ ] CLR discussion
- [ ] Pre-filtering assumptions

**Phase 4: Minor Revisions** (NICE TO HAVE)
- [ ] Text improvements
- [ ] Table improvements
- [ ] Additional validations (if time permits)

---

## Notes
- Remember: This is a NEW submission (not a revision)
- Include previous manuscript number in cover letter
- Mention editor name (Gail Rosen) in cover letter
