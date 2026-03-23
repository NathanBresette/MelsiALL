# MeLSI Project — CLAUDE.md

## What This Project Is

MeLSI (Metric Learning for Statistical Inference) is an R package for microbiome beta-diversity analysis. It learns a weighted distance metric from data and uses it with PERMANOVA for hypothesis testing. The manuscript is under review at **mSystems** (ASM Journals).

**Current status:** Round 3 resubmission in progress (Manuscript ID: mSystems00130-26)

---

## THE MANUSCRIPT TO EDIT

> **Always edit:** `manuscript/MeLSI_Research_Paper_mSystems.md`
>
> **Never edit:**
> - `resubmission_materials/new_manuscript/MeLSI_Research_Paper_mSystems_REVISED.md` — round 2 snapshot, for reference only
> - `resubmission_materials/old_manuscript/MeLSI_Research_Paper_mSystems_ORIGINAL.md` — round 1 snapshot, for reference only

The `.docx` and `.pdf` files in `manuscript/` are exports of the markdown. When revisions are done in the `.md`, they need to be re-exported to `.docx` for journal submission.

---

## Key People

- **Nathan Bresette** — first author, primary contact with journal
- **Aaron C. Ericsson** — co-author, provided feedback on manuscript
- **Carter Woods** — co-author
- **Ai-Ling Lin** — corresponding author (Dr.)

---

## File Structure

```
MeLSI/
├── CLAUDE.md                          ← this file
│
├── manuscript/                        ← ALL MANUSCRIPT WORK HAPPENS HERE
│   ├── MeLSI_Research_Paper_mSystems.md       ← CANONICAL MANUSCRIPT (edit this)
│   ├── MeLSI_Research_Paper_mSystems.docx     ← Word export for submission
│   ├── MeLSI_Manuscript_Resubmission.pdf      ← most recent submitted PDF
│   ├── MeLSI_Research_Paper_mSystems 3.pdf    ← older PDF version
│   ├── MeLSI_Research_Paper_mSystems_wordcount.md
│   ├── word_count.txt
│   ├── SUPPLEMENTARY_TABLES.md                ← supplementary tables (markdown)
│   ├── SUPPLEMENTARY_TABLES.docx              ← supplementary tables (Word)
│   ├── Cover_Letter.md
│   ├── Cover_Letter.pdf
│   ├── count_words_clean.R
│   ├── remove_tables_for_wordcount.R
│   └── figures/                               ← manuscript figures (TIF/PNG)
│       ├── atlas1006_vip_combined.tif         ← Figure 1 (submission-ready)
│       ├── atlas1006_vip_combined.png
│       ├── atlas1006_vip_combined_compressed.tif
│       ├── atlas1006_pcoa.tif
│       ├── atlas1006_pcoa_compressed.tif
│       ├── atlas1006_pcoa.png
│       ├── atlas1006_vip.png
│       ├── dietswap_combined.png              ← Figure 2 combined
│       ├── dietswap_vip_combined.png
│       ├── dietswap_vip.png
│       ├── dietswap_pcoa.png
│       ├── skiome_combined.png                ← Figure 3 combined
│       ├── skiome_combined.ps
│       ├── skiome_omnibus_pcoa_plot.pdf/.ps
│       ├── skiome_omnibus_vip_plot.pdf/.ps
│       ├── skiome_pcoa_plot.png
│       └── skiome_vip_plot.png
│
├── github/                            ← R package source code
│   ├── R/
│   │   ├── melsi_robust.R             ← main MeLSI function (primary algorithm)
│   │   ├── generate_test_data.R       ← test data generation
│   │   └── zzz.R                      ← package startup
│   ├── man/                           ← R documentation files (.Rd)
│   ├── tests/testthat/                ← package tests
│   ├── vignettes/melsi_tutorial.Rmd   ← package vignette
│   ├── inst/CITATION
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── NEWS.md
│   ├── LICENSE / LICENSE.md
│   └── README.md
│
├── extdata/                           ← precomputed results CSVs
│   ├── complete_type1_error.csv
│   ├── complete_power_analysis.csv
│   ├── complete_main_comparison.csv
│   ├── complete_scalability.csv
│   ├── complete_parameter_sensitivity.csv
│   ├── complete_prefiltering_value.csv
│   ├── complete_biological_validation.csv
│   ├── method_comparison_results.csv
│   └── SRP214545_taxonomy_abundances_SSU_v5.0.tsv  ← SKIOME raw data
│
├── figures/                           ← working figure PNGs (top-level copies)
│   ├── atlas1006_vip_combined.png
│   ├── atlas1006_vip_no_directionality.png
│   ├── atlas1006_vip.png
│   └── atlas1006_pcoa.png
│
├── reproducibility_scripts/           ← scripts to regenerate all tables/figures
│   ├── README.md
│   ├── figure_atlas1006_vip_pcoa.R    ← Figure 1 generation
│   └── hellbender/                    ← HPC (Hellbender cluster) job scripts
│       ├── HELLBENDER_GUIDE.md        ← how to run jobs on the cluster
│       ├── table1_type1_error.R / _job.sh
│       ├── table2_power_analysis.R / _job.sh
│       ├── table3_scalability.R / _job.sh
│       ├── table4_parameter_sensitivity.R / _job.sh
│       ├── table5_prefiltering.R / _job.sh        ← NOTE: script named table5 but this is
│       ├── table6_feature_correlation.R / _job.sh ←   manuscript Table 6 (see table note below)
│       ├── skiome_validation.R / _job.sh
│       ├── figure_dietswap_vip_pcoa.R / _job.sh
│       ├── combine_table[1-6]_results.R           ← merge HPC output CSVs
│       ├── combine_skiome_figures.R
│       ├── create_skiome_combined_figure.R
│       ├── figures/                               ← figures generated on cluster
│       │   ├── dietswap_vip_combined.pdf/.ps
│       │   ├── dietswap_pcoa.pdf/.ps
│       │   └── dietswap_vip_no_directionality.pdf/.ps
│       └── *.csv / *.RData                        ← simulation results and backups
│
├── resubmission_materials/            ← Round 2 archive (mSystems01710-25) — DO NOT EDIT
│   ├── Response_to_Reviewers_FORMAL.md / .docx   ← round 2 responses (reference)
│   ├── Cover_Letter.md / .docx
│   ├── REVISION_CHECKLIST.md
│   ├── README.md
│   ├── new_manuscript/
│   │   └── MeLSI_Research_Paper_mSystems_REVISED.md  ← round 2 snapshot — DO NOT EDIT
│   ├── old_manuscript/
│   │   └── MeLSI_Research_Paper_mSystems_ORIGINAL.md ← round 1 snapshot — DO NOT EDIT
│   └── figures_for_submission/        ← TIFs submitted in round 2
│       ├── Figure1_atlas1006_vip_combined.tif
│       ├── Figure2_dietswap_combined.tif
│       └── Figure3_skiome_combined.tif
│
├── resubmission_materials_round3/     ← Round 3 resubmission (mSystems00130-26) — ACTIVE
│   └── Response_to_Reviewers.md      ← point-by-point responses (in progress)
│
├── MeLSI_Manuscript_Resubmission.docx ← root-level copy of submission docx
├── MeLSI_0.99.2.tar.gz               ← built R package tarball
├── mSystems_FINAL_FORMATTING.md
├── mSystems_FINAL_READINESS_REPORT.md
├── mSystems_READINESS_CHECK.md
├── mSystems_formatting_progress.md
└── mSystems_formatting_summary.md
```

---

## Manuscript Tables (what they contain)

> **Note on script naming:** The hellbender scripts `table5_prefiltering.R` and `table6_feature_correlation.R` have names that do NOT match the manuscript table numbers. In the manuscript, Table 5 = Feature Correlation and Table 6 = Pre-filtering.

| Table | Manuscript # | Description |
|-------|-------------|-------------|
| Table 1 | 1 | Type I error control — empirical rejection rates across 600 simulations |
| Table 2 | 2 | Statistical power — empirical power (%) for ALL 6 methods across 450 simulations (restructured round 3) |
| Table S2 | S2 | Per-method power breakdown (supplementary) |
| Table 3 | 3 | Scalability — computation time vs. sample size and feature count |
| Table 4 | 4 | Parameter sensitivity — ensemble size (B) and feature fraction (m_frac) |
| Table S4 | S4 | Standard deviations for Table 4 (supplementary) |
| Table 5 | 5 | Feature correlation robustness — 200 simulations across r = 0, 0.3, 0.6, 0.8 |
| Table 6 | 6 | Pre-filtering — power and time savings across dimensionalities |
| Table S1 | S1 | Signal taxa recovery — Precision@k, Recall@k, Mean Rank, AUC-ROC |

---

## Real Datasets Used

| Dataset | n | Groups | Body site | Figure |
|---------|---|--------|-----------|--------|
| Atlas1006 | 1,114 | Male vs. Female | Gut | Figure 1 |
| DietSwap | 74 (baseline) | African vs. American diet | Gut | Figure 2 |
| SKIOME | 511 | Atopic Dermatitis / Healthy / Psoriasis | Skin | Figure 3 |

---

## Methods Compared Against MeLSI

1. Bray-Curtis
2. Jaccard
3. Unweighted UniFrac
4. Weighted UniFrac
5. Euclidean (CLR)
6. **MeLSI** (proposed)

All comparisons use PERMANOVA (vegan `adonis2`, 999 permutations).

---

## How We Write Responses to Reviewers

- Mirror the format in `resubmission_materials/Response_to_Reviewers_FORMAL.md` (round 2)
- Each reviewer comment gets its own `---`-delimited section
- Reviewer comment is quoted **word for word** under `**Reviewer Comment:**`
- Response goes under `**Response:**`
- Location (manuscript line numbers) goes under `**Location:**`
- Open items use `[TODO]`
- Start the document with a Cover Note and Summary of Major Revisions (fill in last)
- End with Closing Remarks summarizing total simulations run

---

## How We Handle Figures

- Working PNGs: `figures/` (top-level) and `manuscript/figures/`
- Submission-ready TIFs: `manuscript/figures/` and `resubmission_materials/figures_for_submission/`
- Figures generated on the cluster live in `reproducibility_scripts/hellbender/figures/`
- Final submission files must be high-resolution TIFF or EPS for ASM
- Multi-panel figures must be assembled into a single file
- All PCoA plots use **95% confidence ellipses** (changed from 68% in round 2)

---

## HPC / Cluster Notes

- Cluster: Hellbender (SLURM)
- See `reproducibility_scripts/hellbender/HELLBENDER_GUIDE.md` for job submission details
- Each table has a corresponding `_job.sh` submission script and `combine_*.R` merge script
- Results CSVs are stored in `reproducibility_scripts/hellbender/` and mirrored to `extdata/`

---

## Journal Requirements (mSystems / ASM)

- Response to Reviewers: separate file, NOT in cover letter
- Marked-up manuscript: upload as "Marked-Up Manuscript" (no figures)
- Clean manuscript: `.DOC/.DOCX`, remove previous version
- Figures: separate high-res files (TIFF or EPS preferred), multi-panel assembled into one file
- Supplemental material: combine into one file (preferred), max 10 files, legends separate from main manuscript
- Cover letter must include previous manuscript number and editor name (Gail Rosen)
- Submission link: https://msystems.msubmit.net/cgi-bin/main.plex

---

## Round History

| Round | Manuscript ID | Outcome |
|-------|--------------|---------|
| 1 | (original submission) | Reject with invitation to resubmit |
| 2 | mSystems01710-25 | Reject with invitation to resubmit |
| 3 | mSystems00130-26 | In progress |

---

## Active Round 3 Issues — STATUS

See `resubmission_materials_round3/Response_to_Reviewers.md` for the complete point-by-point responses.

1. **Table 2 / Table S2 inconsistency** — FIXED. Root cause: S2 had wrong traditional method F-stats. Table 2 restructured to show all methods' power. S2 corrected from CSV data.
2. **All methods same power** — ADDRESSED. Explained: power converges at large effects (expected), Jaccard/UWF have ~0% power (binary metrics can't detect fold changes).
3. **F-statistic ranking** — FIXED. Rank columns removed from Tables 2, 3, and 5. F-stats from different distance metrics are not comparable — let raw metrics speak for themselves.
4. **"9.1% improvement"** — FIXED. Removed throughout manuscript and conclusions.
5. **Terminology** — FIXED. "Pre-filtering score" ($I_j$) vs "learned metric weight" ($M_{jj}$) distinguished in Methods.
6. **Figures 1-3 redesigned** — IN PROGRESS. Script `figure_all_datasets_round3.R` needs rewrite: must load existing MeLSI results instead of re-running (stochastic = different F-stats each run). SKIOME has saved results (`skiome_omnibus_results_backup.RData`, `skiome_omnibus_distance_matrix.csv`). Atlas1006/DietSwap have no saved RData on Hellbender — script must re-run MeLSI for those but save results afterward. All traditional method distances (Euclidean, Bray-Curtis) are deterministic and fast. SKIOME uses all 3 groups (omnibus). Script must use `set.seed(42)` and save RData so results are reproducible.
7. **Correlation simulation details** — FIXED. Block structure (10 blocks × 20 taxa) described.
8. **Signal taxa recovery discussion** — FIXED. New subsection with benchmarks (random baseline) added.
9. **Weighted distance mechanism** — FIXED. New paragraph in Conclusions explaining signal-to-noise improvement.
10. **Real data F/p values** — VERIFIED. All values confirmed from Hellbender CSVs or original manuscript runs. No estimates remain.
11. **SKIOME interpretation** — FIXED. Manuscript now frames MeLSI vs Euclidean as the apples-to-apples comparison (both CLR-space), acknowledges BC yields larger pseudo-F but on a different scale, emphasizes interpretability as MeLSI's key advantage.

## Verified Real Data F/p Values

All values confirmed from Hellbender runs or original manuscript submissions. No estimates.

| Dataset | Method | F | p | Source |
|---------|--------|---|---|--------|
| Atlas1006 | MeLSI | 5.141 | 0.005 | Round 1 manuscript |
| Atlas1006 | Euclidean | 4.711 | 0.001 | Round 1 manuscript |
| Atlas1006 | Bray-Curtis | 4.442 | 0.001 | Round 1 manuscript |
| Atlas1006 | Jaccard | 1.791 | 0.144 | Round 1 manuscript |
| DietSwap | MeLSI | 2.856 | 0.015 | Hellbender confirmed |
| DietSwap | Bray-Curtis | 2.153 | 0.058 | Round 1 manuscript |
| DietSwap | Jaccard | 1.921 | 0.100 | Round 1 manuscript |
| DietSwap | Euclidean | 1.645 | 0.090 | Round 1 manuscript |
| SKIOME | MeLSI | 4.895 | 0.005 | Hellbender CSV |
| SKIOME | Euclidean | 4.897 | 0.001 | Hellbender CSV |
| SKIOME | Bray-Curtis | 16.275 | 0.001 | Hellbender CSV |
| SKIOME | Jaccard | 11.058 | 0.001 | Hellbender CSV |

Note: Atlas1006/DietSwap UniFrac not evaluated (no phylogenetic tree in public dataset objects).
SKIOME UniFrac not evaluated (no phylogenetic tree).
