# MeLSI Project — CLAUDE.md

## What This Project Is

MeLSI (Metric Learning for Statistical Inference) is an R package for microbiome beta-diversity analysis. It learns a weighted distance metric from data and uses it with PERMANOVA for hypothesis testing. The manuscript is under review at **mSystems** (ASM Journals).

**Current status:** Round 3 resubmission READY TO SUBMIT (Manuscript ID: mSystems00130-26)

**MeLSI R package accepted into Bioconductor** — install via `BiocManager::install("MeLSI")`

---

## THE MANUSCRIPT TO EDIT

> **Always edit:** `manuscript/MeLSI_Research_Paper_mSystems.md`
>
> **Never edit:**
> - `resubmission_materials/new_manuscript/MeLSI_Research_Paper_mSystems_REVISED.md` — round 2 snapshot, for reference only
> - `resubmission_materials/old_manuscript/MeLSI_Research_Paper_mSystems_ORIGINAL.md` — round 1 snapshot, for reference only

The `.docx` files in `manuscript/` are exports of the markdown. When revisions are done in the `.md`, re-export with:
```bash
cd manuscript && pandoc MeLSI_Research_Paper_mSystems.md -o MeLSI_Research_Paper_mSystems.docx
```
Then copy to `resubmissionv3/`.

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
│       ├── atlas1006_figure1.tif              ← Figure 1 — 4-panel, 4200×3000px ✓
│       ├── atlas1006_figure1.pdf
│       ├── dietswap_figure2.tif               ← Figure 2 — 4-panel, 4200×3000px ✓
│       ├── dietswap_figure2.pdf
│       ├── skiome_figure3.tif                 ← Figure 3 — 4-panel, 4200×3000px, y-axis fixed ✓
│       ├── skiome_figure3.pdf
│       └── (legacy individual panel files — for reference only)
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
├── resubmission_materials_round3/     ← Round 3 working files (mSystems00130-26)
│   ├── Response_to_Reviewers.md      ← point-by-point responses (COMPLETE)
│   └── Response_to_Reviewers.docx    ← Word export
│
├── resubmissionv3/                    ← UPLOAD THIS FOLDER TO JOURNAL ✓
│   ├── Cover_Letter.md / .docx
│   ├── MeLSI_Research_Paper_mSystems.docx
│   ├── Response_to_Reviewers.docx
│   ├── SUPPLEMENTARY_TABLES.docx
│   ├── Figure1_atlas1006.tif          ← 4200×3000px, 300dpi
│   ├── Figure2_dietswap.tif           ← 4200×3000px, 300dpi
│   └── Figure3_skiome.tif             ← 4200×3000px, 300dpi (fixed y-axis)
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

- Each reviewer comment gets its own `---`-delimited section
- Reviewer comment is quoted **word for word** under `**Reviewer Comment:**`
- Response goes under `**Response:**`
- Location (manuscript line numbers) goes under `**Location:**`
- Open items use `[TODO]`
- **No cover note, no summary of revisions, no author names** — journal says to put those in the cover letter only; response should jump straight into reviewer comments
- End with Closing Remarks summarizing total simulations run
- After editing `.md`, rebuild `.docx` with pandoc and copy to `resubmissionv3/`

---

## How We Handle Figures

- Submission-ready TIFs: `manuscript/figures/` (atlas1006_figure1, dietswap_figure2, skiome_figure3)
- All figures are 4-panel layout: (A) VIP plot, (B) MeLSI PCoA, (C) Euclidean PCoA, (D) Bray-Curtis PCoA
- All PCoA plots use **95% confidence ellipses**
- Figure scripts run on Hellbender (no X11 display — only PDF output works directly from R)
- TIF conversion from PDF: use ghostscript only: `gs -dNOPAUSE -dBATCH -sDEVICE=tiff24nc -r300 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sOutputFile=out.tif in.pdf`
- **Do NOT use `sips` for PDF→TIF** — it produces blurry output
- Target: 4200×3000px (14"×10" at 300dpi) matching Figure 1 and 2

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
| 3 | mSystems00130-26 | Ready to submit |

---

## Round 3 Issues — ALL RESOLVED ✓

1. **Table 2 / Table S2 inconsistency** — FIXED.
2. **All methods same power** — ADDRESSED.
3. **F-statistic ranking** — FIXED. Rank columns removed from Tables 2, 3, and 5.
4. **"9.1% improvement"** — FIXED. Removed throughout.
5. **Terminology** — FIXED. $I_j$ vs $M_{jj}$ distinguished in Methods.
6. **Figures 1-3 redesigned** — COMPLETE. 4-panel TIFs in `manuscript/figures/` and `resubmissionv3/`. SKIOME y-axis taxon names fixed (regenerated via `figure_skiome_round3.R` on Hellbender).
7. **Correlation simulation details** — FIXED.
8. **Signal taxa recovery discussion** — FIXED.
9. **Weighted distance mechanism** — FIXED.
10. **Real data F/p values** — VERIFIED.
11. **SKIOME interpretation** — FIXED.
12. **Bioconductor acceptance** — Added to Data Availability, Software Availability, cover letter, and Response to Reviewers closing remarks.

## Verified Real Data F/p Values

All values confirmed from reproducible Hellbender runs with set.seed(42). No estimates.

| Dataset | Method | F | p | Source |
|---------|--------|---|---|--------|
| Atlas1006 | MeLSI | 4.841 | 0.005 | Hellbender set.seed(42) |
| Atlas1006 | Euclidean | 4.711 | 0.001 | Hellbender set.seed(42) |
| Atlas1006 | Bray-Curtis | 4.442 | 0.001 | Hellbender set.seed(42) |
| Atlas1006 | Jaccard | 1.791 | 0.144 | Round 1 manuscript |
| DietSwap | MeLSI | 3.063 | 0.015 | Hellbender set.seed(42), DI vs HE groups |
| DietSwap | Bray-Curtis | 2.153 | 0.066 | Hellbender set.seed(42), DI vs HE groups |
| DietSwap | Jaccard | 1.921 | 0.100 | Round 1 manuscript |
| DietSwap | Euclidean | 1.670 | 0.077 | Hellbender set.seed(42), DI vs HE groups |
| SKIOME | MeLSI | 4.972 | 0.005 | skiome_omnibus_results.csv |
| SKIOME | Euclidean | 4.897 | 0.001 | Hellbender CSV |
| SKIOME | Bray-Curtis | 16.275 | 0.001 | Hellbender CSV |
| SKIOME | Jaccard | 11.058 | 0.001 | Hellbender CSV |

Note: Atlas1006/DietSwap UniFrac not evaluated (no phylogenetic tree in public dataset objects).
SKIOME UniFrac not evaluated (no phylogenetic tree).
DietSwap grouping: diet group (DI=high-fiber, HE=Western), NOT nationality. n=74 baseline samples.
