# mSystems Formatting Summary for MeLSI Paper

## Article Type
**Methods and Protocols** - This article describes a major technical and methodological development in systems microbiology (a new bioinformatic method for microbiome beta diversity analysis).

## Word Count
- **Body text (excluding references, tables, figure legends):** 4,470 words
- **Limit:** 5,000 words ✓
- **Abstract:** 245 words (limit: 250 words) ✓
- **Importance:** 149 words (limit: 150 words) ✓

## Major Formatting Changes

### 1. Two-Part Abstract
**Original:** Single abstract section
**mSystems:** Split into two sections:
- **ABSTRACT** (245 words): Describes what was done and found, focusing on methods and results
- **IMPORTANCE** (149 words): Non-technical explanation of why the work matters, written for broad readership

### 2. Title Page Elements Added
- Running title: "Metric Learning for Microbiome Analysis" (54 characters including spaces)
- Correspondent footnote placeholder
- Word count statements for abstract and body text
- Affiliation structure using superscript letters

### 3. Section Headings
**Changed to ALL CAPS for main sections** (per mSystems style):
- ABSTRACT
- IMPORTANCE
- INTRODUCTION
- MATERIALS AND METHODS
- RESULTS
- DISCUSSION
- ACKNOWLEDGMENTS
- REFERENCES

**Subsections:** Sentence case with bold formatting

### 4. Section Reorganization
**Materials and Methods** now appears BEFORE Results (can be after Introduction or after Discussion; we chose the former for logical flow).

### 5. Data Availability Statement
Added a dedicated paragraph at the end of Materials and Methods:
```
Data availability. All code for the MeLSI R package is available at 
https://github.com/NathanBresette/MeLSI. All validation experiments are 
fully reproducible using scripts available in the reproducibility_scripts 
directory. Real datasets (Atlas1006 and DietSwap) are publicly available 
via the microbiome R package and original publications.
```

### 6. References
**Changed to numbered citation-sequence system** (mSystems style):
- In-text citations now use parenthetical numbers: (1), (2, 3), etc.
- References numbered in order of first appearance
- Journal names abbreviated per PubMed standards
- Full author lists (no "et al." truncation)
- Format: Author. Year. Title. Journal volume:pages.

### 7. Content Condensation
To fit within the 5,000 word limit, we:
- **Combined Results and Discussion subsections** where appropriate
- **Removed redundant explanations** between Introduction and Methods
- **Streamlined validation experiment descriptions**
- **Condensed limitations discussion** while retaining key points
- **Integrated conclusions into Discussion** rather than separate section
- **Removed the separate "Conclusions" section** (Section 5 from original)

### 8. Structure Simplification
**Original sections:** 
1. Introduction (4 subsections)
2. Methods (7 subsections with multiple sub-subsections)
3. Results (7 subsections)
4. Discussion (7 subsections)
5. Conclusions

**mSystems format:**
- INTRODUCTION (no numbered subsections, just bold paragraph lead-ins)
- MATERIALS AND METHODS (bold paragraph lead-ins for each component)
- RESULTS (bold paragraph lead-ins)
- DISCUSSION (integrated conclusions, bold paragraph lead-ins)

### 9. Tables and Figures
**Note:** The current manuscript includes table numbers referenced (Tables 1-5) but the actual tables need to be formatted separately:
- Tables should be simple and formatted according to mSystems guidelines
- Figure references removed (need to be added as separate figure files)
- Figure legends should be in a separate section after References

### 10. Mathematical Notation
- Retained LaTeX math notation (mSystems accepts this)
- Ensured all symbols render correctly in PDF

## What Needs to Be Done Before Submission

### Required Additions:
1. **Complete affiliations** for all authors with institutional addresses
2. **Corresponding author email** address
3. **ORCID iDs** for all authors (recommended)
4. **Funding statement** in Acknowledgments
5. **Author contributions** using CRediT taxonomy (optional but recommended)
6. **Conflict of interest statement**

### Tables:
- Format Tables 1-5 as separate files or incorporate into manuscript with proper formatting
- Ensure table footnotes are included
- Verify that table captions are concise

### Figures (if including):
- Create separate high-resolution figure files (TIFF, EPS, or PDF format)
- Write detailed figure legends
- Ensure figures meet mSystems specifications (minimum 300 dpi)
- Consider whether Atlas1006 VIP and PCoA figures should be included as formal figures

### Supplemental Material (if any):
- Determine if any validation tables should be moved to supplemental material
- If including supplemental material, create separate file with own References section

### Final Checks:
1. Verify all citations are numbered correctly and match References
2. Confirm all references have complete information
3. Check that all abbreviations are defined on first use
4. Review Data Availability statement for completeness
5. Ensure GitHub repository is public and up-to-date
6. Add DOIs to key references (ASM will add most automatically during copyediting)

## Files Created
1. **MeLSI_Research_Paper_mSystems.md** - Markdown formatted for mSystems
2. **MeLSI_Research_Paper_mSystems.pdf** - PDF preview

## Comparison with Original
**Original file:** MeLSI_Research_Paper.md (8,120 total words)
**mSystems version:** MeLSI_Research_Paper_mSystems.md (4,470 body words + ~1,500 reference words)

## Next Steps
1. Review the PDF to ensure formatting is correct
2. Complete all author information and affiliations
3. Add funding statement and acknowledgments
4. Decide on table and figure presentation
5. Submit through mSystems Editorial Manager: https://msystems.asm.org/



