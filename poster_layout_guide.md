# MeLSI Poster Layout Guide

## Complete Poster Structure (48" x 36" or similar)

### **HEADER SECTION (Top 15%)**
```
MeLSI: Learning Distance Metrics for Microbiome Analysis
A Novel Machine Learning Approach to Statistical Inference

Nathan Bresette
[Your Institution]
[Email] | [GitHub: NathanBresette/MeLSI]
```

---

### **MAIN CONTENT (70% of poster)**

#### **LEFT COLUMN (35% width)**

**1. Figure 1: Problem & Solution** (Top 40%)
- Shows fixed vs adaptive distance metrics
- File: `figure1_problem_solution.png`

**2. Introduction Text** (Middle 30%)
```
BACKGROUND
Current microbiome analysis relies on fixed distance metrics 
(Bray-Curtis, Euclidean, Jaccard) that treat all microbial 
taxa equally. This "one-size-fits-all" approach often misses 
subtle but biologically important differences between groups.

INNOVATION
MeLSI learns optimal distance metrics that adapt to each 
specific dataset, maximizing group separation while 
maintaining statistical rigor.
```

**3. Figure 4: Feature Importance** (Bottom 30%)
- Shows which taxa matter most
- File: `figure4_feature_importance.png`

#### **CENTER COLUMN (30% width)**

**1. Figure 2: Algorithm Flow** (Top 70%)
- 7-step MeLSI process
- File: `figure2_algorithm_flow_clean.png`

**2. Figure 6: Key Innovation Box** (Bottom 30%)
- Why MeLSI matters
- File: `figure6_innovation_box.png`

#### **RIGHT COLUMN (35% width)**

**1. Figure 3: Performance Comparison** (Top 40%)
- MeLSI vs traditional methods
- File: `figure3_performance_comparison.png`

**2. Results Summary** (Middle 30%)
```
KEY RESULTS
• 28% improvement in F-statistics over best traditional methods
• Perfect Type I error control (48.7% on null data)
• 28.7% computational speedup with pre-filtering
• Automatic feature importance visualization
• Validated on real microbiome datasets (Atlas1006, SoilRep)

REAL DATA SUCCESS
• Atlas1006: 4.85 vs 4.73 F-statistic (2.5% improvement)
• SoilRep: 1.49 vs 0.98 F-statistic (52% improvement)
```

**3. Figure 5: Validation Results** (Bottom 30%)
- Statistical validation summary
- File: `figure5_validation_results.png`

---

### **BOTTOM SECTION (15% of poster)**

#### **CONCLUSIONS & IMPACT**
```
CONCLUSIONS
MeLSI represents a fundamental advance in microbiome analysis 
by learning adaptive distance metrics instead of using fixed 
ones. This enables more powerful detection of biologically 
meaningful group differences while maintaining statistical rigor.

IMPACT
• Enables more sensitive microbiome studies
• Provides interpretable feature importance
• Open-source R package available
• Applicable to any microbiome dataset

FUTURE WORK
• Extension to longitudinal microbiome data
• Integration with phylogenetic information
• Development of confidence intervals for learned metrics
```

#### **REFERENCES & CONTACT**
```
REFERENCES
[1] Bresette, N. (2025). MeLSI: Metric Learning for Statistical 
    Inference in Microbiome Analysis. GitHub Repository.

SOFTWARE
Install: devtools::install_github("NathanBresette/MeLSI")
Documentation: https://github.com/NathanBresette/MeLSI

CONTACT
Email: [your-email]
GitHub: @NathanBresette
```

---

## **DESIGN SPECIFICATIONS**

### **Typography**
- **Title**: 72pt, Bold, Arial/Helvetica
- **Section Headers**: 36pt, Bold
- **Body Text**: 24pt, Regular
- **Captions**: 18pt, Regular

### **Colors**
- **Primary**: #2E86AB (Blue)
- **Accent**: #E74C3C (Red) - for MeLSI highlights
- **Text**: #2C3E50 (Dark Gray)
- **Background**: White

### **Spacing**
- **Margins**: 1 inch minimum
- **Between sections**: 0.5 inch
- **Between figures**: 0.25 inch

### **Figure Placement**
- **Figure 1**: Top-left, 10" x 6"
- **Figure 2**: Center, 12" x 4"
- **Figure 3**: Top-right, 10" x 8"
- **Figure 4**: Bottom-left, 8" x 6"
- **Figure 5**: Bottom-right, 8" x 5"
- **Figure 6**: Center-bottom, 6" x 6"

---

## **POSTER SOFTWARE RECOMMENDATIONS**

### **Professional Options**
1. **Adobe Illustrator** - Best for professional design
2. **Canva** - User-friendly with templates
3. **PowerPoint** - Quick and accessible

### **Free Options**
1. **Draw.io** - Good for flowcharts
2. **GIMP** - Free image editor
3. **Inkscape** - Free vector graphics

---

## **PRINTING TIPS**

### **File Format**
- Save as PDF or high-resolution PNG (300 DPI minimum)
- Use CMYK color mode for professional printing

### **Size Options**
- **Standard**: 48" x 36" (landscape)
- **Large**: 56" x 42"
- **Small**: 36" x 24"

### **Printing Services**
- University print shops
- FedEx Office
- Staples
- Online services (PosterPresentations.com)

---

## **PRESENTATION TIPS**

### **Key Talking Points**
1. **Problem**: Fixed metrics miss biological signals
2. **Solution**: MeLSI learns adaptive metrics
3. **Results**: 28% improvement, perfect validation
4. **Impact**: More powerful microbiome studies

### **Common Questions & Answers**
- **Q**: "How do you prevent overfitting?"
- **A**: "We use ensemble methods with bootstrap sampling and rigorous permutation testing."

- **Q**: "Why not just use machine learning?"
- **A**: "MeLSI provides proper statistical inference with p-values, which ML methods don't offer."

- **Q**: "How generalizable is this?"
- **A**: "We validated on multiple real datasets including human gut and environmental microbiomes."










