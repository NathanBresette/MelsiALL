# Complete MeLSI Poster Layout Plan with Colors

## **POSTER SPECIFICATIONS**
- **Size**: 48" x 36" (landscape)
- **Resolution**: 300 DPI
- **Background**: Pure White (#FFFFFF)

---

## **COLOR PALETTE**
```css
Background: #FFFFFF (Pure White)
Primary Blue: #2E86AB (Headers, algorithm boxes)
Accent Red: #E74C3C (MeLSI results, highlights)
Text Gray: #2C3E50 (Main text)
Light Gray: #95A5A6 (Traditional methods, subtle elements)
Success Green: #27AE60 (Validation checkmarks)
Medium Gray: #7F8C8D (Captions)
Light Blue: #EBF3FD (Subtle backgrounds)
```

---

## **HEADER SECTION (Top 15% - 7.2" height)**

### **Background**: White (#FFFFFF)
### **Content**:
```
MeLSI: Learning Distance Metrics for Microbiome Analysis
[Color: #2E86AB, Size: 72pt, Bold, Arial]

A Novel Machine Learning Approach to Statistical Inference
[Color: #2C3E50, Size: 36pt, Regular, Arial]

Nathan Bresette
[Color: #E74C3C, Size: 24pt, Bold, Arial]

[Your Institution] | [Email] | GitHub: NathanBresette/MeLSI
[Color: #7F8C8D, Size: 18pt, Regular, Arial]
```

---

## **MAIN CONTENT (70% - 33.6" height)**

### **LEFT COLUMN (35% width - 16.8")**

#### **1. Introduction Text (Top 40% - 13.4" height)**
- **Background**: White (#FFFFFF)
- **Border**: 2px solid #2E86AB
- **Padding**: 0.5" all sides
- **Content**:
```
BETA DIVERSITY & DISTANCE METRICS
[Color: #2E86AB, Size: 24pt, Bold]

Beta diversity measures differences between microbial communities 
using distance metrics like Bray-Curtis, Euclidean, and Jaccard. 
These metrics quantify how similar or different samples are based 
on their taxonomic composition.

[Color: #2C3E50, Size: 18pt, Regular]

THE PROBLEM
[Color: #2E86AB, Size: 24pt, Bold]

Traditional distance metrics are fixed and treat all microbial 
taxa equally. This "one-size-fits-all" approach often misses 
subtle but biologically important differences between groups.

[Color: #2C3E50, Size: 18pt, Regular]

THE SOLUTION
[Color: #2E86AB, Size: 24pt, Bold]

MeLSI learns optimal distance metrics that adapt to each 
specific dataset, maximizing group separation while 
maintaining statistical rigor.

[Color: #2C3E50, Size: 18pt, Regular]
```

#### **2. Figure 1: Problem & Solution (Middle 30% - 10.1" height)**
- **File**: `figure1_problem_solution.png`
- **Size**: 14" x 6"
- **Colors**: 
  - Traditional bars: #95A5A6
  - MeLSI bars: #E74C3C
  - Text: #2C3E50
- **Caption**: "Traditional methods treat all taxa equally, MeLSI learns optimal weights"

#### **3. Figure 4: Feature Importance (Bottom 30% - 10.1" height)**
- **File**: `figure4_feature_importance.png`
- **Size**: 12" x 7"
- **Colors**: 
  - Bars: #E74C3C
  - Text: #2C3E50
- **Caption**: "Which taxa matter most for group separation?"

---

### **CENTER COLUMN (30% width - 14.4")**

#### **1. Figure 2: Algorithm Flow (Top 70% - 23.5" height)**
- **File**: `figure2_algorithm_flow_clean.png`
- **Size**: 13" x 5"
- **Colors**: 
  - Boxes: #2E86AB
  - Arrows: #95A5A6
  - Text: White (on blue boxes)
- **Caption**: "Ensemble metric learning with rigorous statistical validation"

#### **2. Figure 6: Key Innovation Box (Bottom 30% - 10.1" height)**
- **File**: `figure6_innovation_box.png`
- **Size**: 10" x 8"
- **Colors**: 
  - Border: #2E86AB
  - Background: #EBF3FD (subtle)
  - Text: #2C3E50
- **Caption**: "Why MeLSI matters"

---

### **RIGHT COLUMN (35% width - 16.8")**

#### **1. Figure 3: Performance Comparison (Top 40% - 13.4" height)**
- **File**: `figure3_performance_comparison.png`
- **Size**: 14" x 9"
- **Colors**: 
  - MeLSI bars: #E74C3C
  - Traditional bars: #95A5A6
  - Significance markers: #27AE60
- **Caption**: "MeLSI outperforms traditional methods with higher F-statistics"

#### **2. Results Summary (Middle 30% - 10.1" height)**
- **Background**: White (#FFFFFF)
- **Border**: 2px solid #E74C3C
- **Padding**: 0.5" all sides
- **Content**:
```
KEY RESULTS
[Color: #E74C3C, Size: 24pt, Bold]

• 28% improvement in F-statistics over best traditional methods
• Perfect Type I error control (48.7% on null data)
• 28.7% computational speedup with pre-filtering
• Automatic feature importance visualization
• Validated on real microbiome datasets (Atlas1006, SoilRep)

[Color: #2C3E50, Size: 18pt, Regular]

REAL DATA SUCCESS
[Color: #E74C3C, Size: 24pt, Bold]

• Atlas1006: 4.85 vs 4.73 F-statistic (2.5% improvement)
• SoilRep: 1.49 vs 0.98 F-statistic (52% improvement)

[Color: #2C3E50, Size: 18pt, Regular]
```

#### **3. Figure 5: Validation Results (Bottom 30% - 10.1" height)**
- **File**: `figure5_validation_results.png`
- **Size**: 12" x 6"
- **Colors**: 
  - Success indicators: #27AE60
  - Background: #95A5A6 (low opacity)
- **Caption**: "Statistical validation across all test scenarios"

---

## **BOTTOM SECTION (15% - 7.2" height)**

### **Background**: White (#FFFFFF)

#### **Conclusions & Impact (Left 50%)**
- **Border**: 2px solid #2E86AB
- **Padding**: 0.5" all sides
- **Content**:
```
CONCLUSIONS
[Color: #2E86AB, Size: 24pt, Bold]

MeLSI represents a fundamental advance in microbiome analysis 
by learning adaptive distance metrics instead of using fixed 
ones. This enables more powerful detection of biologically 
meaningful group differences while maintaining statistical rigor.

[Color: #2C3E50, Size: 18pt, Regular]

IMPACT
[Color: #2E86AB, Size: 24pt, Bold]

• Enables more sensitive microbiome studies
• Provides interpretable feature importance
• Open-source R package available
• Applicable to any microbiome dataset

[Color: #2C3E50, Size: 18pt, Regular]
```

#### **Future Work & Contact (Right 50%)**
- **Border**: 2px solid #E74C3C
- **Padding**: 0.5" all sides
- **Content**:
```
FUTURE WORK
[Color: #E74C3C, Size: 24pt, Bold]

• Extension to longitudinal microbiome data
• Integration with phylogenetic information
• Development of confidence intervals for learned metrics

[Color: #2C3E50, Size: 18pt, Regular]

SOFTWARE & CONTACT
[Color: #E74C3C, Size: 24pt, Bold]

Install: devtools::install_github("NathanBresette/MeLSI")
Documentation: https://github.com/NathanBresette/MeLSI
Email: [your-email] | GitHub: @NathanBresette

[Color: #2C3E50, Size: 18pt, Regular]
```

---

## **DESIGN SPECIFICATIONS**

### **Typography**
- **Font Family**: Arial (clean, professional, widely available)
- **Title**: 72pt, Bold, #2E86AB
- **Section Headers**: 24pt, Bold, #2E86AB or #E74C3C
- **Body Text**: 18pt, Regular, #2C3E50
- **Captions**: 16pt, Regular, #7F8C8D

### **Spacing**
- **Margins**: 1" minimum
- **Between sections**: 0.5"
- **Between figures**: 0.25"
- **Text padding**: 0.5" inside boxes

### **Borders & Boxes**
- **Border width**: 2px
- **Border colors**: #2E86AB (blue) or #E74C3C (red)
- **Border radius**: 5px (slightly rounded corners)
- **Box shadows**: None (clean, flat design)

### **Figure Specifications**
- **Resolution**: 300 DPI minimum
- **Format**: PNG with transparent backgrounds
- **Consistent styling**: All figures use the same color palette
- **Captions**: Below each figure, 16pt, #7F8C8D

---

## **PRINTING CHECKLIST**

### **File Preparation**
- [ ] Save as PDF or high-resolution PNG
- [ ] Use CMYK color mode for professional printing
- [ ] Ensure 300 DPI resolution
- [ ] Test print a small section first

### **Color Verification**
- [ ] All colors are within the specified palette
- [ ] Sufficient contrast for readability
- [ ] No gradients or complex effects
- [ ] Colors look good in grayscale (accessibility)

### **Content Review**
- [ ] All figures are properly sized and positioned
- [ ] Text is readable at poster size
- [ ] No typos or formatting errors
- [ ] Contact information is current

This layout creates a professional, visually appealing poster that tells a clear story: Problem → Solution → Results → Impact, with consistent use of your MeLSI color scheme throughout.
