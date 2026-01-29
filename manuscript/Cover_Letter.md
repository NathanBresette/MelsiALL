---
geometry: margin=1in
fontsize: 11pt
pagestyle: empty
---

\vspace*{0.3in}

\begin{flushleft}
Editor-in-Chief\\
mSystems\\
RE: Resubmission of Manuscript ID mSystems01710-25\\
Title: "MeLSI: Metric Learning for Statistical Inference in Microbiome Community Composition Analysis"
\end{flushleft}

\vspace{0.2in}

Dear Editor,

On behalf of my co-authors at the University of Missouri, I am pleased to resubmit our revised manuscript for consideration in mSystems. We have leveraged this revision period to fundamentally strengthen the statistical foundation of our proposed framework, MeLSI, moving from case-study validation to a rigorous, simulation-based "gold standard" approach.

The primary focus of this revision was addressing the critical need for distributional proof of Type I error control and statistical power. We believe the following upgrades position MeLSI as a statistically robust tool for the microbiome research community:

\textbf{Statistical Validity and Power}

\textbf{Type I Error Control:} We expanded our validation to include 600 independent simulations across multiple sample sizes, demonstrating that MeLSI maintains empirical rejection rates between 3\% and 6\%, aligning perfectly with the nominal 5\% significance level.

\textbf{Comprehensive Power Analysis:} Through 450 new simulations across small, medium, and large effect sizes, we verified that MeLSI provides competitive or superior sensitivity compared to traditional fixed metrics, particularly when signals are distributed across multiple taxa.

\textbf{Robustness to Ecological Correlation:} Addressing a key concern for biological data, we have demonstrated that MeLSI maintains stable F-statistics and power even under high feature correlation ($r=0.8$), suggesting it is uniquely suited for the complex co-occurrence patterns inherent in microbial communities.

\textbf{Interpretability and Software Accessibility}

\textbf{Feature Recovery Metrics:} We now provide quantitative evidence—including Precision at k and AUC-ROC—of MeLSI's ability to accurately recover "true" signal taxa. On real-world datasets like Atlas1006 and DietSwap, MeLSI identified key biological drivers that traditional methods partially obscured.

\textbf{Bioconductor Integration:} Reflecting our commitment to high-quality software standards, the MeLSI R package is currently under review for inclusion in the Bioconductor repository, ensuring long-term stability and interoperability for the bioinformatics community.

By integrating the adaptive nature of metric learning with the rigor of permutation-based inference, MeLSI offers researchers a way to move beyond "one-size-fits-all" distance metrics without sacrificing statistical validity.

We thank the reviewers for their insights, which have significantly improved this work. We have no competing interests to disclose.

\vspace{0.2in}

Sincerely,

\vspace{0.1in}

Nathan Bresette\\
PhD Student, Institute for Data Science and Informatics\\
Roy Blunt NextGen Precision Health, University of Missouri
