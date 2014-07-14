---
title       : Bioconductor RNA-Seq workflow
subtitle    : 
author      : Michael Love
job         : 
framework   : html5slides        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---

# Bioconductor RNA-Seq workflow

1. preparing gene models
2. read counting
3. QA (quality assessment) & EDA (exploratory data analysis)
4. differential expression analysis
5. annotating results

---

# Preparing gene models

[GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

---

# The generalized linear model part 1

\[ K_{ij} \sim \text{NB}( \mu_{ij}, \alpha_i )  \]

- Read count $K_{ij}$ for gene *i* sample *j*.
- 2 parameter count distribution: mean $\mu$, dispersion $\alpha$

\[ \mu_{ij} = s_{ij} q_{ij} \]

- Normalized by size factor $s_{ij}$, with $q_{ij}$ remaining
- Often size factor $s_j$

---

# The generalized linear model part 2

\[ \log q_{ij} = \sum_r x_{jr} \beta_{ir} \]

\[ \begin{array}{c}
\log q_1 \\
\log q_2 \\
\log q_3 \\
\log q_4 \end{array}
= \left( \begin{array}{cc}
1 & 0 \\
1 & 0 \\
1 & 1 \\
1 & 1 \end{array} \right)
\left( \begin{array}{c}
\beta_0 \\
\beta_1 \end{array} \right) \] 

---

# Multigroup comparisons in DESeq2

\[ \begin{array}{c}
\log q_1 \\
\log q_2 \\
\log q_3 \\
\log q_4 \\
\log q_5 \\
\log q_6 \end{array}
= \left( \begin{array}{cccc}
1 & 1 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 \\
1 & 0 & 0 & 1 \\
1 & 0 & 0 & 1 \\ \end{array} \right)
\left( \begin{array}{c}
\beta_0 \\
\beta_1 \\
\beta_2 \\
\beta_3 \end{array} \right) \] 

---

# ReportingTools

<center><img src="reptools.png" width=700/></center>

