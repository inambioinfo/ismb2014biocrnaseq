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

# Preparing gene models part 1

[GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)


```r
library( "GenomicFeatures" )
# takes ~10 min
txdb <- makeTranscriptDbFromBiomart(
  biomart="ensembl",
  dataset="hsapiens_gene_ensembl" )
```

smart to use `saveDb()` to only do this once

---

# Preparing gene models part 2

- `makeTranscriptDbFromGFF()` accepts GTF
- `library(TxDb.Hsapiens.UCSC.hg19.knownGene)` ready to go
- soon also, `AnnotationDbi` will offer ready to go

---

# Super useful


```r
seqlevelsStyle(gr) <- "NCBI"
seqlevelsStyle(gr) <- "UCSC"
```

---

# What's in txdb?

`sqlite3 txdb.sqlite '.tables'`

`cds exon metadata transcript chrominfo gene splicing`

---

# What's in txdb?


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb.Hsapiens.UCSC.hg19.knownGene
```

```
## TranscriptDb object:
## | Db type: TranscriptDb
## | Supporting package: GenomicFeatures
## | Data source: UCSC
## | Genome: hg19
## | Organism: Homo sapiens
## | UCSC Table: knownGene
## | Resource URL: http://genome.ucsc.edu/
## | Type of Gene ID: Entrez Gene ID
## | Full dataset: yes
## | miRBase build ID: GRCh37
## | transcript_nrow: 82960
## | exon_nrow: 289969
## | cds_nrow: 237533
## | Db created by: GenomicFeatures package from Bioconductor
## | Creation time: 2014-03-17 16:15:59 -0700 (Mon, 17 Mar 2014)
## | GenomicFeatures version at creation time: 1.15.11
## | RSQLite version at creation time: 0.11.4
## | DBSCHEMAVERSION: 1.0
```

---

# Extract exons for each gene


```r
# takes ~30 seconds
exonsByGene <- exonsBy( txdb, by="gene" )
```

---

# Extract exons for each gene


```r
# TODO show exonsByGene
```

---

# Read counting


```r
library( "Rsamtools" )
bamLst <- BamFileList( fls, yieldSize=2000000 )
```

---

# Read counting


```r
library( "GenomicAlignments" )
register( MulticoreParam( workers=4 ) )
se <- summarizeOverlaps( features=exonsByGene,
                        reads=bamLst,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
```

---

# SummarizedExperiment

<center><img src="summarizedexperiment.png" width=700/></center>

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

* results are logarithm base 2 fold changes 

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

---

# Support

- mailing list:
    1. package name in title
    2. describe experiment, what's the *biological question*
    3. provide code
    4. sessionInfo()
- `browseVignettes("pkg")`
- `?function`
