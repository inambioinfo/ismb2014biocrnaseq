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
3. EDA (exploratory data analysis)
4. differential expression analysis
5. annotating results

---

# Preparing gene models


```r
library( "GenomicFeatures" )
# takes ~10 min
txdb <- makeTranscriptDbFromBiomart( biomart="ensembl",
                                    dataset="hsapiens_gene_ensembl" )
```

smart to use `saveDb()` to only do this once

- [GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)

---

# Preparing gene models

- `makeTranscriptDbFromGFF()` accepts GTF
- `library(TxDb.Hsapiens.UCSC.hg19.knownGene)` ready to go
- soon also, `AnnotationHub` will offer ready to go

---

# Super useful


```r
seqlevelsStyle(gr) <- "NCBI"
seqlevelsStyle(gr) <- "UCSC"
```

---

# Extract exons for each gene




```r
# takes ~30 seconds
exonsByGene <- exonsBy( txdb, by="gene" )
```


```r
head( exonsByGene[[1]] )
```

```
## GRanges with 6 ranges and 2 metadata columns:
##       seqnames               ranges strand |   exon_id       exon_name
##          <Rle>            <IRanges>  <Rle> | <integer>     <character>
##   [1]        X [99883667, 99884983]      - |    664095 ENSE00001459322
##   [2]        X [99885756, 99885863]      - |    664096 ENSE00000868868
##   [3]        X [99887482, 99887565]      - |    664097 ENSE00000401072
##   [4]        X [99887538, 99887565]      - |    664098 ENSE00001849132
##   [5]        X [99888402, 99888536]      - |    664099 ENSE00003554016
##   [6]        X [99888402, 99888536]      - |    664100 ENSE00003658801
##   ---
##   seqlengths:
##                    1                 2 ...            LRG_99
##            249250621         243199373 ...             13294
```

---

# Read counting


```r
library( "Rsamtools" )
bamLst <- BamFileList( fls, yieldSize=2000000 )
```

- [Rsamtools](http://www.bioconductor.org/packages/release/bioc/html/Rsamtools.html)

---

# Read counting


```r
library( "GenomicAlignments" )
register( MulticoreParam( workers=4 ) )
# takes e.g. ~30 minutes per sample for 40 million PE reads 
se <- summarizeOverlaps( features=exonsByGene,
                        reads=bamLst,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
```

- [GenomicAlignments](http://www.bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
- [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html): python library
- [featureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html): `featureCounts()`

---

# SummarizedExperiment

<center><img src="summarizedexperiment.png" width=700/></center>

---

# Metadata stored without user effort


```r
metadata( rowData( se ) )[[ 1 ]][ 1:6 ]
```

```
## $`Db type`
## [1] "TranscriptDb"
## 
## $`Supporting package`
## [1] "GenomicFeatures"
## 
## $`Data source`
## [1] "BioMart"
## 
## $Organism
## [1] "Homo sapiens"
## 
## $`Resource URL`
## [1] "www.biomart.org:80"
## 
## $`BioMart database`
## [1] "ensembl"
```

---

# Add sample data


```r
samples <- read.csv( "sample_data.csv" )
colData( se ) <- DataFrame( samples )
```

---

# Exploratory data analysis (EDA)


```r
dds <- DESeqDataSet( se, ~ group + condition )
rld <- rlog( dds )
plotPCA( rld )
```

<center><img src="pca.png" width=500/></center>

---

# Differential expression analysis

- [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html)
differential exon usage
- [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) + voom normalization
- [DSS](http://www.bioconductor.org/packages/release/bioc/html/DSS.html)
- [BitSeq](http://www.bioconductor.org/packages/release/bioc/html/BitSeq.html)
transcript expression inference

---

# The generalized linear model

\[ K_{ij} \sim \text{NB}( \mu_{ij}, \alpha_i )  \]

- Read count $K_{ij}$ for gene *i* sample *j*.
- 2 parameter count distribution: mean $\mu$, dispersion $\alpha$

\[ \mu_{ij} = s_{ij} q_{ij} \]

- Normalized by size factor $s_{ij}$, with $q_{ij}$ remaining
- Often size factor $s_j$

---

# The generalized linear model

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

# Differential expression analysis


```r
# takes e.g. ~25 seconds for 8 samples, 35,000 genes
dds <- DESeq( dds )
res <- results( dds )
res <- results( dds, contrast=c("condition","trt","untrt") )
plotMA( res )
```

<center><img src="plotma.png" width=400/></center>

---

# Plot counts for a single gene


```r
gene <- rownames( res )[ order( res$pvalue ) ][ 1 ]
plotCounts( dds, gene, "condition" )
```

<center><img src="plotcounts.png" width=400/></center>

in *DESeq2* >= v1.5

---

# Normalization for sample-specific GC and transcript length

- [cqn](http://www.bioconductor.org/packages/release/bioc/html/cqn.html)
conditional quantile normalization
- [EDASeq](http://www.bioconductor.org/packages/release/bioc/html/EDASeq.html)

<center><img src="cqn.png" width=700/></center>

<center>cqn package vignette</center>

---

# Controlling for unknown batch

- [sva](http://www.bioconductor.org/packages/release/bioc/html/sva.html): `svaseq()`
surrogate variable analysis
- [RUVSeq](http://www.bioconductor.org/packages/release/bioc/html/RUVSeq.html):
remove unwanted variation

returns a matrix with columns which are surrogate variables

---

# Controlling for unknown batch

<center><img src="sva.png"/></center>

[Leek and Storey (2007)](http://dx.doi.org/10.1371/journal.pgen.0030161)

---

# ReportingTools


```r
rprt <- HTMLReport(shortName = "analysis",
                   title = "RNA-Seq analysis",
                   reportDirectory = "./reports")
publish(dds, rprt, pvalueCutoff=0.1,
        annotation.db="org.Hs.eg.db",
        factor = dds$condition,
        reportDir="./reports")
finish(rprt)
```

- [ReportingTools](http://www.bioconductor.org/packages/release/bioc/html/ReportingTools.html)

---

# ReportingTools

<center><img src="reptools.png" width=700/></center>

---

# Manual annotation

- [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html): `select()` function, works with the annotation packages `org.Hs.eg.db`:


```r
tab <- select(org.Hs.eg.db, genes, "SYMBOL", "ENSEMBL")
```

---

# Support

- mailing list:
    1. package name in title
    2. describe experiment, what's the *biological question*
    3. provide code
    4. sessionInfo()
- `browseVignettes("pkg")`
- `?function`

---

# Acknowledgments

- Bioconductor core team
- *DESeq*/*DEXSeq* team
    - Simon Anders
    - Alejandro Reyes
    - Wolfgang Huber
- Rafael Irizarry

