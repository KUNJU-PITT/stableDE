# stableDE: Stability of results from differential expression analysis in RNA sequencing data

## Abstract

As RNA-seq becomes the assay of choice for measuring gene expression levels, differential expression analysis has received extensive attentions of researchers. To date, for the evaluation of DE methods, most attention has been paid on validity. Yet another important aspect of DE methods, stability, is often overlooked. In this study, we empirically show the need of assessing stability of DE methods and propose a stability metric, called Area Under the Correlation curve (AUCOR), that generates the perturbed datasets by a mixture distribution and combines the information of similarities between sets of selected features from these perturbed datasets and the original dataset. Empirical results support that AUCOR can effectively rank the DE methods in terms of stability for given RNA-seq datasets. In addition, we explore how biological or technical factors from experiments and data analysis affect the stability of DE methods.

## Installation Guide

For the installation, the R package devtools is needed.

```{r}
install.packages("devtools")
library(devtools)
```
I recommend to install first the dependencies manually and then stableDE:

```{r}
pkg <- c('edgeR', 'limma', 'DESeq', 'DESeq2', 'EBSeq', 'limma')
if(any(!pkg %in% installed.packages()[, "Package"])) {
  new.pkg = pkg[!pkg %in% installed.packages()[, "Package"]]
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.pkg, dependencies = TRUE, ask = FALSE)
}
```
Package `samr' has been removed from CRAN. You may download the package from R Cran Archive 'https://cran.r-project.org/src/contrib/Archive/samr/' and install it manually. 
            
Then, you can install stableDE directly using devtools:
```{r}
devtools::install_github("linbingqing/stableDE")
```

## One simple example

```{r}
bottomly <- RNAseq_data(sel_size = 3)
ave <- ave_similarities(counts = bottomly$counts_evaluation,
group = bottomly$group_evaluation, DEmethod = "edgeR")
res <- AUCOR_fun(ave$ave_sim)
```
The input counts is a RNA-seq raw read counts matrix.

The other input group is a vector of factor which specifies the two groups in the matrix to be compared, corresponding to the columns in counts.



