# stableDE: Stability of results from differential expression analysis in RNA sequencing data

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
Package `samr' have been removed from CRAN. You may download the package from R Cran Archive 'https://cran.r-project.org/src/contrib/Archive/samr/' and install it manually. 
            
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



