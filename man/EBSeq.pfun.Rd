% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stableDE_fun.R
\name{EBSeq.pfun}
\alias{EBSeq.pfun}
\title{EBSeq.pfun}
\usage{
EBSeq.pfun(counts, group)
}
\arguments{
\item{counts}{a matrix of raw read counts}

\item{group}{factor that giving the experimental condition for each sample}
}
\value{
pval vector of p-values

padj vector of adjusted p-values by Benjamini-Hochberg procedure
}
\description{
Implement EBSeq.pfun with standard pipeline (Leng et al., 2013) with small
modification from
"http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
}
\examples{
  bottomly <- RNAseq_data(sel_size = 3)
  res <- EBSeq.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)

}
\keyword{internal}
