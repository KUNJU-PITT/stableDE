
#' Normalization factors, library size and dispersion
#'
#' Calculate normalization factors, library size and dispersion by
#' edgeR_robust
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return tag_disp genewise dispersion arameter estimate for each gene for
#'         the negative binomial model
#' @return norm.factors normalization factors by TMM for each replicate
#' @return lib.size library size for each replicate
#' @examples
#'  bottomly <- RNAseq_data(sel_size = 3)
#'  counts <- bottomly$counts_evaluation
#'  group <- bottomly$group_evaluation
#'  dn <- disp_norm(counts, group)
#' @export



disp_norm <- function (counts, group) {

  library(edgeR)

  design <- model.matrix(~ group)

  d <- DGEList(counts = counts, group = group )
  d <- calcNormFactors(d)
  dw <- estimateGLMRobustDisp(d, design = design)
  tag_disp = dw$tagwise.dispersion

  list(tag_disp = tag_disp, norm.factors = d$samples$norm.factors,
       lib.size = d$samples$lib.size)
}


#' Average similaries for each considered perturbed size
#'
#' Calculate average similaries between sets of selected features from the generated
#' datasets for each perturbed size and the set of DE features from the
#' original dataset
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @param DEmethod differential expression analysis method: edgeR, edgeR_robust,
#'    DESeq_glm, DESeq2, samr_SAMseq, EBSeq or limma_voom.
#' @param alpha.vec numeric vector of perturbed size, default values is
#'    seq(0.01, 0.1, 0.01).
#' @param cut.point threshold for adjusted P values, default value is 0.05.
#' @param nrep number of repetition for each perturbed size, default value is 5.
#' @return ave_sim matrix of perturbed size and the corresponding average
#'    similaries.
#' @return pval.mat matrix of p-values and adjusted pvalues for the original
#'    dataset
#' @examples
#'  bottomly <- RNAseq_data(sel_size = 3)
#'  res <- ave_similarities(counts = bottomly$counts_evaluation,
#'      group = bottomly$group_evaluation, DEmethod = "edgeR")
#' @export

ave_similarities <- function (counts, group, DEmethod = c("edgeR", "edgeR_robust", "DESeq_glm", "DESeq2", "samr_SAMseq", "EBSeq", "limma_voom"),
                             alpha.vec = seq(0.01, 0.1, 0.01), cut.point = 0.05, nrep = 5) {


  DEmethod <- match.arg(DEmethod)
  DEmethod <- paste0(DEmethod, ".pfun", sel = "")

  stopifnot(cut.point < 1, cut.point > 0)
  stopifnot(nrep > 0)

  dn <- disp_norm(counts, group)
  dispersion <- dn$tag_disp


  lib.size = dn$lib.size *  dn$norm.factors
  counts.normalized <- t(t(counts) / lib.size)

  # normalized mean for two conditions
  A.mean <- rowMeans(counts.normalized[, group == levels(group)[1]])
  B.mean <- rowMeans(counts.normalized[, group == levels(group)[2]])

  # number of samples for each condition
  n.A <- sum(group == levels(group)[1])
  n.B <- sum(group == levels(group)[2])
  G <- nrow(counts)
  n <- n.A + n.B

  # average read counts
  mu <- matrix(A.mean, nrow = length(A.mean), ncol = n.A)
  mu <- cbind(mu, matrix(B.mean, nrow = length(B.mean), ncol = n.B))
  mu <- t(t(mu) * lib.size)

  # design matrix and call DE function to compute p values for origianl dataset

  orig.pvals <- do.call(DEmethod, list(counts = counts, group = group))
  set.DE <- rep(0, G)
  set.DE[orig.pvals[, "padj"] < cut.point] <- 1

  cor.vec = rep(0, length(alpha.vec))

  for (j in 1:length(alpha.vec)) {

    tmp.vec <- rep(0, length = nrep)
    for (k in 1:nrep) {

      counts.new <- counts
      for (i in 1:ncol(mu)) {
        tmp <- rnbinom(n = nrow(counts), mu = mu[,i], size= 1 / dispersion)
        flip.coin <- rbinom(G, 1, alpha.vec[j])
        counts.new[, i] <- counts[, i] * (1 - flip.coin) + flip.coin * tmp
      }

      simu.pvals <- do.call(DEmethod, list(counts = counts.new, group = group))
      set.DE.tmp <- rep(0, G)
      set.DE.tmp[simu.pvals[, "padj"] < cut.point] <- 1

      if(sum(set.DE.tmp) == 0 | sum(set.DE.tmp) == 0){
        tmp.vec[k] = 0
      }

      if(sum(set.DE.tmp) != 0 & sum(set.DE.tmp) != 0){
        tmp.vec[k] = max(0, cor(set.DE, set.DE.tmp))
      }

    }
    cor.vec[j] <- mean(tmp.vec)

  }

  ave_sim <- cbind(alpha.vec, cor.vec)
  colnames(ave_sim) <- c("alpha.vec", "cor")
  list(ave_sim = ave_sim, orig.pvals = orig.pvals )

}


#' AUCOR
#'
#' Calculate AUCOR
#'
#' @param object object from function "ave_similarites"
#'
#' @return smooth_line fitted smooth line for average of similaries against
#'    perturbed sizes
#' @return AUCOR AUCOR value
#' @examples
#'  bottomly <- RNAseq_data(sel_size = 3)
#'  ave <- ave_similarities(counts = bottomly$counts_evaluation,
#'  group = bottomly$group_evaluation, DEmethod = "edgeR")
#'  res <- AUCOR_fun(ave$ave_sim)
#' @export


AUCOR_fun <- function (object) {

  alpha <- object[, "alpha.vec"]
  corrs <- object[, "cor"]

  xx <- seq(0, max(alpha), length.out = 20)
  fit <- lm(corrs ~ poly(alpha, 2, raw=T))
  pre <- cbind(1, poly(xx, 2, raw=T)) %*% coef(fit)

  AUCOR <- 0
  for(i in 2:length(xx)){
    AUCOR <- AUCOR + 0.5 * (xx[i]-xx[i-1]) * (pre[i]+pre[i-1])
  }

  smooth_line <- cbind(x = xx, predict = pre)
  colnames(smooth_line) <- c("x", "predict")

  list(smooth_line = smooth_line, AUCOR = AUCOR)

}

#' Bottomly
#'
#' RNA-seq count data from Bottomly et al. (2011) that compares two
#' genetically homogeneous mice strains, C57BL/6J and DBA/2J. This dataset contains ten
#' and eleven replicates for each condition.
#'
#' @param sel_size number of replicates of each condition for evaluation set,
#'    default value is 3.
#' @return counts_evaluation  matrix of evaluation RNAseq dataset
#' @return group_evaluation factors of evaluation RNAseq dataset
#' @return design_evaluation design matrix of evaluation RNAseq dataset
#' @return counts_verification  matrix of verification RNAseq dataset
#' @return group_verification factors of verification RNAseq dataset
#' @return design_verification design matrix of verification RNAseq dataset
#' @return counts_full  matrix of full RNAseq dataset
#' @return group_full factors of full RNAseq dataset
#' @return design_full design matrix of full RNAseq dataset
#' @examples bottomly <- RNAseq_data(sel_size = 3)
#' @export
#'
RNAseq_data = function (sel_size = 3) {

  data(bottomly)
  counts <- bottomly
  group <- as.factor(c(rep("C57BL/6J", 10), rep("DBA/2J", 11)))

  # filter features with small expression levels
  cpm.counts <- 10^6 * t(t(counts) / colSums(counts))
  counts <- counts[rowSums(cpm.counts > 1) >= (ncol(counts)/2), ]

  # number of replicates for each condition
  n1 <- sum(group == levels(group)[1])
  n2 <- sum(group == levels(group)[2])

  # randomly split datasets into evaluation and verification parts
  index <- c(sample(which(group == levels(group)[1]), sel_size, replace=F),
            sample(which(group == levels(group)[2]), sel_size, replace=F))

  counts_evaluation <- as.matrix(counts[, index])
  group_evaluation <- group[index]
  design_evaluation <- model.matrix(~ group_evaluation)

  counts_verification <- as.matrix(counts[, -index])
  group_verification <- group[-index]
  design_verification <- model.matrix(~ group_verification)

  counts_full <- as.matrix(counts)
  group_full <- group
  design_full <- model.matrix(~ group_full)

  list(counts_evaluation = counts_evaluation[1:1000, ], group_evaluation = group_evaluation,
       design_evaluation = design_evaluation, counts_verification = counts_verification,
       group_verification = group_verification, design_verification = design_verification,
       counts_full = counts_full, group_full = group_full, design_full = design_full)

}


#' DESeq with glm test
#'
#' Implement DESeq method with glm test (Anders and Huber, 2010) with small
#' modification from
#'  "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- DESeq_glm.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal

DESeq_glm.pfun <-
  function(counts, group)
  {
    ## implement DESeq via the NB GLMs ##
    library(DESeq)
    design <- model.matrix(~ group)

    de <- newCountDataSet(counts, group)
    de <- estimateSizeFactors(de)
    de <- estimateDispersions(de)
    fit1 = fitNbinomGLMs(de, count ~ group )
    fit0 = fitNbinomGLMs(de, count ~ 1 )
    pval <- nbinomGLMTest( fit1, fit0 )
    padj <- p.adjust(pval, method="BH" )
    cbind(pval = pval, padj = padj)
  }


#' DESeq2
#'
#' Implement DESeq2 (Love et al., 2014) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- DESeq2.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#'
#' @keywords internal
#'
DESeq2.pfun <-
  function(counts, group)
  {
    ## implement DESeq2 ##
    library(DESeq2)

    colData <- data.frame(group)
    dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
    colData(dse)$group <- as.factor(colData(dse)$group)
    dse <- DESeq(dse, quiet=T)
    res <- results(dse)
    out <- cbind(pval = res$pvalue, padj = res$padj)
    out[is.na(out)] <- 1
    out
  }


#' edgeR
#'
#' Implement edgeR with standard pipeline (Robinson et al., 2010) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- edgeR.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal

edgeR.pfun <-
  function(counts, group, prior.df=10)
  {
    ## edgeR standard pipeline ##
    library(edgeR)
    design <- model.matrix(~ group)
    d <- DGEList(counts = counts, group = group )
    d <- calcNormFactors(d)
    d <- estimateGLMCommonDisp(d, design = design)
    d <- estimateGLMTrendedDisp(d,design=design)
    d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
    f <- glmFit(d, design = design)
    lr <- glmLRT(f, coef=2)
    pval = lr$table$PValue
    padj = p.adjust(pval, "BH")
    cbind(pval = pval, padj = padj)
  }


#' edgeR_robust
#'
#' Implement edgeR_robust with standard pipeline (Zhou et al., 2014) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- edgeR_robust.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal
#'
edgeR_robust.pfun <-
  function(counts, group, prior.df=10)
  {
    ## edgeR-robsut pipeline ##
    library(edgeR)
    design <- model.matrix(~ group)
    d <- DGEList(counts = counts, group = group )
    d <- calcNormFactors(d)
    dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6)
    fw <- glmFit(dw, design=design)
    lrw <- glmLRT(fw,coef=2)
    pval = lrw$table$PValue
    padj = p.adjust(pval, "BH")
    cbind(pval = pval, padj = padj)
  }

#' limma_Voom
#'
#' Implement limma_Voom with standard pipeline (Law et al., 2014) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- limma_voom.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal
#'
limma_voom.pfun <-
  function(counts, group)
  {
    ## limma voom pipeline ##
    library(limma)
    design <- model.matrix(~ group)
    nf <- calcNormFactors(counts)
    y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf)
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    pval <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")$P.Value
    padj <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")$adj.P.Val
    cbind(pval = pval, padj = padj)
  }

#' SAMseq
#'
#' Implement SAMseq with standard pipeline (Li and Tibshirani,2013) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'     bottomly <- RNAseq_data(sel_size = 3)
#'     res <- samr_SAMseq.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal
#'
samr_SAMseq.pfun <-
  function(counts, group)
  {
    ## SAMseq pipeline ##
    library(samr)
    f <- SAMseq(counts, group, resp.type = "Two class unpaired", fdr.output = 1)
    f.table = rbind(f$siggenes.table$genes.up, f$siggenes.table$genes.lo)
    fdr = rep(NA, nrow(counts)) #contains NA value
    fdr[as.numeric(f.table[, "Gene Name"])] <- as.numeric(f.table[, "q-value(%)"])/100
    padj <- fdr
    cbind(pval = pval, padj = padj) #pval, padj are identical
  }

#' EBSeq.pfun
#'
#' Implement EBSeq.pfun with standard pipeline (Leng et al., 2013) with small
#' modification from
#' "http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/".
#'
#' @param counts a matrix of raw read counts
#' @param group factor that giving the experimental condition for each sample
#' @return pval vector of p-values
#' @return padj vector of adjusted p-values by Benjamini-Hochberg procedure
#' @examples
#'   bottomly <- RNAseq_data(sel_size = 3)
#'   res <- EBSeq.pfun(bottomly$counts_evaluation, bottomly$group_evaluation)
#'
#' @keywords internal
#'
 EBSeq.pfun <-
  function(counts, group)
  {
    ## EBSeq pipeline ##
    library(EBSeq)
    group <- as.factor(group)
    sizes = MedianNorm(counts)
    f <- EBTest(Data = counts, Conditions = group, sizeFactors = sizes, maxround = 5)
    pp = GetPPMat(f)
    padj <- fdr <- 1 - pp[, "PPDE"]
    if(!length(padj) == nrow(counts)) #check rm 0 counts
    {
      i <- match(names(padj), rownames(counts))
      Padj <- rep(NA, nrow(counts))
      Padj[i] <- padj
      padj <- Padj
    }
    cbind(pval = padj, padj = padj) #pval, padj are identical
  }


