\name{RUVr-methods}
\docType{methods}
\alias{RUVr}
\alias{RUVr-methods}
\alias{RUVr,matrix,ANY,numeric,matrix-method}
\alias{RUVr,SeqExpressionSet,character,numeric,matrix-method}

\title{
Remove Unwanted Variation Using Residuals
}

\description{
This function implements the RUVr method of Risso et al. (2014).
}

\usage{
RUVr(x, cIdx, k, residuals, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8)
}

\arguments{
\item{x}{
Either a genes-by-samples numeric matrix or a
\linkS4class{SeqExpressionSet} object containing the read counts.
}

\item{cIdx}{
A character, logical, or numeric vector indicating the subset of genes to be used as negative controls in the estimation of the factors of unwanted variation.
}

\item{k}{
The number of factors of unwanted variation to be estimated from the data.
}

\item{residuals}{
A genes-by-samples matrix of residuals obtained from a first-pass regression of the counts on the covariates of interest, usually the negative binomial deviance residuals obtained from \pkg{\link{edgeR}} with the \code{\link[RUVSeq:residuals.DGEGLM]{residuals}} method.
}

\item{center}{
If \code{TRUE}, the residuals are centered, for each gene, to have mean zero across samples.
}

\item{round}{
If \code{TRUE}, the normalized measures are rounded to form pseudo-counts.
}

 \item{epsilon}{
A small constant (usually no larger than one) to be added to the counts prior to the log transformation to avoid problems with log(0).
}

\item{tolerance}{
Tolerance in the selection of the number of positive singular values, i.e., a singular value must be larger than \code{tolerance} to be considered positive.
}
}

\section{Methods}{
\describe{
\item{\code{signature(x = "matrix", cIdx = "ANY", k = "numeric", residuals = "matrix")}}{
It returns a list with 
\itemize{
\item
A samples-by-factors matrix with the estimated factors of unwanted variation (\code{W}).
\item
The genes-by-samples matrix of normalized expression measures (possibly
rounded) obtained by removing the factors of unwanted variation from the
original read counts (\code{normalizedCounts}).
}
}

\item{\code{signature(x = "SeqExpressionSet", cIdx = "character", k="numeric",
    residuals = "matrix")}}{
It returns a \linkS4class{SeqExpressionSet} with 
\itemize{
\item
The normalized counts in the \code{normalizedCounts} slot.
\item
The estimated factors of unwanted variation as additional columns of the
\code{phenoData} slot.
}
}
}
}

\details{
The RUVr procedure performs factor analysis on residuals, such as deviance
residuals from a first-pass GLM regression of the counts on the
covariates of interest using \pkg{\link{edgeR}}. The counts may be either unnormalized or
normalized with a method such as upper-quartile (UQ) normalization. 
 }

\references{
D. Risso, J. Ngai, T. P. Speed, and S. Dudoit. 
Normalization of RNA-seq data using factor analysis of control genes or samples.
 \emph{Nature Biotechnology}, 2014. (In press).

D. Risso, J. Ngai, T. P. Speed, and S. Dudoit. The role of spike-in
standards in the normalization of RNA-Seq. In D. Nettleton and S. Datta,
editors, \emph{Statistical Analysis of Next Generation Sequence
  Data}. Springer, 2014. (In press).
}

\author{
Davide Risso
}

\seealso{
\code{\link{RUVg}}, \code{\link{RUVs}}, \code{\link[RUVSeq:residuals.DGEGLM]{residuals}}.
}

\examples{
library(edgeR)
library(zebrafishRNASeq)
data(zfGenes)

## run on a subset of genes for time reasons 
## (real analyses should be performed on all genes)
genes <- rownames(zfGenes)[grep("^ENS", rownames(zfGenes))]
spikes <- rownames(zfGenes)[grep("^ERCC", rownames(zfGenes))]
set.seed(123)
idx <- c(sample(genes, 1000), spikes)
seq <- newSeqExpressionSet(as.matrix(zfGenes[idx,]))

# Residuals from negative binomial GLM regression of UQ-normalized
# counts on covariates of interest, with edgeR
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
design <- model.matrix(~x)
y <- DGEList(counts=counts(seq), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# RUVr normalization (after UQ)
seqUQ <- betweenLaneNormalization(seq, which="upper")
controls <- rownames(seq)
seqRUVr <- RUVr(seqUQ, controls, k=1, res)

pData(seqRUVr)
head(normCounts(seqRUVr))
}
