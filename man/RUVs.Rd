\name{RUVs-methods}
\docType{methods}
\alias{RUVs}
\alias{RUVs-methods}
\alias{RUVs,matrix,ANY,numeric,matrix-method}
\alias{RUVs,SeqExpressionSet,character,numeric,matrix-method}

\title{
Remove Unwanted Variation Using Replicate/Negative Control Samples
}

\description{
This function implements the RUVs method of Risso et al. (2014).
}

\usage{
RUVs(x, cIdx, k, scIdx, round=TRUE, epsilon=1, tolerance=1e-8)
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

\item{scIdx}{
A numeric matrix specifying the replicate samples for which to
compute the count differences used to estimate the factors of unwanted
variation (see details).
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
\item{\code{signature(x = "matrix", cIdx = "ANY", k = "numeric", scIdx = "matrix")}}{
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

\item{\code{signature(x = "SeqExpressionSet", cIdx = "character", k="numeric", scIdx = "matrix")}}{
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
The RUVs procedure performs factor analysis on a matrix of count
differences for replicate/negative control samples, for which the
biological covariates of interest are constant.
  
Each row of \code{scIdx} should correspond to a set of replicate
samples. The number of columns is the size of the largest set of
replicates; rows for smaller sets are padded with -1 values.

For example, if the sets of replicate samples are
{(1,11,21),(2,3),(4,5),(6,7,8)}, then \code{scIdx} should be

 1 11 21\cr
 2 3  -1\cr
 4 5  -1\cr
 6 7   8
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
Davide Risso (building on a previous version by Laurent Jacob).
}

\seealso{
\code{\link{RUVg}}, \code{\link{RUVr}}.
}

\examples{
library(zebrafishRNASeq)
data(zfGenes)

## run on a subset of genesfor time reasons 
## (real analyses should be performed on all genes)
genes <- rownames(zfGenes)[grep("^ENS", rownames(zfGenes))]
spikes <- rownames(zfGenes)[grep("^ERCC", rownames(zfGenes))]
set.seed(123)
idx <- c(sample(genes, 1000), spikes)
seq <- newSeqExpressionSet(as.matrix(zfGenes[idx,]))

# RUVs normalization
controls <- rownames(seq)
differences <- matrix(data=c(1:3, 4:6), byrow=TRUE, nrow=2)
seqRUVs <- RUVs(seq, controls, k=1, differences)

pData(seqRUVs)
head(normCounts(seqRUVs))

}
