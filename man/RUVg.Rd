\name{RUVg-methods}
\docType{methods}
\alias{RUVg}
\alias{RUVg-methods}
\alias{RUVg,matrix,ANY,numeric-method}
\alias{RUVg,SeqExpressionSet,character,numeric-method}

\title{
Remove Unwanted Variation Using Control Genes
}

\description{
This function implements the RUVg method of Risso et al. (2014).
}

\usage{
RUVg(x, cIdx, k, drop=0, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8)
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

\item{drop}{
 The number of singular values to drop in the estimation of the factors
 of unwanted variation. This number is usually zero, but might be set to
 one if the first  singular value captures the effect of interest. It
 must be less than \code{k}.
 }

\item{center}{
If \code{TRUE}, the counts are centered, for each gene, to have mean
zero across samples. This is important to ensure that the first singular
value does not capture the average gene expression.
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
\item{\code{signature(x = "matrix", cIdx = "ANY", k = "numeric")}}{
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

\item{\code{signature(x = "SeqExpressionSet", cIdx = "character", k="numeric")}}{
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
The RUVg procedure performs factor analysis of the read counts based on
a suitably-chosen subset of negative control genes known a priori not
be differentially expressed (DE) between the samples under consideration.

Several types of controls can be used, including housekeeping genes,
spike-in sequences (e.g., ERCC), or ``in-silico'' empirical controls
(e.g., least significantly  DE genes based on a DE analysis performed
prior to RUV  normalization).

Note that one can relax the negative control gene assumption by
  requiring instead the identification of a set of positive or negative
  controls, with a priori known expression fold-changes between
  samples. RUVg can then simply be applied to control-centered log counts, as
  detailed in the vignette. 
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
\code{\link{RUVr}}, \code{\link{RUVs}}.
}

\examples{
library(zebrafishRNASeq)
data(zfGenes)

## run on a subset of genes for time reasons 
## (real analyses should be performed on all genes)
genes <- rownames(zfGenes)[grep("^ENS", rownames(zfGenes))]
spikes <- rownames(zfGenes)[grep("^ERCC", rownames(zfGenes))]
set.seed(123)
idx <- c(sample(genes, 1000), spikes)
seq <- newSeqExpressionSet(as.matrix(zfGenes[idx,]))

# RUVg normalization
seqRUVg <- RUVg(seq, spikes, k=1)

pData(seqRUVg)
head(normCounts(seqRUVg))

plotRLE(seq, outline=FALSE, ylim=c(-3, 3))
plotRLE(seqRUVg, outline=FALSE, ylim=c(-3, 3))

barplot(as.matrix(pData(seqRUVg)), beside=TRUE)
}
