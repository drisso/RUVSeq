\name{residuals.DGEGLM}
\alias{residuals.DGEGLM}

\title{
Deviance and Pearson Residuals for the Negative Binomial Model of \pkg{\link{edgeR}}
}

\description{
This function implements the \code{\link[stats]{residuals}} method for the \pkg{\link{edgeR}} function \code{\link{glmFit}}.
}

\usage{
\method{residuals}{DGEGLM}(object, type = c("deviance", "pearson"), ...)
}

\arguments{
  \item{object}{
  An object of class \code{\link[=DGEGLM-class]{DGEGLM}} as created by the \code{\link{glmFit}} function of \pkg{\link{edgeR}}.
}
  \item{type}{
  Compute deviance or Pearson residuals.
}
  \item{...}{
  Additional arguments to be passed to the generic function.
}
}

\value{
  A genes-by-samples numeric matrix with the negative binomial residuals
  for each gene and sample.
}

\references{
McCullagh P, Nelder J (1989). \emph{Generalized Linear Models}. Chapman and Hall, New York.

Venables, W. N. and Ripley, B. D. (1999). \emph{Modern Applied Statistics with S-PLUS}. Third Edition. Springer.
}

\author{
 Davide Risso
}

\examples{
library(edgeR)
library(zebrafishRNASeq)
data(zfGenes)

## run on a subset genes for time reasons 
## (real analyses should be performed on all genes)
genes <- rownames(zfGenes)[grep("^ENS", rownames(zfGenes))]
spikes <- rownames(zfGenes)[grep("^ERCC", rownames(zfGenes))]
set.seed(123)
idx <- c(sample(genes, 1000), spikes)
seq <- newSeqExpressionSet(as.matrix(zfGenes[idx,]))

x <- as.factor(rep(c("Ctl", "Trt"), each=3))
design <- model.matrix(~x)
y <- DGEList(counts=counts(seq), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
head(res)
}
