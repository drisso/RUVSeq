\name{RUVSeq-package}
\alias{RUVSeq-package}
\alias{RUVSeq}
\docType{package}

\title{
Remove Unwanted Variation from RNA-Seq Data
}

\description{
This package implements the remove unwanted variation (RUV) methods of
Risso et al. (2014) for the normalization of RNA-Seq read counts between
samples.
}

\details{
\tabular{ll}{
Package: \tab RUVSeq\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2014-04-15\cr
License: \tab Artistic-2.0\cr
}

The \code{\link{RUVg}} function implements the RUVg normalization
procedure of Risso et al. (2014), by using control genes to remove
unwanted variation from the RNA-Seq read counts.

See also \code{\link{RUVr}} and \code{\link{RUVs}} for the "residual"
and "sample" methods, based, respectively, on residuals (e.g., deviance
residuals from a first-pass GLM regression of the unnormalized counts
on the covariates of interest) and replicate/negative control
samples for which the covariates of interest are constant. 
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
Davide Risso and Sandrine Dudoit
  
Maintainer: Davide Risso <\email{risso.davide@gmail.com}>
}

\seealso{
\code{\link{RUVg}}, \code{\link{RUVr}}, \code{\link{RUVs}}
}


\keyword{package}
