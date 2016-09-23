#' Make a matrix suitable for use with RUVSeq methods such as \code{\link{RUVs}}.
#'
#' Each row in the returned matrix corresponds to a set of replicate samples.
#' The number of columns is the size of the largest set of replicates; rows for
#' smaller sets are padded with -1 values.
#'
#' @param xs A vector indicating membership in a group.
#' @author Kamil Slowikowski
#' @seealso \code{\link{RUVs}}
#' @examples
#'  makeGroups(c("A","B","B","C","C","D","D","D","A"))
makeGroups <- function(xs) {
  xs <- factor(xs)
  groups <- matrix(-1, nrow = length(levels(xs)), ncol = max(table(xs)))
  for (i in 1:length(levels(xs))) {
    idxs <- which(xs == levels(xs)[i])
    groups[i,1:length(idxs)] <- idxs
  }
  groups
}

.isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
    !is.na(x) & abs(x - round(x)) < tol
}
