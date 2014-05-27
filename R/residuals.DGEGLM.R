residuals.DGEGLM <- function(object, type=c("deviance", "pearson"), ...) {
  y <- as.matrix(object$counts)
  mu <- as.matrix(object$fitted.values)
  theta <- 1/object$dispersion
  if(is.null(object$weights)) {
    wts <- rep(1, ncol(object$counts))
  } else {
    wts <- as.matrix(object$weights)
  }
  type <- match.arg(type)
  ymut <- cbind(y, mu, theta)

  res <- t(apply(ymut, 1, function(x) {
    yy <- as.vector(x[1:ncol(y)])
    mm <- as.vector(x[(ncol(y)+1):(ncol(y)+ncol(mu))])
    t <- x[length(x)]
    if(type=="deviance") {
      if(t==Inf) {
        d.res <- sqrt(pmax((poisson()$dev.resids)(yy, mm, wts), 0))
      } else {
        d.res <- sqrt(pmax((negative.binomial(theta=t)$dev.resids)(yy, pmax(mm, 1e-8), wts), 0))
      }
      return(ifelse(yy > mm, d.res, -d.res))
    } else if(type=="pearson") {
      if(t==Inf) {
        return((yy - mm) * sqrt(wts) / pmax(sqrt(poisson()$variance(mm)), 1))
      } else {
        return((yy - mm) * sqrt(wts) / pmax(sqrt(negative.binomial(theta=t)$variance(mm)), 1))
      }
    }
  }))
  return(res)
}
