setMethod(
          f = "RUVg",
          signature = signature(x="matrix", cIdx="ANY", k="numeric"),
          definition = function(x, cIdx, k, drop=0, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE) {

            if(!isLog && !all(.isWholeNumber(x))) {
                warning(paste0("The expression matrix does not contain counts.\n",
                               "Please, pass a matrix of counts (not logged) or set isLog to TRUE to skip the log transformation"))
            }

            if(isLog) {
                Y <- t(x)
            } else {
                Y <- t(log(x+epsilon))
            }

            if (center) {
              Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale=FALSE))
            } else {
              Ycenter <- Y
            }
            if (drop >= k) {
              stop("'drop' must be less than 'k'.")
            }
            m <- nrow(Y)
            n <- ncol(Y)
            svdWa <- svd(Ycenter[, cIdx])
            first <- 1 + drop
            k <- min(k, max(which(svdWa$d > tolerance)))
            W <- svdWa$u[, (first:k), drop = FALSE]
            alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
            correctedY <- Y - W %*% alpha
            if(!isLog && all(.isWholeNumber(x))) {
                if(round) {
                    correctedY <- round(exp(correctedY) - epsilon)
                    correctedY[correctedY<0] <- 0
                } else {
                    correctedY <- exp(correctedY) - epsilon
                }
            }
            colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
            return(list(W = W, normalizedCounts = t(correctedY)))
          }
          )

setMethod(
          f = "RUVg",
          signature = signature(x="SeqExpressionSet", cIdx="character", k="numeric"),
          definition = function(x, cIdx, k, drop=0, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE) {
            if(isLog) {
              stop("SeqExpressionSet cannot contain log counts. Please, set isLog=FALSE.")
            }
            if(!all(cIdx %in% rownames(x))) {
              stop("'cIdx' must contain gene names present in 'x'")
            }
            if(k >= ncol(x)) {
              stop("'k' must be less than the number of samples in 'x'")
            }
            if(all(is.na(normCounts(x)))) {
              counts <- counts(x)
            } else {
              counts <- normCounts(x)
            }
            retval <- RUVg(counts, cIdx, k, drop, center, round, epsilon, tolerance, isLog=FALSE)
            newSeqExpressionSet(counts = counts(x),
                                normalizedCounts = retval$normalizedCounts,
                                phenoData = cbind(pData(x), retval$W)
                                )

          }
          )
