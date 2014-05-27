setMethod(
          f = "RUVs",
          signature = signature(x="matrix", cIdx="ANY", k="numeric", scIdx="matrix"),
          definition = function(x, cIdx, k, scIdx, round=TRUE, epsilon=1, tolerance=1e-8) {
            Y <- t(log(x+epsilon))
            scIdx <- scIdx[rowSums(scIdx > 0) >= 2, , drop = FALSE]
            Yctls <- matrix(0, prod(dim(scIdx)), ncol(Y))
            m <- nrow(Y)
            n <- ncol(Y)
            c <- 0
            for (ii in 1:nrow(scIdx)) {
              for (jj in 1:(ncol(scIdx))) {
                if (scIdx[ii, jj] == -1) 
                  next
                c <- c + 1
                Yctls[c, ] <- Y[scIdx[ii, jj], , drop = FALSE] - 
                  colMeans(Y[scIdx[ii, (scIdx[ii, ] > 0)], , drop = FALSE])
              }
            }
            Yctls <- Yctls[rowSums(Yctls) != 0, ]
            Y <- rbind(Y, Yctls)
            sctl <- (m + 1):(m + nrow(Yctls))
            svdRes <- svd(Y[sctl, ], nu = 0, nv = k)
            k <- min(k, max(which(svdRes$d > tolerance)))
            a <- diag(as.vector(svdRes$d[1:k]), ncol=k, nrow=k) %*% t(as.matrix(svdRes$v[, 1:k]))
            colnames(a) <- colnames(Y)
            W <- Y[, cIdx] %*% t(solve(a[, cIdx, drop = FALSE] %*% t(a[, cIdx, drop = FALSE]), a[, cIdx, drop = FALSE]))
            Wa <- W %*% a
            correctedY <- Y[1:m, ] - W[1:m, ] %*% a

            if(round) {
              correctedY <- round(exp(correctedY))
            } else {
              correctedY <- exp(correctedY)
            }
            W <- as.matrix(W[1:m,])
            colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
            return(list(W = W, normalizedCounts = t(correctedY)))
          }
          )

setMethod(
          f = "RUVs",
          signature = signature(x="SeqExpressionSet", cIdx="character", k="numeric", scIdx="matrix"),
          definition = function(x, cIdx, k, scIdx, round=TRUE, epsilon=1, tolerance=1e-8) {

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
            retval <- RUVs(counts, cIdx, k, scIdx, round, epsilon, tolerance)
            newSeqExpressionSet(counts = counts(x),
                                normalizedCounts = retval$normalizedCounts,
                                phenoData = cbind(pData(x), retval$W)
                                )

          }
          )
