library(RUVSeq)

mat <- matrix(data=rpois(100, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")

differences <- matrix(data=c(1:3, 4:6), byrow=TRUE, nrow=2)

es <- newSeqExpressionSet(mat)

## dimension of W
ks <- 1:4

## matrix
r1 <- lapply(ks, function(k) RUVs(mat, 1:10, k, differences))

print(sapply(r1, function(x) dim(x$W)))
stopifnot(all(lapply(r1, function(x) dim(x$W)[2])==ks))

## already logged data
r1b <- lapply(ks, function(k) RUVs(log(mat+1), 1:10, k, differences))
r1c <- lapply(ks, function(k) RUVs(log(mat+1), 1:10, k, differences, isLog=TRUE))

stopifnot(all(sapply(ks, function(i) all(r1[[i]]$W==r1b[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1[[i]]$normalizedCounts==r1b[[i]]$normalizedCounts))))

stopifnot(all(sapply(ks, function(i) all(r1c[[i]]$W==r1b[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1c[[i]]$normalizedCounts==r1b[[i]]$normalizedCounts))))

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVs(es, rownames(es)[1:10], k, differences))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of log counts
r3 <- RUVs(mat, 1:10, k=1, differences)
r4 <- RUVs(log(mat+1), 1:10, k=1, differences)

print(r3$W)
print(r4$W)

stopifnot(identical(r3, r4))

