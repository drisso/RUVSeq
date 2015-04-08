library(RUVSeq)

mat <- matrix(data=rpois(100, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")

es <- newSeqExpressionSet(mat)

## dimension of W
ks <- 1:5

## matrix
r1 <- lapply(ks, function(k) RUVg(mat, 1:10, k))

print(sapply(r1, function(x) dim(x$W)))
stopifnot(all(lapply(r1, function(x) dim(x$W)[2])==ks))

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVg(es, rownames(es)[1:10], k))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of log counts
r3 <- RUVg(mat, 1:10, k=1)
r4 <- RUVg(log(mat+1), 1:10, k=1)

print(r3$W)
print(r4$W)

stopifnot(identical(r3, r4))

