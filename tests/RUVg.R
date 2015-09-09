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

## already logged data
r1b <- lapply(ks, function(k) RUVg(log(mat+1), 1:10, k))
r1c <- lapply(ks, function(k) RUVg(log(mat+1), 1:10, k, isLog=TRUE))
r1d <- lapply(ks, function(k) RUVg(mat, 1:10, k, round=FALSE))

stopifnot(all(sapply(ks, function(i) all(r1[[i]]$W==r1b[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1d[[i]]$W==r1b[[i]]$W))))

stopifnot(all(sapply(ks, function(i) all(log(r1d[[i]]$normalizedCounts+1)-r1b[[i]]$normalizedCounts<1e-8))))

stopifnot(all(sapply(ks, function(i) all(r1c[[i]]$W==r1b[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1c[[i]]$normalizedCounts==r1b[[i]]$normalizedCounts))))

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVg(es, rownames(es)[1:10], k))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of zeros
mat <- matrix(data=rpois(100, lambda=2), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")
r3 <- RUVg(mat, 1:10, k=1)
print(table(mat==0))
print(table(r3$normalizedCounts==0))
