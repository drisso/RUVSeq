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
r1d <- lapply(ks, function(k) RUVs(mat, 1:10, k, differences, round=FALSE))

stopifnot(all(sapply(ks, function(i) all(r1[[i]]$W==r1c[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1d[[i]]$W==r1c[[i]]$W))))

stopifnot(all(sapply(ks, function(i) all(log(r1d[[i]]$normalizedCounts+1)-r1c[[i]]$normalizedCounts<1e-8))))

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVs(es, rownames(es)[1:10], k, differences))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of zeros
mat <- matrix(data=rpois(100, lambda=2), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")
r3 <- RUVs(mat, 1:10, k=1, differences)
print(table(mat==0))
print(table(r3$normalizedCounts==0))

## make groups
factor1 <- rep(c("a", "b", "c"), each=3)
factor2 <- c(rep("a", 4), rep("b", 2), rep("c", 3))
factor3 <- rep(1:6, each=2)

makeGroups(factor1)
makeGroups(as.factor(factor1))
makeGroups(factor2)
makeGroups(as.factor(factor2))
makeGroups(factor3)
makeGroups(as.factor(factor3))


