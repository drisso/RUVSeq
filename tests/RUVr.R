library(RUVSeq)
library(edgeR)

mat <- matrix(data=rpois(100, lambda=10), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")

es <- newSeqExpressionSet(mat)

## compute edgeR residuals
x <- as.factor(rep(c("Ctl", "Trt"), each=5))
design <- model.matrix(~x)
y <- DGEList(counts=mat, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

## dimension of W
ks <- 1:5

## matrix
r1 <- lapply(ks, function(k) RUVr(mat, 1:10, k, res))

print(sapply(r1, function(x) dim(x$W)))
stopifnot(all(lapply(r1, function(x) dim(x$W)[2])==ks))

## already logged data
r1b <- lapply(ks, function(k) RUVr(log(mat+1), 1:10, k, res))
r1c <- lapply(ks, function(k) RUVr(log(mat+1), 1:10, k, res, isLog=TRUE))
r1d <- lapply(ks, function(k) RUVr(mat, 1:10, k, res, round=FALSE))

stopifnot(all(sapply(ks, function(i) all(r1[[i]]$W==r1c[[i]]$W))))
stopifnot(all(sapply(ks, function(i) all(r1d[[i]]$W==r1c[[i]]$W))))

stopifnot(all(sapply(ks, function(i) all(log(r1d[[i]]$normalizedCounts+1)-r1c[[i]]$normalizedCounts<1e-8))))

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVr(es, rownames(es)[1:10], k, res))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of zeros
mat <- matrix(data=rpois(100, lambda=2), ncol=10)
rownames(mat) <- paste("gene", 1:10, sep="")
r3 <- RUVr(mat, 1:10, k=1, res)
print(table(mat==0))
print(table(r3$normalizedCounts==0))
