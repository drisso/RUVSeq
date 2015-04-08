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

## SeqExpressionSet
r2 <- lapply(ks, function(k) RUVr(es, rownames(es)[1:10], k, res))

print(sapply(r2, function(x) dim(pData(x))))
stopifnot(all(lapply(r2, function(x) dim(pData(x))[2])==ks))

## check handling of log counts
r3 <- RUVr(mat, 1:10, k=1, res)
r4 <- RUVr(log(mat+1), 1:10, k=1, res)

print(r3$W)
print(r4$W)
r3 <- RUVr(mat, 1:10, k=1, res)
r4 <- RUVr(log(mat+1), 1:10, k=1, res)

print(r3$W)
print(r4$W)

stopifnot(identical(r3, r4))
