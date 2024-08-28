rm(list = ls())
library(future)
library(future.apply)
library(survival)
library(limma)


load(file = "shared_genes.Rdata")

exprSet <- read.csv(file = "exprSet.csv", row.names = 1)
exprSet <- exprSet[which(rownames(exprSet) %in% shared_genes),]

FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exprSet_tpm <- apply(exprSet, 2, FPKM2TPM)
exprSet <- as.data.frame(t(exprSet_tpm))
exprSet[1:5,1:5]
dim(exprSet)

quantileNormalizationPositive <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  x_normalized <- ((x - Q1) / (Q3 - Q1))  - min((x - Q1) / (Q3 - Q1)) + 1
  return(x_normalized)
}

exprSet_normalized <- apply(exprSet, 2, quantileNormalizationPositive)
exprSet <- as.data.frame(exprSet_normalized)
exprSet[1:5,1:5]

exprSetNorm <- exprSet
save(exprSetNorm, file = "exprSetNorm.Rdata")
