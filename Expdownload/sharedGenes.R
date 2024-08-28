rm(list=ls())

library(magrittr)
library(plyr)


setwd("CGGA1")
expr_CGGA1 <- read.csv("exprSet.csv", row.names = 1)

setwd("CGGA2")
expr_CGGA2 <- read.csv("exprSet.csv", row.names = 1)

setwd("GSE16011")
expr_GSE16011 <- read.csv("exprSet.csv", row.names = 1)

setwd("GSE43378")
expr_GSE43378 <- read.csv("exprSet.csv", row.names = 1)

setwd("GSE83300")
expr_GSE83300 <- read.csv("exprSet.csv", row.names = 1)

setwd("TCGAGBM")
expr_TCGA <- read.csv("exprSet.csv", row.names = 1)

datasets <- list(D1 = expr_CGGA1, 
                 D2 = expr_CGGA2, 
                 D3 = expr_GSE16011, 
                 D4 = expr_GSE43378,
                 D5 = expr_GSE83300,
                 D6 = expr_TCGA)

shared_genes <- NULL

for (df in datasets) {
  genes <- rownames(df)
  
  if (is.null(shared_genes)) {
    shared_genes <- genes
  } else {
    shared_genes <- intersect(shared_genes, genes)
  }
}

length(shared_genes)

expr_CGGA1 <- expr_CGGA1[which(rownames(expr_CGGA1) %in% shared_genes),]
expr_CGGA2 <- expr_CGGA2[which(rownames(expr_CGGA2) %in% shared_genes),]
expr_GSE16011 <- expr_GSE16011[which(rownames(expr_GSE16011) %in% shared_genes),]
expr_GSE43378 <- expr_GSE43378[which(rownames(expr_GSE43378) %in% shared_genes),]
expr_GSE83300 <- expr_GSE83300[which(rownames(expr_GSE83300) %in% shared_genes),]
expr_TCGA <- expr_TCGA[which(rownames(expr_TCGA) %in% shared_genes),]


save(shared_genes, file = "shared_genes.Rdata")

length(shared_genes)
