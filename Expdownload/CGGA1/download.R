rm(list = ls()) 


library(tibble)
library(dplyr)
library(tidyr)
library(Biobase)
library(data.table)
library(stringr)
options(stringsAsFactors = FALSE)
library(GEOquery)
library(Biobase)

exprSet <- read.table('CGGA.mRNAseq_693.RSEM-genes.20200506.txt', 
                      header = TRUE, 
                      row.names = 1,
                      sep = "\t")
pheno <- read.table('CGGA.mRNAseq_693_clinical.20200506.txt', 
                    header = TRUE, 
                    sep = "\t")

table(pheno$Histology)

pheno <- subset(pheno, pheno$Histology == 'GBM')

table(pheno$PRS_type)

names(pheno)

pheno <- subset(pheno, select = c("CGGA_ID","Gender","Age", "Censor..alive.0..dead.1.","OS"))

str(pheno)

names(pheno) <- c("Sample",  "Gender","Age", "OS", "OS_Time")

pd <- pheno

head(pd)
str(pd)

table(pd$Gender)
pd$Gender  <- tools::toTitleCase(tolower(gsub(" ", "", pd$Gender)))
table(pd$Gender)

pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)

pd <- subset(pd, pd$OS_Time > 0)

pd <- na.omit(pd)

pd <- subset(pd, select=c("Sample", "Gender","Age","OS", "OS_Time"))

head(pd)

str(pd)

colSums(is.na(pd))

commSamples <- intersect(colnames(exprSet), pd$Sample)

exprSet[1:5, 1:5]
exprSet <- exprSet[, which(colnames(exprSet) %in% commSamples)]
pd <- pd[which(pd$Sample %in% commSamples),]

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = 'pd.csv')
