
rm(list = ls())


library(ggplot2)
library(tidyr)
library(dplyr)


load(file = "selected_variables.Rdata")


TME1 <- read.csv("MCPresults.csv", row.names = 1)   
TME1 <- as.data.frame(t(TME1))
TME1$genes <- rownames(TME1)

pd1 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd1$dataset <- "CGGA1"

TME2 <- read.csv("MCPresults.csv", row.names = 1)   
TME2 <- as.data.frame(t(TME2))
TME2$genes <- rownames(TME2)

pd2 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd2$dataset <- "CGGA2"

TME3 <- read.csv("MCPresults.csv", row.names = 1)  
TME3 <- as.data.frame(t(TME3))
TME3$genes <- rownames(TME3)

pd3 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd3$dataset <- "GSE16011"

TME4 <- read.csv("MCPresults.csv", row.names = 1)  
TME4 <- as.data.frame(t(TME4))
TME4$genes <- rownames(TME4)

pd4 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd4$dataset <- "GSE43378"

TME5 <- read.csv("MCPresults.csv", row.names = 1)  
TME5 <- as.data.frame(t(TME5))
TME5$genes <- rownames(TME5)

pd5 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd5$dataset <- "GSE83300"

TME6 <- read.csv("MCPresults.csv", row.names = 1)  
TME6 <- as.data.frame(t(TME6))
TME6$genes <- rownames(TME6)

pd6 <- read.csv("pd.csv", header = TRUE, row.names = 1)
pd6$dataset <- "TCGA"


TME <- merge(TME1, TME2, by = "genes")
TME <- merge(TME, TME3, by = "genes")
TME <- merge(TME, TME4, by = "genes")
TME <- merge(TME, TME5, by = "genes")
TME <- merge(TME, TME6, by = "genes")

rownames(TME) <- TME$genes
TME$genes <- NULL

TME <- as.data.frame(t(TME))
TME$Sample <- rownames(TME)


dt1 <- read.csv(file = "trainSet.csv", row.names = 1)
dt2 <- read.csv(file = "validSet.csv", row.names = 1)
dt3 <- read.csv(file = "testSet.csv", row.names = 1)

exprset <- rbind(dt1, dt2, dt3)
exprset$OS <- NULL
exprset$OS_Time <- NULL
exprset$Sample <- rownames(exprset)


pd <- rbind(pd1, pd2, pd3, pd4, pd5, pd6)
pd$Age <- pd$Age / 10

data <- merge(pd, exprset, by = "Sample")
data <- merge(data, TME, by = "Sample")

rownames(data) <- data$Sample
data$Sample <- NULL


multi_data <- data
save(multi_data, file = "multi_data.Rdata")
