rm(list=ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(sva)

setwd("CGGA1")
load("exprSetNorm.Rdata")
expr1 <- as.data.frame(t(exprSetNorm))
expr1$genes <- rownames(expr1)                       
pd1 <- read.csv("pd.csv", header = T, row.names = 1)
pd1$dataset <- "CGGA1"

setwd("CGGA2")
load("exprSetNorm.Rdata")
expr2 <- as.data.frame(t(exprSetNorm))
expr2$genes <- rownames(expr2)                       
pd2 <- read.csv("pd.csv", header = T, row.names = 1)
pd2$dataset <- "CGGA2"

setwd("GSE16011")
load("exprSetNorm.Rdata")
expr3 <- as.data.frame(t(exprSetNorm))
expr3$genes <- rownames(expr3)                       
pd3 <- read.csv("pd.csv", header = T, row.names = 1)
pd3$dataset <- "GSE16011"

setwd("GSE43378")
load("exprSetNorm.Rdata")
expr4 <- as.data.frame(t(exprSetNorm))
expr4$genes <- rownames(expr4)                       
pd4 <- read.csv("pd.csv", header = T, row.names = 1)
pd4$dataset <- "GSE43378"

setwd("GSE83300")
load("exprSetNorm.Rdata")
expr5 <- as.data.frame(t(exprSetNorm))
expr5$genes <- rownames(expr5)                       
pd5 <- read.csv("pd.csv", header = T, row.names = 1)
pd5$dataset <- "GSE83300"

setwd("TCGAGBM")
load("exprSetNorm.Rdata")
expr6 <- as.data.frame(t(exprSetNorm))
expr6$genes <- rownames(expr6)                       
pd6 <- read.csv("pd.csv", header = T, row.names = 1)
pd6$dataset <- "TCGA"


exprset <- merge(expr1, expr2, by="genes")
exprset <- merge(exprset, expr3, by="genes")
exprset <- merge(exprset, expr4, by="genes")
exprset <- merge(exprset, expr5, by="genes")
exprset <- merge(exprset, expr6, by="genes")

rownames(exprset) <- exprset$genes
exprset$genes <- NULL
exprset <- as.data.frame(t(exprset))
exprset$Sample <- rownames(exprset)


pd <- rbind(pd1, pd2, pd3, pd4, pd5, pd6)
pd <- pd %>%
  dplyr::select(Sample, dataset, Age)

median_age <- median(pd$Age)
pd$Age <- ifelse(pd$Age <= median_age, "low", "high")

data <- merge(pd, exprset, by="Sample")
data[1:5,1:5]
rownames(data) <- data$Sample
data$Sample <- NULL
data[1:5,1:5]


pca_result <- prcomp(as.matrix(data[,-c(1:2)]), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$dataset <- data$dataset

custom_colors <- c("#6EACDA", "#000000", "#2ca02c", "#FB773C", "#BEC6A0", "#8c564b")
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "PCA before Combat Batch Correction", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = custom_colors)
print(p)

setwd("p4exp")
pdf(file = "PCA_before.pdf", height = 5, width = 5)
print(p)
dev.off()


expr_data <- as.matrix(data[, -(1:2)])
batch_info <- data$dataset
gender_info <- data$Age
mod <- model.matrix(~gender_info)
expr_data_t <- t(expr_data)
combat_data <- ComBat(dat = expr_data_t, batch = batch_info, mod = mod, par.prior = TRUE, prior.plots = FALSE)
combat_data_t <- t(combat_data)


pca_result <- prcomp(combat_data_t, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$dataset <- data$dataset
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(size = 2) +
  labs(title = "PCA after Combat Batch Correction", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = custom_colors)
print(p)

pdf(file = "PCA_after.pdf", height = 5, width = 5)
print(p)
dev.off()


combat_data_t <- as.data.frame(combat_data_t)
combat_data_t$Sample <- rownames(combat_data_t)
clinicalDT <- rbind(pd1, pd2, pd3, pd4, pd5, pd6)
median_age <- median(clinicalDT$Age)
clinicalDT$Age <- ifelse(clinicalDT$Age <= median_age, "low", "high")
DF_com <- merge(clinicalDT, combat_data_t, by="Sample")
DF_com[1:7,1:7]
rownames(DF_com) <- DF_com$Sample
DF_com$Sample <- NULL

save(DF_com, file = "DF_com.Rdata")
