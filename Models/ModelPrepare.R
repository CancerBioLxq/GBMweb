rm(list=ls())
setwd("p4exp/data")

library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(caret)
library(randomForestSRC)
library(fastDummies)

trainData <- read.csv("trainSet.csv", row.names = 1)
validData <- read.csv("validSet.csv", row.names = 1)
testData <- read.csv("testSet.csv", row.names = 1)

data <- rbind(trainData, validData, testData)

yearnumvalue <- 0.5

data_alive <- data %>%
  dplyr::filter(OS_Time > yearnumvalue) %>%
  mutate(across(everything(), as.numeric))

data_alive$Group <- "Over"

data_dead <- data %>%
  dplyr::filter(OS == 1, OS_Time < yearnumvalue) %>%
  mutate(across(everything(), as.numeric))

data_dead$Group <- "Under"

data2 <- rbind(data_alive, data_dead) %>%
  na.omit() %>%
  dplyr::select(Group, everything()) %>%
  mutate(Group = as.factor(Group))

data2$OS <- NULL
data2$OS_Time <- NULL

train_data <- data2[which(rownames(data2) %in% rownames(trainData)),]
valid_data <- data2[which(rownames(data2) %in% rownames(validData)),]
test_data <- data2[which(rownames(data2) %in% rownames(testData)),]

setwd("p4exp/eML6")

write.csv(train_data, file = "train_data.csv")
write.csv(valid_data, file = "valid_data.csv")
write.csv(test_data, file = "test_data.csv")
