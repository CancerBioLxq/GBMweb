rm(list = ls())
library(future)
library(future.apply)
library(survival)
library(limma)
library(sva)
library(dplyr)
library(caret)
library(compareGroups)

setwd("p4exp")

load(file = "DF_com.Rdata")

data <- DF_com


set.seed(1234)
trainIndex <- createDataPartition(data$OS, p = 0.5, list = FALSE)
trainData <- data[trainIndex, ]
tempData <- data[-trainIndex, ]
validIndex <- createDataPartition(tempData$OS, p = 0.5, list = FALSE)
validData <- tempData[validIndex, ]
testData <- tempData[-validIndex, ]

setwd("p4exp/data")

save(trainData, file = "trainData.Rdata")
save(validData, file = "validData.Rdata")
save(testData, file = "testData.Rdata")

trainData$resource <- "0trainData"
validData$resource <- "1validData"
testData$resource <- "2testData"

data1 <- rbind(trainData, validData, testData)

str(data1)

table(data1$OS)

data1$OS[data1$OS == 0] <- "alive"
data1$OS[data1$OS == 1] <- "dead"

table(data1$OS)

data1[1:8,1:8]
data1 <- data1 %>%
  dplyr::select(Gender, Age, OS, OS_Time, dataset, resource)

setwd("p4exp")

base_tab <- descrTable(resource ~ OS + OS_Time + Age + Gender + dataset,
                       data = data1)

export2word(base_tab, file='info.docx')
