
rm(list = ls())


trainsamples <- read.csv("trainSet.csv", header = TRUE, row.names = 1)
validsamples <- read.csv("validSet.csv", header = TRUE, row.names = 1)
testsamples <- read.csv("testSet.csv", header = TRUE, row.names = 1)


load(file = "multi_data.Rdata")


data <- multi_data %>%
  dplyr::select(-c("dataset")) %>%
  dplyr::select(OS, OS_Time, everything())


data$Gender <- ifelse(data$Gender == "Female", 1, 0)


trainData <- data[rownames(data) %in% rownames(trainsamples), ]
validData <- data[rownames(data) %in% rownames(validsamples), ]
testData  <- data[rownames(data) %in% rownames(testsamples), ]


write.csv(trainData, file = "trainSet.csv")
write.csv(validData, file = "validSet.csv")
write.csv(testData, file = "testSet.csv")
