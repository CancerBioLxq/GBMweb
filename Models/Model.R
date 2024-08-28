
rm(list = ls())
library(randomForest)
library(caret)
library(ROCR)
library(dplyr)
library(pROC)
library(e1071)    
library(class)   
library(rpart)  
library(MASS)
library(kknn)
library(glmnet)


train_data <- read.csv(file = "train_data.csv", row.names = 1)
valid_data <- read.csv(file = "valid_data.csv", row.names = 1)
test_data <- read.csv(file = "test_data.csv", row.names = 1)


train_data$Group <- factor(train_data$Group, levels = c("Over", "Under"))
valid_data$Group <- factor(valid_data$Group, levels = c("Over", "Under"))
test_data$Group <- factor(test_data$Group, levels = c("Over", "Under"))

finalData <- test_data %>%
  dplyr::select(Group)


best_auc <- 0
best_model <- NULL

for (lambda in c(0.01, 0.1, 1)) {
  glm_model <- glmnet(as.matrix(train_data[, -which(names(train_data) == "Group")]), 
                      train_data$Group, 
                      family = "binomial", alpha = 0, lambda = lambda)
  
  valid_prob <- as.vector(predict(glm_model, as.matrix(valid_data[, -which(names(valid_data) == "Group")]), type = "response"))
  
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for lambda=', lambda, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- glm_model
  }
}

finalData$LM <- as.vector(predict(best_model, as.matrix(test_data[, -which(names(test_data) == "Group")]), type = "response"))


best_auc <- 0
best_model <- NULL

shrinkage_values <- c(0, 0.01, 0.1)  

for (shrinkage in shrinkage_values) {
  cov_mat <- cov(train_data[, -which(names(train_data) == "Group")])
  regularized_cov_mat <- (1 - shrinkage) * cov_mat + shrinkage * diag(ncol(cov_mat))
  
  lda_model <- lda(Group ~ ., data = train_data, method = "moment", tol = shrinkage)
  
  valid_prob <- predict(lda_model, valid_data)$posterior[, 2]
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for shrinkage =', shrinkage, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- lda_model
  }
}

finalData$LDA <- predict(best_model, test_data)$posterior[, 2]


best_auc <- 0
best_model <- NULL

for (laplace in c(0, 1, 2)) {
  nb_model <- naiveBayes(Group ~ ., data = train_data, laplace = laplace)
  
  valid_prob <- predict(nb_model, valid_data, type = "raw")[, 2]
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for Laplace=', laplace, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- nb_model
  }
}

finalData$NB <- predict(best_model, test_data, type = "raw")[, 2]


best_auc <- 0
best_model <- NULL

for (cp in c(0.01, 0.05, 0.1)) {
  dt_model <- rpart(Group ~ ., data = train_data, control = rpart.control(cp = cp))
  
  valid_prob <- predict(dt_model, valid_data, type = "prob")[, 2]
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for cp=', cp, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- dt_model
  }
}

finalData$DT <- predict(best_model, test_data, type = "prob")[, 2]


best_auc <- 0
best_model <- NULL

for (n_trees in c(50, 100, 200)) {
  rf_model <- randomForest(Group ~ ., data = train_data, ntree = n_trees)
  
  valid_prob <- predict(rf_model, valid_data, type = "prob")[, 2]
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for', n_trees, 'trees is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- rf_model
  }
}

finalData$RF <- predict(best_model, test_data, type = "prob")[, 2]

# SVM - Linear Kernel
best_auc <- 0
best_model <- NULL

for (C in c(0.01, 0.1, 1)) {
  svm_linear_model <- svm(Group ~ ., data = train_data, kernel = "linear", cost = C, probability = TRUE)
  
  valid_prob <- attr(predict(svm_linear_model, valid_data, probability = TRUE), "probabilities")[, 2]
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for Linear SVM with C=', C, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- svm_linear_model
  }
}

test_prob <- attr(predict(best_model, test_data, probability = TRUE), "probabilities")[, 2]
finalData$SVM <- test_prob


performance_metrics <- data.frame(Model = character(), 
                                  AUC = numeric(), 
                                  Sensitivity = numeric(), 
                                  Specificity = numeric(), 
                                  PPV = numeric(), 
                                  NPV = numeric(), 
                                  stringsAsFactors = FALSE)

model_names <- colnames(finalData)[-1]

for (model_name in model_names) {
  roc_obj <- roc(finalData$Group, finalData[[model_name]])
  auc_value <- auc(roc_obj)
  
  cat("AUC for", model_name, "is", auc_value, "\n")
  
  cutff <- mean(finalData[[model_name]])
  
  print(cutff)
  
  predictions <- ifelse(finalData[[model_name]] > cutff, "Under", "Over")
  predictions <- factor(predictions, levels = c("Over", "Under"))
  
  confusion_matrix <- confusionMatrix(predictions, finalData$Group, positive = "Under")
  
  sensitivity_value <- confusion_matrix$byClass['Sensitivity']
  specificity_value <- confusion_matrix$byClass['Specificity']
  ppv_value <- confusion_matrix$byClass['Pos Pred Value']
  npv_value <- confusion_matrix$byClass['Neg Pred Value']
  
  performance_metrics <- rbind(performance_metrics, data.frame(
    Model = model_name,
    AUC = auc_value,
    Sensitivity = sensitivity_value,
    Specificity = specificity_value,
    PPV = ppv_value,
    NPV = npv_value
  ))
}

rownames(performance_metrics) <- performance_metrics$Model
performance_metrics$Model <- NULL

performance_metrics <- performance_metrics %>%
  dplyr::arrange(desc(AUC))

print(performance_metrics)

performance_metrics <- round(performance_metrics, 3)

write.csv(performance_metrics, file = "performance_metrics.csv")
