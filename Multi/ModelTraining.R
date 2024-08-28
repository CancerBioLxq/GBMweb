#=======================================================

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

#=======================================================

train_data <- read.csv(file = "train_data.csv", row.names = 1)
valid_data <- read.csv(file = "valid_data.csv", row.names = 1)
test_data <- read.csv(file = "test_data.csv", row.names = 1)


train_data$Group <- factor(train_data$Group, levels = c("Over", "Under"))
valid_data$Group <- factor(valid_data$Group, levels = c("Over", "Under"))
test_data$Group <- factor(test_data$Group, levels = c("Over", "Under"))

finalData <- test_data %>%
  dplyr::select(Group)


clinical_columns <- names(train_data)[2:3]
rna_columns <- names(train_data)[4:14]
tme_columns <- names(train_data)[15:24]


train_clinical <- train_data[, clinical_columns]
valid_clinical <- valid_data[, clinical_columns]
test_clinical <- test_data[, clinical_columns]

train_rna <- train_data[, rna_columns]
valid_rna <- valid_data[, rna_columns]
test_rna <- test_data[, rna_columns]

train_tme <- train_data[, tme_columns]
valid_tme <- valid_data[, tme_columns]
test_tme <- test_data[, tme_columns]

train_early <- cbind(train_clinical, train_rna, train_tme)
valid_early <- cbind(valid_clinical, valid_rna, valid_tme)
test_early <- cbind(test_clinical, test_rna, test_tme)

#=======================================================

best_auc <- 0
best_model <- NULL

for (lambda in c(0.01, 0.1, 1)) {
  cli_model <- glmnet(as.matrix(train_clinical), 
                      train_data$Group, 
                      family = "binomial", alpha = 0, lambda = lambda)
  
  valid_prob <- as.vector(predict(cli_model, as.matrix(valid_clinical), type = "response"))
  
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for lambda=', lambda, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- cli_model
  }
}

finalData$clinical <- as.vector(predict(best_model, as.matrix(test_clinical), type = "response"))

#=======================================================

best_auc <- 0
best_model <- NULL

for (lambda in c(0.01, 0.1, 1)) {
  rna_model <- glmnet(as.matrix(train_rna), 
                      train_data$Group, 
                      family = "binomial", alpha = 0, lambda = lambda)
  
  valid_prob <- as.vector(predict(rna_model, as.matrix(valid_rna), type = "response"))
  
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for lambda=', lambda, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- rna_model
  }
}

finalData$rna <- as.vector(predict(best_model, as.matrix(test_rna), type = "response"))

#=======================================================

best_auc <- 0
best_model <- NULL

for (lambda in c(0.01, 0.1, 1)) {
  tme_model <- glmnet(as.matrix(train_tme), 
                      train_data$Group, 
                      family = "binomial", alpha = 0, lambda = lambda)
  
  valid_prob <- as.vector(predict(tme_model, as.matrix(valid_tme), type = "response"))
  
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for lambda=', lambda, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- tme_model
  }
}

finalData$tme <- as.vector(predict(best_model, as.matrix(test_tme), type = "response"))

#=======================================================

best_auc <- 0
best_model <- NULL

for (lambda in c(0.01, 0.1, 1)) {
  early_model <- glmnet(as.matrix(train_early), 
                        train_data$Group, 
                        family = "binomial", alpha = 0, lambda = lambda)
  
  valid_prob <- as.vector(predict(early_model, as.matrix(valid_early), type = "response"))
  
  roc_obj <- roc(valid_data$Group, valid_prob)
  valid_auc <- auc(roc_obj)
  
  cat('AUC for lambda=', lambda, 'is', valid_auc, '\n')
  
  if (valid_auc > best_auc) {
    best_auc <- valid_auc
    best_model <- early_model
  }
}

finalData$early <- as.vector(predict(best_model, as.matrix(test_early), type = "response"))

multi_m6 <- best_model

save(multi_m6, file = "multi_m6.Rdata")

#=======================================================

results_tt <- data.frame(h = numeric(),
                         i = numeric(),
                         j = numeric(),
                         AUC = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         PPV = numeric(),
                         NPV = numeric())

num1 <- seq(0, 1, by = 0.05)
num2 <- seq(0, 1, by = 0.05)
num3 <- seq(0, 1, by = 0.05)

for (h in num1) {
  for (i in num2) {
    for (j in num3) {
      valid_prob_cli <- as.vector(predict(cli_model, as.matrix(valid_clinical), type = "response"))
      valid_prob_rna <- as.vector(predict(rna_model, as.matrix(valid_rna), type = "response"))
      valid_prob_tme <- as.vector(predict(tme_model, as.matrix(valid_tme), type = "response"))
      
      valid_prob <- h * valid_prob_cli + i * valid_prob_rna + j * valid_prob_tme
      
      roc_obj <- roc(valid_data$Group, valid_prob)
      valid_auc <- auc(roc_obj)
      
      predictions <- ifelse(valid_prob > 0.5, "Over", "Under")
      
      confusion_matrix <- confusionMatrix(as.factor(predictions), as.factor(valid_data$Group), positive = "Over")
      
      sensitivity_value <- confusion_matrix$byClass['Sensitivity']
      specificity_value <- confusion_matrix$byClass['Specificity']
      ppv_value <- confusion_matrix$byClass['Pos Pred Value']
      npv_value <- confusion_matrix$byClass['Neg Pred Value']
      
      results_tt <- rbind(results_tt, data.frame(h = h, i = i, j = j, AUC = valid_auc, Sensitivity = sensitivity_value, Specificity = specificity_value, PPV = ppv_value, NPV = npv_value))
    }
  }
}

best_AUC_results <- results_tt[which.max(results_tt$AUC), ]
print(best_AUC_results)

finalData$late <- best_AUC_results$h * finalData$clinical + best_AUC_results$i * finalData$rna + best_AUC_results$j * finalData$tme

cols_to_normalize <- finalData[, 2:6]
normalized_cols <- as.data.frame(lapply(cols_to_normalize, function(x) (x - min(x)) / (max(x) - min(x))))
finalData[, 2:6] <- normalized_cols

#=======================================================

performance_metrics <- data.frame(Model = character(), AUC = numeric(), Sensitivity = numeric(), Specificity = numeric(), PPV = numeric(), NPV = numeric(), stringsAsFactors = FALSE)

model_names <- colnames(finalData)[-1]

for (model_name in model_names) {
  roc_obj <- roc(finalData$Group, finalData[[model_name]], levels = c("Under", "Over"))
  auc_value <- auc(roc_obj)
  
  cutffvalue <- 0.5
  
  predictions <- ifelse(finalData[[model_name]] > cutffvalue, "Over", "Under")
  
  confusion_matrix <- confusionMatrix(as.factor(predictions), as.factor(finalData$Group), positive = "Over")
  
  sensitivity_value <- confusion_matrix$byClass['Sensitivity']
  specificity_value <- confusion_matrix$byClass['Specificity']
  ppv_value <- confusion_matrix$byClass['Pos Pred Value']
  npv_value <- confusion_matrix$byClass['Neg Pred Value']
  
  performance_metrics <- rbind(performance_metrics, data.frame(Model = model_name, AUC = auc_value, Sensitivity = sensitivity_value, Specificity = specificity_value, PPV = ppv_value, NPV = npv_value))
}

auc_values <- subset(performance_metrics, select = c("Model", "AUC"))

rownames(performance_metrics) <- performance_metrics$Model
performance_metrics$Model <- NULL

print(performance_metrics)

performance_metrics <- round(performance_metrics, 3)
write.csv(performance_metrics, file = "performance_metrics.csv")

print(auc_values)

auc_values$Model <- factor(auc_values$Model, levels = c("clinical", "rna", "tme", "early", "late"))

custom_colors <- c("clinical" = "#4daf4b", "rna" = "#6EACDA", "tme" = "#C75B7A", "early" = "#984ea3", "late" = "#F4DEB3")

p <- ggplot(auc_values, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(AUC, digits = 3)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "", y = "AUC") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.position = "none")

print(p)

pdf(file = "auc.pdf", height = 5, width = 5)
print(p)
dev.off()
