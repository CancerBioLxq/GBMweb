rm(list = ls()) 
library(future)
library(future.apply)
library(survival)
library(glmnet)
library(dplyr)

setwd("p4exp/data")

load(file = "trainData.Rdata")
load(file = "validData.Rdata")
load(file = "testData.Rdata")

trainData$Gender <- NULL
trainData$Age <- NULL
trainData$dataset <- NULL

validData$Gender <- NULL
validData$Age <- NULL
validData$dataset <- NULL

testData$Gender <- NULL
testData$Age <- NULL
testData$dataset <- NULL

survival_dat <- trainData %>% dplyr::select("OS", "OS_Time", everything())
colnames(survival_dat) <- gsub("\\-", "_", colnames(survival_dat))
colnames(survival_dat) <- gsub("\\.", "_", colnames(survival_dat))
colnames(survival_dat) <- gsub("\\:", "_", colnames(survival_dat))
colnames(survival_dat) <- gsub("\\@", "", colnames(survival_dat))
colnames(survival_dat) <- gsub("\\/", "", colnames(survival_dat))
colnames(survival_dat) <- gsub(" ", "", colnames(survival_dat))

covariates <- as.character(colnames(survival_dat))[-(1:2)]
univ_formulas <- sapply(covariates, function(x) {
  as.formula(paste('Surv(OS_Time, OS)~', x))
})
univ_models <- future_lapply(univ_formulas, function(x) {
  coxph(x, data = survival_dat)
})
univ_results <- lapply(univ_models, function(x) {
  x <- summary(x)
  p.value <- signif(x$wald["pvalue"], digits = 2)
  beta <- signif(x$coef[1], digits = 2)
  HR <- signif(x$coef[2], digits = 2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  res <- c(beta, HR, p.value)
  names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
  return(res)
})
res <- as.data.frame(t(do.call(cbind, univ_results)))
res$p.value <- as.numeric(as.character(res$p.value))
res$coef <- as.numeric(as.character(res$coef))

res$p_adjusted <- p.adjust(res$p.value, method = "BH")

write.csv(res, file = "res.csv")

res <- subset(res, res$p_adjusted < 0.05)
seleVarai <- rownames(res)
seleVarai <- c("OS", "OS_Time", seleVarai)

survival_dat <- survival_dat[,which(colnames(survival_dat) %in% seleVarai)]


x <- as.matrix(survival_dat[, -c(1, 2)])
y <- survival_dat$OS
time <- survival_dat$OS_Time

surv_object <- Surv(time, y)
cv_fit <- cv.glmnet(x, surv_object, family = "cox", alpha = 1)

lambda_values <- cv_fit$lambda

for (lambda in lambda_values) {
  model <- glmnet(x, surv_object, family = "cox", alpha = 1, lambda = lambda)
  coef_matrix <- as.matrix(coef(model))
  non_zero_genes <- sum(coef_matrix != 0)
  
  if (non_zero_genes >= 10 && non_zero_genes <= 20) {
    selected_lambda <- lambda
    break
  }
}

final_model <- glmnet(x, surv_object, family = "cox", alpha = 1, lambda = selected_lambda)
coefficients_matrix <- as.matrix(coef(final_model))
selected_variables <- rownames(coefficients_matrix)[coefficients_matrix != 0]

selected_variables <- c("OS", "OS_Time", selected_variables)

setwd("p4exp")

save(selected_variables, file="selected_variables.Rdata")
write.csv(selected_variables, file = "selected_variables.csv")

trainSet <- trainData[,which(colnames(trainData) %in% selected_variables)]
validSet <- validData[,which(colnames(validData) %in% selected_variables)]
testSet <-  testData[,which(colnames(testData) %in% selected_variables)]

setwd("p4exp/data")

write.csv(trainSet, file = "trainSet.csv")
write.csv(validSet, file = "validSet.csv")
write.csv(testSet, file = "testSet.csv")
