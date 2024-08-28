rm(list=ls())
library(MCPcounter)

load("exprSetNorm.Rdata")
exprSet <- as.data.frame(t(exprSetNorm))
exprSet[1:6,1:6]

results <- MCPcounter.estimate(expression = exprSet, featuresType = "HUGO_symbols")

results <- as.data.frame(t(results))
results[1:6,1:6]

quantileNormalizationPositive <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  x_normalized <- ((x - Q1) / (Q3 - Q1))  - min((x - Q1) / (Q3 - Q1)) + 1
  return(x_normalized)
}

results_normalized <- apply(results, 2, quantileNormalizationPositive)
results_normalized <- as.data.frame(results_normalized)
colnames(results_normalized) <- gsub(" ", "", colnames(results_normalized))
results_normalized[1:5,1:5]

write.csv(results_normalized, file = "MCPresults.csv")
