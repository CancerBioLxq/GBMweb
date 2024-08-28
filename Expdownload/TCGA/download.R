rm(list = ls()) 
library(dplyr)
library(tidyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)
library(openxlsx)

cancer <- TCGAbiolinks:::getGDCprojects()$project_id
cancer <- str_subset(cancer, "TCGA")
cancer <- sort(cancer)
cancer_select <- "TCGA-GBM"

query_gbm <- GDCquery(
  project = cancer_select,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_gbm, method = "api", files.per.chunk = 50)

expr <- GDCprepare(query = query_gbm)
TPM <- as.data.frame(assay(expr, i = "tpm_unstrand"))

anno <- as.data.frame(expr@rowRanges@elementMetadata)
anno <- subset(anno, select = c('gene_name', 'gene_id'))
TPM$gene_id <- rownames(TPM)
TPM <- merge(anno, TPM, by = 'gene_id')
TPM$gene_id <- NULL

TPM[which(is.na(TPM), arr.ind = TRUE)] <- 0
TPM[1:5,1:5]

exprSet <- aggregate(x = TPM[, 2:ncol(TPM)], 
                     by = list(TPM$gene_name), 
                     FUN = max)
exprSet[1:5,1:5]

names(exprSet)[1] <- "ID"
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]
max(exprSet)

metad <- data.frame(names(exprSet))
colnames(metad)[1] <- 'sample'
metad$id <- substr(metad$sample, start = 1, stop = 12)
metad$tape <- substr(metad$sample, start = 14, stop = 16)
metad <- metad[order(metad$id), ]
metad <- subset(metad, metad$tape == "01A")
metad <- metad %>% distinct(id, .keep_all = TRUE)
exprSet <- exprSet[, which(colnames(exprSet) %in% metad$sample)]
colnames(exprSet) <- substr(colnames(exprSet), start = 1, stop = 12)
max(exprSet)

rowMeans_exprSet <- rowMeans(exprSet, na.rm = TRUE)
rowMedians_exprSet <- apply(exprSet, 1, median, na.rm = TRUE)

results <- data.frame(RowMeans = rowMeans_exprSet, RowMedians = rowMedians_exprSet)

results <- subset(results, results$RowMedians > 0.1)
results <- subset(results, results$RowMeans > 0.1)

exprSet <- exprSet[which(rownames(exprSet) %in% rownames(results)),]

my_data <- read.xlsx("mmc1.xlsx")
my_data <- subset(my_data, my_data$type == "GBM")
my_data <- my_data %>%
  dplyr::select("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis",
                "gender", "OS", "OS.time")

names(my_data) <- c("Sample", "Age", "Gender", "OS", "OS_Time")

my_data$Gender  <- tools::toTitleCase(tolower(gsub(" ", "", my_data$Gender)))
my_data$OS_Time <- my_data$OS_Time / 365
my_data$OS_Time <- round(my_data$OS_Time, 2)

pd <- subset(my_data, select=c("Sample", "Gender", "Age", "OS", "OS_Time"))

commSamples <- intersect(colnames(exprSet), pd$Sample)
exprSet <- exprSet[, which(colnames(exprSet) %in% commSamples)]
pd <- pd[which(pd$Sample %in% commSamples),]

colnames(exprSet) <- chartr(old = '-', new = '_', x=colnames(exprSet))
pd$Sample <- chartr(old = '-', new = '_', x=pd$Sample)

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = 'pd.csv')
