rm(list = ls()) 


library(tibble)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(GEOquery)
library(org.Hs.eg.db)
library(readxl)


gsename <- "GSE16011"
gse <- getGEO(gsename, destdir = ".") 
gpl <- getGEO('GPL8542', destdir = ".") 
gse <- getGEO(filename = 'GSE16011_series_matrix.txt.gz')

gpl <- gpl@dataTable@table
gpl <- gpl[!duplicated(gpl$ORF),]

hs <- org.Hs.eg.db
my.symbols <- as.character(gpl$ORF)
dt <- select(hs, keys = my.symbols,
             columns = c("ENTREZID", "SYMBOL"),
             keytype = "ENTREZID")
names(dt)[1] <- 'ORF'
gpl <- merge(gpl, dt, by = 'ORF')
gpl <- gpl %>%
  dplyr::select('ID', 'SYMBOL')

exprSet <- as.data.frame(exprs(gse))
exprSet$ID <- rownames(exprSet)
express <- merge(x = gpl, y = exprSet, by = "ID")
express$ID <- NULL
express[which(is.na(express), arr.ind = TRUE)] <- 0
express[1:5,1:5]
exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$SYMBOL),
                     FUN = max)
exprSet <- as.data.frame(exprSet)[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL

pd <- pData(gse)
pd <- subset(pd, select = c('title', 'geo_accession', 
                            "age at diagnosis:ch1" ,
                            'histology:ch1'))

stabs_1_6 <- read_excel("stabs_1-6.xlsx")
stabs_1_6 <- subset(stabs_1_6, select = c("Database number", "Gender", "Alive", "Survival (years)"))
stabs_1_6$title <- paste0('glioma ', stabs_1_6$`Database number`)
pd <- merge(pd, stabs_1_6, by = 'title')

names(pd) <- c("title", "Sample", "Age", "Grade", "number", "Gender", "OS", "OS_Time" )
pd$title <- NULL
pd <- na.omit(pd)

pd$Age <- as.numeric(as.character(pd$Age))
pd$Age <- round(pd$Age, 0)

table(pd$Grade)
pd <- subset(pd, pd$Grade == 'GBM (grade IV)')
pd$Grade <- NULL
pd$number <- NULL

table(pd$Gender)
pd$Gender  <- tools::toTitleCase(tolower(gsub(" ", "", pd$Gender)))

str(pd)

pd$OS <- ifelse(pd$OS == 'Dead', 1, 0)
pd$OS <- as.numeric(as.character(pd$OS))

pd$OS_Time <- chartr(old = ',', new = '.', x = pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))

pd <- subset(pd, pd$OS_Time > 0)

pd <- subset(pd, select=c("Sample", "Gender","Age","OS", "OS_Time"))

commSamples <- intersect(colnames(exprSet), pd$Sample)

exprSet[1:5, 1:5]
exprSet <- exprSet[, which(colnames(exprSet) %in% commSamples)]
pd <- pd[which(pd$Sample %in% commSamples),]

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = 'pd.csv')
