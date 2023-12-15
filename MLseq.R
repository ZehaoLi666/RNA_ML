library(MLSeq)
library(S4Vectors)


# add classification info by S4Vectors 
class_1<-DataFrame(condition=factor(rep(c("R","RT","M","T","S","G"), c(20,5,5,13,30,5))))
class_1 

# load the data and transfer the class of values from character to numeric 
selected_genes <- readRDS("selected_genes.rds")
head(selected_genes)
selected_genes <- sapply(selected_genes,as.numeric)
rownames_data <- rownames(selected_genes)
rownames(selected_genes) <- rownames_data
head(selected_genes) 
class(selected_genes[1,1])

# split data into training and testing sets
library(DESeq2)
set.seed(2128)

vars<-sort(apply(selected_genes,1,var,na.rm=TRUE),decreasing=TRUE)
data_1<-selected_genes[names(vars)[1:100],] 
nTest<-ceiling(ncol(data_1)* 0.3)
ind_1 <- sample(ncol(data_1),nTest,FALSE)


data.train_1 <- as.matrix(data_1[ ,-ind])
data.test_1 <- as.matrix(data_1[ ,ind])
classtr_1 <- DataFrame(condition = class_1[-ind,])
classts_1 <- DataFrame(condition = class_1[ind,])
head(data.train_1)

# generate DES matrix for next analysis 
data.trainS4_1=DESeqDataSetFromMatrix(countData=data.train_1, colData=classtr_1, design= formula(~condition))
data.testS4_1 = DESeqDataSetFromMatrix(countData = data.test_1, colData = classts_1,  design = formula(~condition))


# set up model parateters 
fit_1 <- classify(data = data.trainS4_1, method = "svmRadial",preProcessing = "deseq-rlog", ref = "T",  control = trainControl(method = "repeatedcv", number = 2, repeats = 2, classProbs = TRUE))
fit.svm <- classify(data = data.trainS4_1, method = "knn",
preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
control = trainControl(method = "repeatedcv", number = 5,
repeats = 10, classProbs = TRUE))


show(fit.svm)





library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)








filepath <- system.file("extdata/cervical.txt", package = "MLSeq")

cervical <- read.table(filepath, header=TRUE)
head(cervical[ ,1:10])
class <- DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))
set.seed(2128)
vars <- sort(apply(cervical, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- cervical[names(vars)[1:100], ]
nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)
data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
class(data.train[1,1])
data.test
class(
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,design = formula(~condition))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,design = formula(~condition))
head(data.train)



fit.svm <- classify(data = data.trainS4, method = "svmRadial",
preProcessing = "deseq-vst", ref = "T", tuneLength = 10,
control = trainControl(method = "repeatedcv", number = 5,
repeats = 10, classProbs = TRUE))
show(fit.svm)





















