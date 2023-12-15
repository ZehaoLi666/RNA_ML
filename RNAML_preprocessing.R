# import data 
pairEnd100 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd100_countDFeByg.csv", header=TRUE,sep=",")

pairEnd200 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd200_countDFeByg.csv", header=TRUE,sep=",")
pairEnd200 <- pairEnd200[,-1]

pairEnd300 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd300_countDFeByg.csv", header=TRUE,sep=",")
pairEnd300 <- pairEnd300[,-1]

pairEnd400 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd400_countDFeByg.csv", header=TRUE,sep=",")
pairEnd400 <- pairEnd400[,-1]

pairEnd401 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd401-450_countDFeByg.csv", header=TRUE,sep=",")
pairEnd401 <- pairEnd401[,-1]

pairEnd451 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd451-470_countDFeByg.csv", header=TRUE,sep=",")
pairEnd451 <- pairEnd451[,-1]

pairEnd481 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd481-500_countDFeByg.csv", header=TRUE,sep=",")
pairEnd481 <- pairEnd481[,-1]

pairEnd610 <- read.csv("~/lab/RNA_machine_learning/csv/pairEnd610_countDFeByg.csv", header=TRUE,sep=",")
pairEnd610 <- pairEnd610[,-1]

pairEnd <- cbind(pairEnd100,pairEnd200,pairEnd300,pairEnd400,pairEnd401,pairEnd451,pairEnd481,pairEnd610)

singleEnd <- read.csv("~/lab/RNA_machine_learning/csv/singleEnd_countDFeByg.csv", header=TRUE,sep=",")
singleEnd <- singleEnd[,-1]

All_data <- cbind(pairEnd,singleEnd)
rownames_Alldata <- All_data$X
All_data <- All_data[,-1]
All_data <- sapply(All_data,as.integer)
rownames(All_data) <- rownames_Alldata


# Desq2 
library(DESeq2)
library(S4Vectors)


All_data_matrix <- as.matrix(All_data)



coldata <- read.csv("~/lab/RNA_machine_learning/csv/anno.csv", row.names=1)
#coldata <- coldata[,c("condition")]

coldata$condition <- factor(coldata$condition)

data_des=DESeqDataSetFromMatrix(countData=All_data_matrix, colData=coldata, design= formula(~condition))

#rownames(coldata) <- sub("A", "1", rownames(coldata))
#rownames(coldata) <- sub("M", "2", rownames(coldata))
#rownames(coldata) <- sub("RT", "3", rownames(coldata))
#rownames(coldata) <- sub("R", "4", rownames(coldata))
#rownames(coldata) <- sub("S", "5", rownames(coldata))
#rownames(coldata) <- sub("T", "6", rownames(coldata))
all(rownames(coldata) %in% colnames(All_data_matrix))
all(rownames(coldata) == colnames(All_data_matrix))

dds <- DESeqDataSetFromMatrix(countData = All_data_matrix,
                              colData = coldata,
                              design = ~ condition)

# do the gene_wise_variances to select Top 100 gens 

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
gene_wise_variances <- assay(vsd)

vars<-sort(apply(gene_wise_variances,1,var,na.rm=TRUE),decreasing=TRUE)
data<-gene_wise_variances[names(vars)[1:100],]

# Divide data into training set and testing set 
data<-as.data.frame(data)
rownames_Selected_data <- rownames(data)
data <- sapply(data,as.integer)
rownames(data) <- rownames_Selected_data
data <- as.matrix(data)

nTest<-ceiling(ncol(data)*0.3)
ind<-sample(ncol(data),nTest, FALSE)

data.train<-data[,-ind]

data.test<-data[,ind]

classtr<-DataFrame(condition= coldata[-ind,])
classts<-DataFrame(condition= coldata[ind,])
data.trainS4=DESeqDataSetFromMatrix(countData=data.train,colData=classtr, design=formula(~condition))
data.testS4DESq=DESeqDataSetFromMatrix(countData=data.test,colData=classts, design=formula(~condition))



# run the training model 
#fit_1 <- classify(data = data.trainS4_1, method = "svmRadial",preProcessing = "deseq-rlog", ref = "T",  control = trainControl(method = "repeatedcv", number = 2, repeats = 2, classProbs = TRUE))
library(MLSeq)
set.seed(2128)
fit.voomNSC <- classify(data = data.trainS4, method = "voomNSC",
                    preProcessing = "deseq", ref = "T", 
                    control = voomControl(method="repeatedcv",number=100,repeats=100,tuneLength=20))
show(fit.voomNSC)
trained(fit.voomNSC)

### SVM1
fit.svm<-classify(data=data.trainS4,method="svmRadial", preProcessing= "deseq-vst",ref="S",tuneLength=10, control=trainControl(method="repeatedcv",number=5, repeats=10,classProbs=TRUE))
show(fit.svm)
trained(fit.svm)

# Define your parameter grid
### SVM_2 
param_grid <- list(C = c(0.1, 1, 10, 100), sigma = c(0.1, 0.5, 1))
fit.svm_2 <- classify(
  data = data.trainS4,
  method = "svmRadial",
  preProcessing = "deseq-vst",
  ref = "S",
  tuneParams = param_grid,
  control = trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    classProbs = TRUE
  )
)
show(fit.svm_2)
trained(fit.svm_2)
png("/bigdata/lerochlab/zli529/RNA_machine_learning/RESULTS/fit.svm_2.png")
plot(fit.svm_2)
dev.off()

### nnet 
fit.nn <- classify(
  data = data.trainS4,
  method = "nnet",
  preProcessing = "deseq-vst",
  ref = "S",
  control = trainControl(
    method = "repeatedcv",
    number = 100,
    repeats = 100,
    classProbs = TRUE
  )
)
show(fit.nn)
trained(fit.nn)
png("/bigdata/lerochlab/zli529/RNA_machine_learning/RESULTS/fit_nn.png")
plot(fit.nn)
dev.off()

### Training model (nn_tuned)
param_grid <- expand.grid(size = c(5, 10, 15), decay = c(0, 0.01, 0.1))
fit.nn_tuned <- classify(
  data = data.trainS4,
  method = "nnet",
  preProcessing = "deseq-vst",
  ref = "S",
  tuneParams = param_grid,
  control = trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10,
    classProbs = TRUE
  )
)
show(fit.nn_tuned)
trained(fit.nn_tuned)
png("/bigdata/lerochlab/zli529/RNA_machine_learning/RESULTS/fit.nn_tuned.png")
plot(fit.nn_tuned)
dev.off()
selectedGenes(fit.nn_tuned)



pred.nn_tuned <- predict(fit.nn_tuned, data.testS4)
pred.svm<-relevel(pred.svm,ref="T")
actual<-relevel(classts$condition,ref="T")
tbl<-table(Predicted=pred.svm,Actual=actual)
confusionMatrix(tbl,positive="T")




# using caret library 
## the data for training data.train , the data for testing : data.test
library(caret)
data.train_t <- t(data.train)
data.test_t <- t(data.test)
conditions_training <- colData(data.trainS4)$condition
conditions_testing <- colData(data.testS4)$condition

## combine training data with condition 
Caret_train <- data.frame(data.train_t, condition = as.factor(conditions_training))

## combine testing data with condition 
Caret_test <- data.frame(data.test_t, condition = as.factor(conditions_testing))

## run the nnet model 
### set up parameters 
tuneGrid <- expand.grid(size = c(1,5), decay = c(0.2,0.1, 0.01))
fitControl <- trainControl(method = "repeatedcv", number = 10,repeats = 20,classProbs = TRUE)
model <- train(condition ~ ., data = Caret_train, method = "nnet",tuneGrid = tuneGrid,trControl=fitControl )

predictions <- predict(model, newdata = Caret_test)

conf_matrix <- confusionMatrix(predictions, Caret_test$condition)


## algorithm : avNNet

fitControl_avNNet <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 20,
  classProbs = TRUE  # For classification tasks, enable class probabilities
)

tuneGrid_avNNet <- expand.grid(
  size = c(1, 3, 5),      # Specify different values for hidden units
  decay = c(0.1, 0.01, 0.001),  # Specify different values for weight decay
  bag = c(0.5, 0.7, 1)      # Specify different bagging fractions
)

model_avNNet <- train(condition ~ ., data = Caret_train ,method = "avNNet",trControl=fitControl_avNNet,tuneGrid =tuneGrid_avNNet )
p_bayesglm <- predict(model_glmboost, newdata = Caret_test)
conf_matrix <- confusionMatrix(p_bayesglm, Caret_test$condition)

## algorithm : pls
ctrl <- trainControl(
  method = "repeatedcv", 
  repeats = 100,
  classProbs = TRUE, 
  summaryFunction = defaultSummary
)

set.seed(123)
plsFit <- train(
  condition ~ .,
  data = Caret_train,
  method = "pls",
  preProc = c("center", "scale"),
  tuneLength = 20,
  trControl = ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(ncomp = c(20))
)

p_plsFit <- predict(plsFit, newdata = Caret_test)
conf_matrix_3 <- confusionMatrix(p_plsFit, Caret_test$condition)

## only 3 class R,S,T


condition_subset_1 <- Caret_train$condition %in% c('R', 'S', 'T')
Caret_train_1 <- Caret_train[condition_subset_1, ]
class_train <- as.character(Caret_train$condition[condition_subset_1])
Caret_train_1 <- Caret_train_1[,-101]
Caret_train_1 <- data.frame(Caret_train_1, condition = as.factor(class_train))

condition_subset_2 <- Caret_test$condition %in% c('R', 'S', 'T')
Caret_test_1 <- Caret_test[condition_subset_2, ]
class_test <-  as.character(Caret_test$condition[condition_subset_2])
Caret_test_1 <- Caret_test_1[,-101]
Caret_test_1 <- data.frame(Caret_test_1, condition = as.factor(class_test))
#condition_subset_3 <-conditions_training %in% c('R', 'S', 'T') 
#conditions_training <- conditions_training[condition_subset_3]

#condition_subset_4 <-conditions_testing %in% c('R', 'S', 'T')
#conditions_testing <- conditions_testing[condition_subset_4]


tuneGrid <- expand.grid(size = c(1,3,5), decay = c(0.3,0.25,0.2,0.15,0.1, 0.01,0.001)) 
fitControl <- trainControl(method = "repeatedcv", number = 10,repeats = 20,classProbs = TRUE)       
model <- train(condition ~ ., data = Caret_train_1, method = "nnet",tuneGrid = tuneGrid,trControl=fitControl )                                         
predictions <- predict(model, newdata = Caret_test_1)
conf_matrix <- confusionMatrix(predictions, Caret_test_1$condition)
png("/bigdata/lerochlab/zli529/RNA_machine_learning/RESULTS/Caret_train_1.png") 
plot(model)
dev.off() 

importance <- varImp(model)
png("/bigdata/lerochlab/zli529/RNA_machine_learning/RESULTS/importance.png")
plot(importance)
dev.off()
