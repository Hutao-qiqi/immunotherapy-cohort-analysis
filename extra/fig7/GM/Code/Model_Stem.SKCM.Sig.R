
###===Analysis of melanoma specific signatures===###

###===Model Prediction===###

rm (list=ls())
gc()
###===loading data===###
library(stringr)
library(gridExtra)
library(future)
library(sva)
library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)


load('data/ios/ios.Rdata');names(ios)
load('data/sig/Stem.SKCM.Sig.Rdata')

ios = ios[grepl('SKCM',names(ios))]
###===remove batch effect===### 

for(i in 1:length(ios)){
  if(i == 1){
    cm = colnames(ios[[1]])
  } else{
    cm = intersect(cm,colnames(ios[[i]]))
  }
}
ios_bc <- lapply(ios,function(z){ return(z[,cm])}) %>% do.call(rbind,.)
ios_bc <- cbind(rownames(ios_bc),ios_bc)
colnames(ios_bc)[1] <- 'batch'
ios_bc$batch <- str_split_fixed(ios_bc$batch,'\\.',2)[,1]
rownames(ios_bc) <- ios_bc$ID
edata <- t(ios_bc[,9:ncol(ios_bc)])
combat <- ComBat(dat = edata, 
                 batch = ios_bc$batch
)
combat <- as.data.frame(t(combat))
combat <- cbind(ios_bc[,1:8],combat)
combat <- combat[!is.na(combat$response),] ##remove patient with unknown response status

###########################################################################################################
####===============================Model Training with Stem.Sig========================================####
###########################################################################################################


#5 datasets for model training and validation
#"Liu_SKCM_pre_aPD1_combo"  "Gide_SKCM_pre_combo" "Riaz_SKCM_pre_aPD1_combo" 

#5 independent datasets for model testing
#"Hugo_SKCM_pre_aPD1" "Van_SKCM_pre_aPD1"

data <- combat
grp <- unique(data$batch);grp
combat <- data[!data$batch %in% grp[c(1,5)],] # comat for training and validation set


# 80% as training and 20% as testing
set.seed(13942)
trainIndex <- createDataPartition(combat$response, p = .8,  ## 80% training set; 20% validation set
                                  list = FALSE, 
                                  times = 1)

training <- combat[trainIndex,]
validation  <- combat[-trainIndex,]
test <- data[data$batch %in% grp[c(1,5)],]


training$response <- ifelse(training$response=='0','NR','R') %>% factor(.,levels = c('R','NR'))
validation$response <- ifelse(validation$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
test$response <- ifelse(test$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))


#Training set: train model and tune parameters (10 times repeated, 5 folds cross-validation, 10 tuneLenth for each parameter)
#Validation set: Compare model performance and pick the best one as the final model
#Test set: Independent evaluation of the performance of final model

CompareModel <- function(training, validation, method,sig){
  
  training <- training[,colnames(training) %in% c('response', sig)]
  validation  <- validation[,colnames(validation) %in% c('response',sig)]
  
  #7 models adpoted in this study as followings: 
  #'nb': navie bayes
  #'svmRadialWeights': Support Vector Machines with Class Weights
  #'rf': random forest
  #'kknn': k-Nearest Neighbors
  #'adaboost':AdaBoost Classification Trees
  #'LogitBoost':Boosted Logistic Regressions
  #'cancerclass': cancerclass
  
  #Grid search for parameter tuning
  Grid <- list( nb = expand.grid(fL =  c(0,0.5,1,1.5,2.0), usekernel = TRUE, adjust = c(0.5,0.75,1,1.25,1.5)),
                svmRadialWeights = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01 ,0.05),C = c( 1 ,3 ,5 ,10 ,20), Weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
                rf = expand.grid(mtry = c(2,42,83,124,165,205,246,287,328,369)),
                kknn = expand.grid(kmax = c(5,7,9,11,13), distance = 2 , kernel = 'optimal'),
                adaboost = expand.grid(nIter = c(50,100,150,200,250) ,method= c('Adaboost.M1','Real adaboost')),
                LogitBoost = expand.grid(nIter = c(11,21,31,41,51,61,71,81,91,101) )
  )
  
  
  TuneLength =  list( nb = nrow(Grid[['nb']]),
                      svmRadialWeights = nrow(Grid[['svmRadialWeights']]) ,
                      rf = nrow(Grid[['rf']]),
                      kknn =nrow(Grid[['kknn']]) ,
                      adaboost = nrow(Grid[['adaboost']]),
                      LogitBoost =  nrow(Grid[['LogitBoost']])
  )
  
  
  ls_model <- lapply(method,function(m){
    if(m == 'cancerclass'){ #cancerclass is not avaliable in caret
      pData <- data.frame(class = training$response, sample = rownames(training),row.names = rownames(training))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      Sig.Exp <- t(training[,-1])
      Sig.Exp.train <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
      predictor <- fit(Sig.Exp.train, method = "welch.test") 
      model.tune <- predictor
    } else{
      
      f = 5  # f folds resampling
      r = 10 # r repeats
      n = f*r
      
      # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
      seeds <- vector(mode = "list", length = n + 1)
      #the number of tuning parameter
      for(i in 1:n) seeds[[i]] <- sample.int(n=1000, TuneLength[[m]])
      
      #for the last model
      seeds[[n+1]]<-sample.int(1000, 1)
      
      
      ctrl <- trainControl(method="repeatedcv",
                           number = f, ## 5-folds cv
                           summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                           classProbs=TRUE,
                           repeats = r, ## 10-repeats cv,
                           seeds = seeds
      )
      
      
      
      model.tune <- train(response ~ .,
                          data = training,
                          method = m,
                          metric="ROC",
                          trControl=ctrl,
                          tuneGrid = Grid[[m]]
      )
      
    }
    print(m)
    return(model.tune)
  }
  )
  
  
  auc <- lapply(ls_model,function(model.tune){
    if(class(model.tune) == 'predictor'){
      pData <- data.frame(class = validation$response, sample = rownames(validation),row.names = rownames(validation))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      Sig.Exp <- t(validation[,-1])
      Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
      prediction <- predict(model.tune, Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
      roc <- roc(response  = prediction@prediction[,'class_membership'],
                 predictor = as.numeric(prediction@prediction[,'z'])
      )
      roc_result <- coords(roc, "best")
      auc <- data.frame(ROC=as.numeric(roc$auc), Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
    }else {
      prob <- predict(model.tune,validation[,-1],type = "prob")
      pre <- predict(model.tune,validation[,-1])
      test_set <- data.frame(obs = validation$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
      auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
    }
    
    return(auc)
  }) %>% do.call(rbind,.)
  
  rownames(auc) <- method
  
  res <- list()
  
  res[['model']] <- ls_model
  res[['auc']] <- auc
  
  
  return(res)
  
}

#parallel processing
cl <- makePSOCKcluster(60)
registerDoParallel(cl)

res <- CompareModel(training = training,
                    validation = validation,
                    method = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                    sig = Stem.SKCM.Sig)

stopCluster(cl)
print(res[['auc']]) 

class(res$model[[which.max(res$auc$ROC)]]$finalModel) # best model
res$model[[which.max(res$auc$ROC)]]$finalModel$tuneValue ##parameters


rocplot <- function(data){
  cohort.name <- deparse(substitute(data))
  prob <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1],type = "prob") # use 'nb' model
  pre <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1]) # use 'nb' model
  test_set <- data.frame(obs = data$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
  roc <- ROCit::rocit(score = test_set$NR,
                      class = test_set$obs,
                      negref = 'R')
  plot(roc,legend=F,YIndex = F)
  title(cohort.name)
  text(x=0.5,y=0.6,labels = paste0("AUC: ",round(roc$AUC,2), " (95%CI: ",round(as.numeric(ciAUC(roc)[5]),2),'-',round(as.numeric(ciAUC(roc)[6]),2),")"))
  
}

rocplot(validation) # validation
rocplot(test) # testing







