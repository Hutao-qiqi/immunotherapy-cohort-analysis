###===Model Prediction===###

rm (list=ls())
gc()

###########################################################################################################
####=======================================1. loading data=============================================####
###########################################################################################################

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

load('data/ios/ios.Rdata');names(ios) #10 ICT RNA-Seq cohorts
load('data/sig/Stem.Sig.Rdata')


###########################################################################################################
####=====================================2. remove batch effect========================================####
###########################################################################################################

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
####===============================3. Model Training and Validation with Stem.Sig======================####
###########################################################################################################


#5 datasets for model training and validation
#"Bruan_RCC_pre_aPD1_combo_tpm" "Mariathasan_UC_pre_aPDL1_combo_tpm"  "Liu_SKCM_pre_aPD1_combo"  "Gide_SKCM_pre_combo" "Riaz_SKCM_pre_aPD1_combo" 

#5 independent datasets for model testing
#"Hugo_SKCM_pre_aPD1" "Van_SKCM_pre_aPD1" "Kim_GC_pre_aPD1" "Zhao_GBM_pre_aPD1" "Snyder_UC_pre_aPDL1" 

data <- combat
grp <- unique(data$batch);grp
combat <- data[!data$batch %in% grp[c(3,7,8,9,10)],] # comat for training and validation set


# 80% as training and 20% as testing
set.seed(13942)
trainIndex <- createDataPartition(combat$response, p = .8,  ## 80% training set; 20% validation set
                                  list = FALSE, 
                                  times = 1)

training <- combat[trainIndex,]
validation  <- combat[-trainIndex,]
test <- data[data$batch %in% grp[c(3,7,8,9,10)],]


training$response <- ifelse(training$response=='0','NR','R') %>% factor(.,levels = c('R','NR'))
validation$response <- ifelse(validation$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
test$response <- ifelse(test$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))

# save(validation,file='results/Machine Learning/validation_set.Rdata')
# save(training,file='results/Machine Learning/training_set.Rdata')
# save(test,file='results/Machine Learning/test_set.Rdata')

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
  

  ##model training with different algorithms
  ls_model <- lapply(method,function(m){
    if(m == 'cancerclass'){ # cancerclass is not avaliable in caret
      pData <- data.frame(class = training$response, sample = rownames(training),row.names = rownames(training))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      Sig.Exp <- t(training[,-1])
      Sig.Exp.train <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
      predictor <- fit(Sig.Exp.train, method = "welch.test") 
      model.tune <- predictor
    } else{ # other algorithms were calculated using R package caret
      
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
  
  ##model validation
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
      auc <- data.frame(ROC=roc$auc, Sens = roc_result$sensitivity, Spec = roc_result$specificity)
      
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
                    sig = Stem.Sig)

# stopCluster(cl) 
print(res[['auc']]) ## 'nb' achieves best performance

class(res$model[[which.max(res$auc$ROC)]]$finalModel) # best model
res$model[[which.max(res$auc$ROC)]]$finalModel$tuneValue ##parameters

# save(res,file='results/Machine Learning/res.Rdata')


###########################################################################################################
####===============================4. Model Testing with Stem.Sig======================================####
###########################################################################################################


sig=Stem.Sig
AUC_Stem.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)],'SKCM'), function(g){
  if(g == 'validation') {
    new = validation
  }else if (g == 'testing'){
    new <- test
  } else if(g == 'SKCM') { #SKCM cohort: Hugo 2016 + Van Allen 2015
    new = test[test$batch %in% grp[c(3,7)], colnames(test) %in% c('response',sig)]
  } else{
    new <- test[test$batch == g, colnames(test) %in% c('response',sig)]
  }
  
  prob <- predict(res[['model']][[which.max(res$auc$ROC)]],new[,-1],type = "prob") # use 'nb' model
  pre <- predict(res[['model']][[which.max(res$auc$ROC)]],new[,-1]) # use 'nb' model
  test_set <- data.frame(obs = new$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
  
  roc <- ROCit::rocit(score = test_set$NR,
                      class = test_set$obs,
                      negref = 'R')
  auc <- roc$AUC
  
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)],'SKCM'))  %>% `colnames<-`('AUC')

AUC_Stem.Sig



###########################################################################################################
####==================5. Compare Stem.Sig with other Signatures in testing datasets====================####
###########################################################################################################

##==loading signatures==##
source('data/sig/sig_others.R')
ImmmunCells.Sig <- readRDS('data/sig/ImSig.rds')
load('data/sig/TcellExc.Sig.Rdata')
load('data/sig/FINAL_Tcell_Exclusion.sig.RData')

ls_sig <- list('Stem.Sig' = Stem.Sig,
               'ImmmunCells.Sig' = ImmmunCells.Sig,
               'TcellExc.Sig' = TcellExc.Sig,
               'IFNG.Sig' = IFNG.Sig,
               'T.cell.inflamed.Sig' = T.cell.inflamed.Sig,
               'Cytotoxic.Sig' = Cytotoxic.Sig,
               'NLRP3.Sig' = NLRP3.Sig,
               'LRRC15.CAF.Sig' = LRRC15.CAF.Sig,
               'PDL1.Sig' = PDL1.Sig,
               'CRMA.Sig' = CRMA.Sig,
               'TRS.Sig' = TRS.Sig
)

#####=============================5.1 pan-cancer signatures==============================================##### 

####5.1.1 IFNG.Sig

#IFNG.Sig use average gene expression as prediction socres
# Reference: Ayers M, Lunceford J, Nebozhyn M, Murphy E, Loboda A, Kaufman DR, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–40. 
# Available from: Available from: https://www.jci.org/articles/view/91190


sig = ls_sig[['IFNG.Sig']];sig


AUC_INFG.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = validation[,colnames(training) %in% c('response', sig)]
  }else if (g == 'testing'){
    new <- test[,colnames(training) %in% c('response', sig)]
  } else{
    new <- test[test$batch == g, colnames(test) %in% c('response',sig)]
  }
  
  roc <- ROCit::rocit(score =rowMeans(new[,2:7]),
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')

AUC_INFG.Sig


####5.1.2 T.cell.inflamed.Sig

#T.cell.inflamed.Sig use average gene expression as prediction socres
# Reference: Ayers M, Lunceford J, Nebozhyn M, Murphy E, Loboda A, Kaufman DR, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–40. 
# Available from: Available from: https://www.jci.org/articles/view/91190

sig = ls_sig[['T.cell.inflamed.Sig']];sig


AUC_T.cell.inflamed.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = validation[,colnames(training) %in% c('response', sig)]
  }else if (g == 'testing'){
    new <- test[,colnames(training) %in% c('response', sig)]
  } else{
    new <- test[test$batch == g, colnames(test) %in% c('response',sig)]
  }
  
  roc <- ROCit::rocit(score =rowMeans(new[,2:7]),
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')

AUC_T.cell.inflamed.Sig

####5.1.3 PDL1.Sig

#PDL1.Sig was the expression of PD-L1
#Reference: Topalian SL, Hodi FS, Brahmer JR, Gettinger SN, Smith DC, McDermott DF, et al. Safety, Activity, and Immune Correlates of Anti–PD-1 Antibody in Cancer. N Engl J Med. 2012;366:2443–54
#Available from: doi: 10.1056/NEJMoa1200690

sig = ls_sig[['PDL1.Sig']];sig


AUC_PDL1.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = validation[,colnames(training) %in% c('response', sig)]
  }else if (g == 'testing'){
    new <- test[,colnames(training) %in% c('response', sig)]
  } else{
    new <- test[test$batch == g, colnames(test) %in% c('response',sig)]
  }
  roc <- ROCit::rocit(score =new$PDCD1,
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')

AUC_PDL1.Sig



####5.1.4 LRRC15.CAF.Sig

#LRRC15.CAF.Sig use eigenWeightedMean method in multiGSEA package to calculate the prediction scores
#Reference: Dominguez CX, Müller S, Keerthivasan S, Koeppen H, Hung J, Gierke S, et al. Single-Cell RNA Sequencing Reveals Stromal Evolution into LRRC15 + Myofibroblasts as a Determinant of Patient Response to Cancer Immunotherapy. Cancer Discov [Internet]. 2020;10:232–53. 
#Available from: http://cancerdiscovery.aacrjournals.org/lookup/doi/10.1158/2159-8290.CD-19-0644

# BiocManager::install("lianos/multiGSEA.shiny@develop")
library(multiGSEA) #please install multiGSEA using local package in "renv/local/lianos-multiGSEA-prebioc-61-g85719e2.tar.gz" the "eigenWeightedMean" method is implanted in this old version of multiGSEA
sig = ls_sig[['LRRC15.CAF.Sig']];sig


AUC_LRRC15.CAF.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = validation[,colnames(training) %in% c('response', sig)]
  }else if (g == 'testing'){
    new <- test[,colnames(training) %in% c('response', sig)]
  } else{
    new <- test[test$batch == g, colnames(test) %in% c('response',sig)]
  }
  expr <- as.data.frame(t(new[,-1]))
  
  scores <- eigenWeightedMean(expr)$score
  
  roc <- ROCit::rocit(score = scores,
                      class = new$response,
                      negref = 'R')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')

AUC_LRRC15.CAF.Sig

####5.1.5 NLRP3.Sig
#NLRP3.Sig was calculated using ssGSEA
#Reference: Ju M, Bi J, Wei Q, et al. Pan-cancer analysis of NLRP3 inflammasome with potential implications in prognosis and immunotherapy in human cancer. Brief Bioinform 2020; 00: 1–16.
#Available from: DOI: 10.1093/bib/bbaa345

library(GSVA)
library(GSEABase)
sig <- NLRP3.Sig

gmt <- getGmt('data/gmt/NLRP3.Sig.gmt')
all(duplicated(names(gmt))) #no duplicated ids


AUC_NLRP3.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = data[rownames(validation), ]
  } else if (g == 'testing'){
    new <- data[data$batch %in% grp[c(3,7,8,9,10)], ]
  } else{
    new <- data[data$batch == g,]
  }
  
  new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
  
  expr <- as.matrix(t(new[,-c(1:8)]))
  
  gsva <- gsva(expr,gmt,method='ssgsea',parallel.sz=1)
  
  score <- as.numeric(gsva)
  
  roc <- ROCit::rocit(score = score,
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')


AUC_NLRP3.Sig


####5.1.6 Cytotoxic.Sig

#Cytotoxic.Sig use use geometric mean of gene expression as prediction socres
#Reference:Rooney MS, Shukla SA, Wu CJ, Getz G, Hacohen N. Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell [Internet]. Elsevier Inc.; 2015;160:48–61. 
#Available from: http://dx.doi.org/10.1016/j.cell.2014.12.033

sig = ls_sig[['Cytotoxic.Sig']];sig


Cytotoxic_combat<- ios_bc[,9:ncol(ios_bc)]
#log2 changes should be removed for the following geometric mean calculation
Cytotoxic_combat <- Cytotoxic_combat^2
Cytotoxic_combat <- cbind(ios_bc[,1:8],Cytotoxic_combat)
rownames(Cytotoxic_combat) <- ios_bc$ID
Cytotoxic_combat <-  Cytotoxic_combat[!is.na(Cytotoxic_combat$response),] ##remove patient with unknown response status

Cytotoxic_data <- Cytotoxic_combat

Cytotoxic_combat <- Cytotoxic_data[!Cytotoxic_data$batch %in% grp[c(3,7,8,9,10)],]



AUC_Cytotoxic.Sig <- lapply(c('validation', 'testing',grp[c(3,7,8,9,10)]), function(g){
  if(g == 'validation') {
    new = Cytotoxic_data[rownames(validation) ,colnames(Cytotoxic_data) %in% c('response',sig)]
  } else if (g == 'testing'){
    new <- Cytotoxic_data[Cytotoxic_data$batch %in% grp[c(3,7,8,9,10)], colnames(Cytotoxic_data) %in% c('response',sig)]
  } else{
    new <- Cytotoxic_data[Cytotoxic_data$batch == g, colnames(Cytotoxic_data) %in% c('response',sig)]
  }
  
  new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
  
  roc <- ROCit::rocit(score = compositions::geometricmeanRow(new[,2:ncol(new)]),
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}) %>% do.call(rbind,.) %>% `rownames<-`(c('validation', 'testing',grp[c(3,7,8,9,10)]))  %>% `colnames<-`('AUC')

AUC_Cytotoxic.Sig






#####=============================5.2 pan-cancer signatures==============================================##### 



####5.2.1 ImmuneCells.Sig

##ImmuneCells.Sig Model was constructed using "cancerclass" according to the original article
#Reference:Xiong D, Wang Y, You M. A gene expression signature of TREM2hi macrophages and γδ T cells predicts immunotherapy response. Nat Commun. 2020;11:1–12. 
#Available from: http://dx.doi.org/10.1038/s41467-020-18546-x

sig = ls_sig[['ImmmunCells.Sig']]

res <- CompareModel(training = training,
                    validation = validation,
                    method = 'cancerclass',
                    sig = sig)

getImmmunCells.Sig <- function(){
  new = test[test$batch %in% grp[c(3,7)], colnames(test) %in% c('response',sig)]#SKCM cohort: Hugo 2016 + Van Allen 2015
 
  pData <- data.frame(class =new$response, sample = rownames(new),row.names = rownames(new))
  phenoData <- new("AnnotatedDataFrame",data=pData)
  Sig.Exp <- t(new[,-1])
  Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
  prediction <- predict(res$model[[1]], Sig.Exp.test,positive="NR", ngenes=nrow(Sig.Exp), dist = "cor")
  roc <- ROCit::rocit(score =as.numeric(prediction@prediction[,'z']),
                      class = new$response,
                      negref = 'R')
  auc <- roc$AUC
}

AUC_ImmmunCells.Sig <- getImmmunCells.Sig()
AUC_ImmmunCells.Sig

####5.2.2 TcellExc.Sig

# TcellExc.Sig use overall expression as prediction scores 
# Reference: Jerby-Arnon L, Shah P, Cuoco MS, Rodman C, Su M-J, Melms JC, et al. A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade. Cell.2018;175:984-997.e24.
# Available from: https://linkinghub.elsevier.com/retrieve/pii/S0092867418311784

source('data/sig/IMPRES/ImmRes_source.R') ## 'ImmRes_OE.R' was downloaded from https://github.com/livnatje/ImmuneResistance
library(recipes)

sig = ls_sig[['TcellExc.Sig']]

gene.sign = list(TcellExc.Sig = TcellExc.Sig)

prd <- function(df){
  r <- list(tpm = t(df[,-1]), genes = colnames(df[,-1]))
  OE <- get.OE.bulk(r = r, 
                    gene.sign = gene.sign)   
  
  roc <- ROCit::rocit(score = OE[,'TcellExc.Sig'],
                      class = df$response,
                      negref = 'R')
  auc <- roc$AUC
  print(auc)
}


getTcellExc.Sig <-  function(g){
    new <- test[test$batch %in% grp[c(3,7)], colnames(test) %in% c('response',sig)] #SKCM cohort: Hugo 2016 + Van Allen 2015
  
  return(prd(new))
}
AUC_TcellExc.Sig <- getTcellExc.Sig()
AUC_TcellExc.Sig

####5.2.3 CRMA.Sig

#CRMA.Sig use geometric mean of gene expression as prediction socres
#Shukla SA, Bachireddy P, Schilling B, Galonska C, Zhan Q, Bango C, et al. Cancer-Germline Antigen Expression Discriminates Clinical Outcome to CTLA-4 Blockade. Cell [Internet]. Elsevier Inc.; 2018;173:624-633.e8. 
#Available from: https://doi.org/10.1016/j.cell.2018.03.026


sig = ls_sig[['CRMA.Sig']];sig


CRMA_combat<- ios_bc[,9:ncol(ios_bc)]
#log2 changes should be removed for the following geometric mean calculation
CRMA_combat <- CRMA_combat^2
CRMA_combat <- cbind(ios_bc[,1:8],CRMA_combat)
rownames(CRMA_combat) <- ios_bc$ID
CRMA_combat <-  CRMA_combat[!is.na(CRMA_combat$response),] ##remove patient with unknown response status

CRMA_data <- CRMA_combat

CRMA_combat <- CRMA_data[!CRMA_data$batch %in% grp[c(3,7,8,9,10)],]

res <- list()
p <- list()
AUC <- list()



getCRMA.Sig <- function(){
  new <- CRMA_data[CRMA_data$batch %in% grp[c(3,7)], colnames(CRMA_data) %in% c('response',sig)]#SKCM cohort: Hugo 2016 + Van Allen 2015
  new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
  roc <- ROCit::rocit(score = compositions::geometricmeanRow(new[,2:ncol(new)]),
                      class = new$response,
                      negref = 'NR')
  auc <- roc$AUC
  return(auc)
}
AUC_CRMA.Sig <- getCRMA.Sig()
AUC_CRMA.Sig




####5.2.4 IMPRES.Sig

#IMPRES.Sig uses the comparision scores between 15 gene pairs as prediction socres. 
#Reference:Auslander N, Zhang G, Lee JS, Frederick DT, Miao B, Moll T, et al. Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma. Nat Med [Internet]. Springer US; 2018;24:1545–9. 
#Available from: http://dx.doi.org/10.1038/s41591-018-0157-9


gp <- data.table::fread('data/sig/IMPRES.csv')

gp <- as.data.frame(gp)
gp1 <- gp
gp2 <- gp

for(i in c('A','B')){
  gp2[,i] <- str_replace(gp2[,i],'C10orf54','VSIR') # C10orf54 and VSIR are interchangable 
}
g1 <- unique(c(gp1$A,gp1$B))
g2 <- unique(c(gp2$A,gp2$B))

IMPRES <- function(pt){
  pt <- as.data.frame(t(pt))
  IMPRES.scores = 0
  for(kk in 1:nrow(g)){
    A = g[kk,'A']
    B = g[kk,'B']
    try(if(pt[,A] > pt[,B]) {IMPRES.scores = IMPRES.scores + 1})
  }
  return(IMPRES.scores)
}


##number of avalible IMPRES.Sig genes for each cohort

for(i in names(ios)){
  io.g <- colnames(ios[[i]])
  if('C10orf54' %in% io.g ){
    print(paste(i, sum(g1 %in% io.g),'C10orf54'))
    ios[[i]]$VSIR = ios[[i]]$C10orf54 #add 'VISR' in cohort with "C10orf54"
  } else if ('VSIR' %in% io.g) {
    print(paste(i, sum(g2 %in% io.g),'VSIR'))
  } else {
    print(paste(i, sum(g1 %in% io.g)))
  }
}

ls_IMPRES <- list()
IMPRES_ios <- ios #remove Kim_GC and Synder_UC. Validation cohort was not affected.


for(i in names(IMPRES_ios)){
  g3 <- g2[g2 %in%  colnames(IMPRES_ios[[i]])]
  io <- IMPRES_ios[[i]][,g3]
  g <- as.data.frame(gp2)
  g <- g[g$A %in% g3 & g$B %in% g3 ,]
  s <- apply(io,1,IMPRES)
  s <- data.frame(IMPRES = s)
  s$ID <- IMPRES_ios[[i]]$ID
  ls_IMPRES[[i]] = s
}


getIMPRES.Sig <- function(){
  
    IMPRES_scores <- ls_IMPRES[grp[c(3,7)]] %>% do.call(rbind,.) #SKCM cohort: Hugo 2016 + Van Allen 2015

  
  IMPRES_response <- data[data$ID %in% IMPRES_scores$ID,'response'] 
  IMPRES_response <-  ifelse(IMPRES_response == 0,'NR','R') %>% factor(.,levels = c('R','NR'))
  
  roc <- ROCit::rocit(IMPRES_scores$IMPRES,IMPRES_response, negref = 'NR')
  
  auc <- roc$AUC
  return(auc)
}

AUC_IMPRES.Sig <- getIMPRES.Sig()
AUC_IMPRES.Sig





####5.2.5 IPRES.Sig

#IPRES.Sig used mean z-scores of GSVA scores calculated from 73 IPRES datasets. 
#Reference:Hugo W, Zaretsky JM, Sun L, Johnson DB, Ribas A, Lo RS, et al. Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma Article Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell [Internet]. Elsevier Inc.; 2016;1–10. 
#Available from: http://dx.doi.org/10.1016/j.cell.2016.02.065

#Original article does not provide the get sets. We manually collected datasets from supplemtary of original article and http://www.gsea-msigdb.org 

###=== GSVA ===###
library(GSVA)
library(GSEABase)

gmt <- getGmt('data/gmt/IPRES.symbols.gmt')
all(duplicated(names(gmt))) #no duplicated ids


getIPRES.Sig <- function(){
  new <- data[data$batch %in% grp[c(3,7)], ] #SKCM cohort: Hugo 2016 + Van Allen 2015
  
  new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
  
  expr <- as.matrix(t(new[,-c(1:8)]))
  
  gsva <- gsva(expr,gmt,method='gsva',parallel.sz=60)
  
  score <- scale(t(gsva))
  
  score <- rowMeans(score)
  
  roc <- ROCit::rocit(score = score,
                      class = new$response,
                      negref = 'R')
  auc <- roc$AUC
  return(auc)
}

AUC_IPRES.Sig <- getIPRES.Sig()
AUC_IPRES.Sig


####5.2.6 TRS.Sig
#TRS.Sig was calculated using GSVA
#Reference: Yan M, Hu J, Ping Y, Xu L, Liao G, Jiang Z, Pang B, Sun S, Zhang Y, Xiao Y, Li X. Single-Cell Transcriptomic Analysis Reveals a Tumor-Reactive T Cell Signature Associated With Clinical Outcome and Immunotherapy Response In Melanoma. Front Immunol. 2021 Nov 5;12:758288. 
#Available from: doi: 10.3389/fimmu.2021.758288.
library(GSVA)
library(GSEABase)
sig <- TRS.Sig

gmt <- getGmt('data/gmt/TRS.Sig.gmt')
all(duplicated(names(gmt))) #no duplicated ids


getTRS.Sig <- function(){
new <- data[data$batch %in% grp[c(3,7)], ] #SKCM cohort: Hugo 2016 + Van Allen 2015
new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
expr <- as.matrix(t(new[,-c(1:8)]))
gsva <- gsva(expr,gmt,method='gsva',parallel.sz=1)
score <- as.numeric(gsva)
roc <- ROCit::rocit(score = score,
                      class = new$response,
                      negref = 'NR')
auc <- roc$AUC
  return(auc)
}

AUC_TRS.Sig <- getTRS.Sig()

AUC_TRS.Sig

####5.2.7 IMS.Sig
#IMS.Sig was calculated as ratio of (IMS scores / IFN-γ scores)
#IMS scores = Average gene expression of 18 immunosupression genes
#IFN-γ scores = Average gene expression of 6 INF-γ genes
#Reference: Cui, C., Xu, C., Yang, W., Chi, Z., Sheng, X., Si, L., ... & Kong, Y. (2021). Ratio of the interferon-γ signature to the immunosuppression signature predicts anti-PD-1 therapy response in melanoma. NPJ genomic medicine, 6(1), 1-12.
#Available from: https://www.nature.com/articles/s41525-021-00169-w
sig1 <- IMS.Sig
sig2 <- IFNG.Sig


getIMS.Sig <- function(){
  new <- data[data$batch %in% grp[c(3,7)], ] #SKCM cohort: Hugo 2016 + Van Allen 2015
  new$response <- ifelse(new$response=='0','NR','R')%>% factor(.,levels = c('R','NR'))
  new1 <- new[,colnames(new) %in% c('response', sig1)]
  new2 <- new[,colnames(new) %in% c('response', sig2)]
  
  score1 <- rowMeans(new1[,-1])
  score2 <- rowMeans(new2[,-1])
  
  score <- score1/score2  
  
  roc <- ROCit::rocit(score = score,
                      class = new$response,
                      negref = 'R')
  auc <- roc$AUC
  return(auc)
}

AUC_IMS.Sig <- getIMS.Sig()

AUC_IMS.Sig






####5.12 compare AUC

AUC <- list(pancancer = list(),SKCM = list())

AUC[['pancancer']] <- list(Stem.Sig = data.frame(AUC = AUC_Stem.Sig[-8,]),
                           INFG.Sig = AUC_INFG.Sig,
                           T.cell.inflamed.Sig = AUC_T.cell.inflamed.Sig,
                           Cytotoxic.Sig = AUC_Cytotoxic.Sig,
                           PDL1.Sig = AUC_PDL1.Sig,
                           LRRC15.CAF.Sig = AUC_LRRC15.CAF.Sig,
                           NLRP3.Sig = AUC_NLRP3.Sig
                          )
AUC[['pancancer']] <- lapply(AUC[['pancancer']],as.data.frame)

AUC[['SKCM']] <- list(Stem.Sig = AUC_Stem.Sig[8,],          
                      IMPRES.Sig = AUC_IMPRES.Sig,
                      CRMA.Sig = AUC_CRMA.Sig,
                      ImmmunCells.Sig = AUC_ImmmunCells.Sig,
                      TcellExc.Sig =  AUC_TcellExc.Sig,
                      IPRES.Sig = AUC_IPRES.Sig,
                      TRS.Sig  = AUC_TRS.Sig,
                      IMS.Sig = AUC_IMS.Sig)


# save(AUC,file='results/Machine Learning/AUC.Rdata')












