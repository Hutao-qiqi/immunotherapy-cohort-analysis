###===Figure 4===###
rm (list=ls())
gc()

##===loading data===###
library(ROCit)
library(caret)
library(dplyr)
library(survminer)
library(survival)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ROCit)
library(ggsci)
library(cancerclass)

load('results/Machine Learning/validation_set.Rdata')
load('results/Machine Learning/training_set.Rdata')
load('results/Machine Learning/test_set.Rdata')
load('results/Machine Learning/res.Rdata')
load('data/sig/Stem.Sig.Rdata')
load('results/Machine Learning/AUC.Rdata')


###===Figure 4A===###

# Flow chart


###===Figure 4B===###


models <-  c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
roc <- lapply(1:7,function(i){
  if (!i == 7) {
    prob <- predict(res[['model']][[i]],validation[,-1],type = "prob") # use 'nb' model
    pre <- predict(res[['model']][[i]],validation[,-1]) # use 'nb' model
    test_set <- data.frame(obs = validation$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
    roc <- rocit(score = test_set$NR,
                 class = test_set$obs,
                 negref = 'R')
  } else{
    vali <- validation[,colnames(validation) %in% c('response',Stem.Sig)]
    pData <- data.frame(class = vali$response, sample = rownames(vali),row.names = rownames(vali))
    phenoData <- new("AnnotatedDataFrame",data=pData)
    Sig.Exp <- t(vali[,-1])
    Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
    prediction <- predict(res[['model']][[i]], Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
    roc <- rocit(as.numeric(prediction@prediction[,'z']),
                 prediction@prediction[,'class_membership'])
  }
})

colset <- pal_lancet("lanonc", alpha = 0.7)(7)
plot(roc[[1]], col = colset[1], 
     legend = FALSE, YIndex = F)

for(i in 2:7){
  lines(roc[[i]]$TPR~roc[[i]]$FPR, 
        col = colset[i], lwd = 2)
}

legend("bottomright", 
       col = colset[1:7],
       paste(models, 'AUC', round(res[['auc']]$ROC,2)), 
       lwd = 2)






###===Figure 4C===###

#ROC plot

rocplot <- function(data){
prob <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1],type = "prob") # use 'nb' model
pre <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1]) # use 'nb' model
test_set <- data.frame(obs = data$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
roc <- ROCit::rocit(score = test_set$NR,
                    class = test_set$obs,
                    negref = 'R')
plot(roc,legend=F,YIndex = F)
title(deparse(substitute(data)))
text(x=0.5,y=0.6,labels = paste0("AUC: ",round(roc$AUC,2), " (95%CI: ",round(as.numeric(ciAUC(roc)[5]),2),'-',round(as.numeric(ciAUC(roc)[6]),2),")"))

}


rocplot(validation) # validation
rocplot(test) # testing


###===Figure 4D===###

#Survival Plot

survplot <- function(data){
  if(all(is.na(data$OS))){return(NA)} ## OS data is unavaliable in "Kim_GC_pre_aPD1"
  
  model.name = deparse(substitute(data))
  pre <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1]) # use 'nb' model

  data$pre <- pre
  data <- data[!is.na(data$OS),] #remove patient without OS data
  data$OS <- as.numeric(data$OS)
  data$status <- as.numeric(data$status)
  data$pre <- factor(data$pre,levels = c('R','NR'))

  surv <- Surv(data$OS, data$status)
  fit <- surv_fit(surv~pre,data =data)
  medianOS <- surv_median(fit)
  survdiff <- survdiff(surv~pre,data = data)
  cox <- coxph(surv~pre, data=data) #


  pic <- ggsurvplot(fit, data = data,
                    title= model.name,
                    legend.labs = c("Low risk","High risk"),
                    break.x.by = 4,ylab="OS(probability)", xlab = " Time (Months)",
                    palette = c("#E7B800","#2E9FDF"),
                    censor.shape = 3,censor.size = 1,
                    legend.title='',
                    surv.median.line = "hv"
  )

  pic$plot <- pic$plot + ggplot2::annotate("text",x = 20, y = 0.68, size =3,
                                           label = paste("HR :",format(round(summary(cox)$conf.int[1],2),nsmall=2))) +
    ggplot2::annotate("text",x = 20, y = 0.62, size = 3,
                      label = paste("(","95%CI: ", format(round(summary(cox)$conf.int[1,3],2),nsmall = 2),"-",format(round(summary(cox)$conf.int[1,4],2),nsmall = 2),")",sep = ""))+
    ggplot2::annotate("text",x = 20, y = 0.95, size =3,
                      label = paste("Median OS"))+
    ggplot2::annotate("text",x = 20, y = 0.89, size =3,
                      label = paste(round(medianOS$median[1],2),"months")) +
    ggplot2::annotate("text",x = 20, y = 0.82, size =3,
                      label = paste(round(medianOS$median[2],2),"months"))+
    ggplot2::annotate("text",x = 20, y = 0.56, size =3,
                      label = paste("Log-rank p:",round(1-pchisq(survdiff$chisq,1),4)))+
    ggplot2::theme(legend.text=element_text(size=8,face = "bold"))+
    theme(axis.line.x = element_line(size=0.8),
          axis.line.y = element_line(size=0.8),
    )+
    font("xy.title",face = "bold")+
    scale_x_continuous(breaks = seq(0, max(data$OS,na.rm=T), by = 12))

  return(pic$plot)
}



survplot(validation)
survplot(test)


