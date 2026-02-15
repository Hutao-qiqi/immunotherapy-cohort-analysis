###===Figure S2===###

rm (list=ls())
gc()

##===loading data===###
library(ROCit)
library(caret)
library(dplyr)
library(survminer)
library(survival)
library(reshape2)
library(ggpubr)
library(ggplot2)

setwd("~/Projects/GM/")

load('results/Machine Learning/test_set.Rdata')
load('results/Machine Learning/res.Rdata')
load('data/ios/ls_meta.Rdata')


###===Fig.S2 A===###

#ROC plot 
rocplot <- function(data){
  prob <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1],type = "prob") # use 'nb' model
  pre <- predict(res[['model']][[which.max(res$auc$ROC)]],data[,-1]) # use 'nb' model
  test_set <- data.frame(obs = data$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
  roc <- ROCit::rocit(score = test_set$NR,
                      class = test_set$obs,
                      negref = 'R')
  plot(roc,legend=F,YIndex = F)
  title(unique(data$batch))
  text(x=0.5,y=0.6,labels = paste0("AUC: ",round(roc$AUC,2), " (95%CI: ",round(as.numeric(ciAUC(roc)[5]),2),'-',round(as.numeric(ciAUC(roc)[6]),2),")"))
  
}


cohort <- unique(test$batch)
cohort

par(mfrow = c(2, 3))
rocplot(test[test$batch == cohort[1],]) #Hugo_SKCM_pre_aPD1
rocplot(test[test$batch == cohort[2],]) #Van_SKCM_pre_aPD1
rocplot(test[test$batch == cohort[3],]) #Kim_GC_pre_aPD1
rocplot(test[test$batch == cohort[4],]) #Zhao_GBM_pre_aPD1
rocplot(test[test$batch == cohort[5],]) #Snyder_UC_pre_aPDL1



###===Fig.S2 B===###

#Survival Plot



survplot <- function(data){
  g = unique(data$batch) #cohort name
  
  if(all(is.na(data$OS))){
    return(NA)
  } else { ## OS data is unavaliable in "Kim_GC_pre_aPD1"
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
    
    ## adjusted available confounding factors for

    if(cohort %in% names(ls_meta)){
    data <- dplyr::left_join(data,ls_meta[[g]],by='ID')
    }
    
    #adjusting avaliable confoudning factors for cox surviaval analysis
    if(g == 'Hugo_SKCM_pre_aPD1'){
    cox <- coxph(surv~pre+Age+Sex+Purity+TMB,data = data)
    } else if (g == 'Van_SKCM_pre_aPD1'){
    cox <- coxph(surv~pre+Age+Sex+TMB,data = data)
    } else if (g == 'Zhao_GBM_pre_aPD1'){
    cox <- coxph(surv~pre+Age+Sex,data = data)
    } else if (g == 'Snyder_UC_pre_aPDL1'){
    cox <- coxph(surv~pre+Age+Sex,data = data)
    } else {
      
    }
  
    
    
    pic <- ggsurvplot(fit, data = data,
                      title= g,
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
      ggplot2::annotate("text",x = 20, y = 0.50, size =3,
                        label = paste("adjust p:",format(round(summary(cox)$coefficients[1,5],2),nsmall=2)))+
      theme(axis.line.x = element_line(size=0.8),
            axis.line.y = element_line(size=0.8),
      )+
      font("xy.title",face = "bold")+
      scale_x_continuous(breaks = seq(0, max(data$OS,na.rm=T), by = 12))
    
    return(pic$plot)
  }
}

ls_survplot <- list()

for(i in cohort[c(1,2,4,5)]){
ls_survplot[[i]] <- survplot(test[test$batch == i,]) #Hugo_SKCM_pre_aPD1

}


ggarrange(plotlist= ls_survplot)





