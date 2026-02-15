
rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(meta)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(survivalROC)
  library(forestploter)
})

## Input tables
patient_annotation<- read.table("18_cancer_unicox_8q24_Amp_vs_WT.txt",header = T,row.names=1,sep = '\t',check.names = F)
unigene <- read.table("multicox_2025_2_10_TMB_version2.txt",header = T,row.names=1,sep = '\t',check.names = F)
metagen1<- read.table("uniresult6.txt",header = T,sep = '\t',check.names = F)

######## Forest plot (meta-analysis) ########
################################ Combine meta effect sizes ################################
### Calculate SE
se<-(log(metagen1$H95CI)-log(metagen1$L95CI))/(2*1.96)
metagen1<-cbind(metagen1,se)


m<-metagen(log(metagen1$HR),metagen1$se, sm="HR",studlab=paste(metagen1$gene),comb.fixed=FALSE,data=metagen1)
summary(m)

#settings.meta('JAMA')
settings.meta('RevMan5')
settings.meta('meta4')
forest(m,family='sans',col.diamond.random="red",col.square="black",
       col.square.lines="black",size.square=0.1)


## Publication bias
funnel(m,comb.fixed=T)

## Sensitivity analysis
metainf(m)
forest(metainf(m), comb.fixed=F)

############################# Uni-/multi-variate Cox regression ########################################

### Single-gene univariate Cox regression
td<-read.table("28_LIHC+Sor_Nomogram_clinical_info.txt",header=T) 
library(survival)


pFilter=0.05 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("OS","OS_time") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(td[,4:ncol(td)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(OS_time, OS) ~ td[,i], data = td)#开始逐一循环cox分析
  tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  #if(pvalue<pFilter){ # 这里我们不需要提取所有基因的数据，只需要有意义的gene的结果，所以设置如果pvalue<0.05才提取数据
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                    cbind(id=i,#合并列，是每个基因的统计数据
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
    )
  #}
}


write.table(outResult,file="UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)

###################绘制单基因森林图################################
# 计算标准误差（SE），它在绘图的时候会表示正方形的大小
outResult$H95CI<-as.numeric(outResult$H95CI)
outResult$L95CI<-as.numeric(outResult$L95CI)
outResult$HR<-as.numeric(outResult$HR)
outResult$se <- (log(outResult$H95CI) - log(outResult$HR))/1.96

# 为森林图添加空白列，为了产生一个绘图区间，用于显示CI
outResult$` ` <- paste(rep(" ", 12), collapse = " ")
outResult<-outResult[order(outResult$HR,decreasing = T),]


# forest_theme函数可以对森林图细节进行调整
tm <- forest_theme(
  base_size = 10,        # 设置文本的基础大小
  
  # 设置可信区间的外观
  ci_pch = 15,           # 可信区间点的形状
  ci_col = "blue4",    # 可信区间的边框颜色
  ci_fill = "blue4",      # 可信区间的填充颜色
  ci_alpha = 0.8,        # 可信区间的透明度
  ci_lty = 1,            # 可信区间的线型
  ci_lwd = 1.5,          # 可信区间的线宽
  ci_Theight = 0.2,      # 设置T字在可信区间末端的高度，默认是NULL
  
  # 设置参考线的外观
  refline_lwd = 1,         # 参考线的线宽
  refline_lty = "dashed",  # 参考线的线型
  refline_col = "grey20",  # 参考线的颜色
  
  # 设置垂直线的外观
  vertline_lwd = 1,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
  vertline_lty = "dashed",  # 垂直线的线型
  vertline_col = "grey20",  # 垂直线的颜色
  
  # 设置脚注的字体大小、字体样式和颜色
  footnote_cex = 0.6,            # 脚注字体大小
  footnote_fontface = "italic",  # 脚注字体样式
  footnote_col = "red4"          # 脚注文本的颜色
)


p <- forest(outResult[,c(1, 2, 7, 2, 5)],  # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
            est = outResult$HR,            # 效应值，也就是HR列
            lower = outResult$L95CI,    # 置信区间下限
            upper = outResult$H95CI,    # 置信区间上限
            sizes = outResult$se,          # 黑框框的大小
            ci_column = 3,               # 在第3列（可信区间列）绘制森林图
            ref_line = 1,                # 添加参考线
            arrow_lab = c("Low risk", "High Risk"),   # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
            xlim = c(-1, 10),            # 设置x轴范围
            ticks_at = c(-0.5, 1, 3, 5, 7),   # 在指定位置添加刻度
            theme = tm,                  # 添加自定义主题
            footnote = "This is the demo data. Please feel free to change\nanything you want.")   # 添加脚注信息
p











### 泛癌单因素COX回归分析----
#设置p值的阈值


ESCA <- patient_annotation[patient_annotation$"cancer type abbreviation"== "ESCA",]
OV <- patient_annotation[patient_annotation$"cancer type abbreviation"== "OV",]
UCS <- patient_annotation[patient_annotation$"cancer type abbreviation"== "UCS",]
BRCA <- patient_annotation[patient_annotation$"cancer type abbreviation"== "BRCA",]
LIHC <- patient_annotation[patient_annotation$"cancer type abbreviation"== "LIHC",]
UVM <- patient_annotation[patient_annotation$"cancer type abbreviation"== "UVM",]
STAD <- patient_annotation[patient_annotation$"cancer type abbreviation"== "STAD",]
PAAD <- patient_annotation[patient_annotation$"cancer type abbreviation"== "PAAD",]
HNSC <- patient_annotation[patient_annotation$"cancer type abbreviation"== "HNSC",]
LUAD <- patient_annotation[patient_annotation$"cancer type abbreviation"== "LUAD",]
UCEC <- patient_annotation[patient_annotation$"cancer type abbreviation"== "UCEC",]
LUSC <- patient_annotation[patient_annotation$"cancer type abbreviation"== "LUSC",]
BLCA <- patient_annotation[patient_annotation$"cancer type abbreviation"== "BLCA",]
CHOL <- patient_annotation[patient_annotation$"cancer type abbreviation"== "CHOL",]
READ <- patient_annotation[patient_annotation$"cancer type abbreviation"== "READ",]
PRAD <- patient_annotation[patient_annotation$"cancer type abbreviation"== "PRAD",]
LGG <- patient_annotation[patient_annotation$"cancer type abbreviation"== "LGG",]
SKCM <- patient_annotation[patient_annotation$"cancer type abbreviation"== "SKCM",]


####16种癌合并unicox
unigene2<-rbind(OV,ESCA,UCS, BRCA,LIHC,UVM,STAD,PAAD,HNSC,LUAD,UCEC,LUSC,BLCA,PRAD,
                LGG,SKCM)
unicox <- coxph(Surv(time = os_time, event = os_status) ~ unigene2[,'8q24'], data = unigene2)  
unisum<- summary(unicox)   
pvalue <- round(unisum$coefficients[,5],3) 
#if(pvalue<pfilter)
  uniresult <- cbind(#gene=names(cancer_18types_8q24_df)[i],
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     )




  uniresult


###18种癌分别unicox

cancer_18types_8q24_df <- list(OV,ESCA,UCS, BRCA,LIHC,UVM,STAD,PAAD,HNSC,LUAD,UCEC,LUSC,BLCA,CHOL,
                                  READ,PRAD,LGG,SKCM)
pfilter <- 0.05 
#新建空白list
uniresult <- data.frame()
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
for (i in seq_along(cancer_18types_8q24_df)) {
  unigene1 <- cancer_18types_8q24_df[[i]]
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ unigene1[,'8q24'], data = unigene1)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  #if(pvalue<pfilter)
  
  {
    uniresult <- rbind(uniresult,
                       cbind(gene=names(cancer_18types_8q24_df)[i],
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}



#新建空白数据框
uniresult <- data.frame()
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(unigene[,3:ncol(unigene)])){   
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ unigene[,i], data = unigene)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  #if(pvalue<pfilter)
  { 
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}   
#保存单因素COX回归分析结果

write.table(uniresult,'uniresult7.txt',sep='\t')  





### 多因素COX回归分析----  
multicox <- coxph(Surv(time = OS_time,event = OS) ~ ., data = unigene) 
multisum <- summary(multicox)
#提取所有基因的多因素COX回归分析结果至multiresult对象中
gene <- colnames(unigene)[3:ncol(unigene)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(#gene=gene,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult <- multiresult[multiresult$pvalue<pfilter,]
#保存多因素COX回归分析结果

write.table(multiresult,'multiresult9.txt',sep='\t')  




ggforest(model =multicox,data = unigene )




