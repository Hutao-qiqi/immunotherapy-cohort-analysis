###### Trend line plot ######

rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

# Input file (tab-delimited)
input_file <- "TIL_map_PRAD.TMPRSS2-ERG Fusion.txt"
TIL <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

#exp_data <- read.table("TARGET-OS.htseq_fpkm-uq.txt",header=T,row.names=1,sep = '\t',check.names = F)

#aggregate(x=TIL$til_percentage, by=list(TIL$"8q24",TIL$"cancer type abbreviation"),mean)

TIL_sum<-aggregate(x=TIL$"Global_Pattern", by=list(TIL$"8q24",TIL$"cancer type abbreviation",TIL$"Global_Pattern"),length)



TIL_mean<-aggregate(x=TIL$til_percentage, by=list(TIL$"8q24",TIL$"cancer type abbreviation"),mean)
#colnames(TIL_mean)[1]<-"E8q24"
#colnames(TIL_mean)[2]<-"cancer_type"
#colnames(TIL_mean)[3]<-"TIL"


TIL_sd<-aggregate(x=TIL$til_percentage, by=list(TIL$"8q24",TIL$"cancer type abbreviation"),sd)
TIL_mean_sd<-cbind(TIL_sd,TIL_mean[,3])
colnames(TIL_mean_sd)[1]<-"E8q24"
colnames(TIL_mean_sd)[2]<-"cancer_type"
colnames(TIL_mean_sd)[3]<-"TIL_sd"
colnames(TIL_mean_sd)[4]<-"TIL_mean"


#fix(TIL_mean)


TIL_mean <- TIL_mean%>%mutate(
  E8q24=factor(E8q24),
  cancer_type=factor(cancer_type))



TIL_mean <- TIL_mean%>%mutate(
E8q24=factor(E8q24,labels=c("WT","gain","Amp")),
cancer_type=factor(cancer_type,labels = c("UCEC","STAD","LUSC","LUAD")))


# Example: filter with dplyr if needed
ggplot(data=TIL_mean_sd,mapping = aes(x=E8q24))+
# Trend lines
geom_line(mapping=aes(y=TIL_mean,group=factor(cancer_type),color=factor(cancer_type)),stat = "identity",size=1.5)+
#geom_errorbar(aes(ymin = TIL_mean - TIL_sd, ymax = TIL_mean + TIL_sd), 
                #width = 0.5,cex=1)+
# 利用 annotate语句添加图中注释
#annotate("text", x = 10, y = 3000, label = "全省合计",size=6,fontface="bold")+
# Axis labels
xlab("8q24 status")+
ylab("Percentage of TIL")+
# Y-axis range
scale_y_continuous(breaks = c(5,7.5,10,12.5,15),limits = c(5,16.98))+
theme_classic()+
# Theme
theme(
axis.text.x = element_text(angle = 90, hjust = 1),
legend.title=element_blank(),
legend.position = c(0.7,0.83),
legend.text=element_text(size=10),
plot.caption = element_text(hjust=0.5, size=15),
axis.text=element_text(size=10),
# 设置轴标题文字大小和文字加粗
axis.title=element_text(size=12,face="bold")
   )+
  theme_bw()+
  theme(axis.text=element_text(family="sans",color="black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 1))->p1
 # 把作好的图存入向量p1


p1


##################趋势百分比堆积柱状图##################################################

library(tidyverse)

# 构造虚拟数据 -- 无实际意义
data <- data.frame(LUAD = sample(x = 1:100, size = 12),
                   STAD = sample(x = 1:100, size = 12),
                   SKCM = sample(x = 1:100, size = 12))
# data$x <- rep(c("WT", "LOH", "Loss"), each = 4)
data$x <- rep(c(1:3), each = 4)
data$group <- rep(c(LETTERS[1:4]), 3)

# 长宽数据转换：
data_long <- data %>% 
  pivot_longer(cols = !c(x, group),
               names_to = "group2",
               values_to = "values") #%>% 
#mutate(x = factor(x, levels = c("WT", "LOH", "Loss")))

head(data_long)



TIL_sum[TIL_sum == "WT"] = 1
TIL_sum[TIL_sum == "gain"] = 2
TIL_sum[TIL_sum == "Amp"] = 3
TIL_sum$x<-as.numeric(TIL_sum$x)

c("#FF7B61", "#FFBC70", "#68DEFC", "#5F90FA")

c("#FF7575", "#FFB27F", "#73D0FF", "#6E73FF")

c("#4552a9", "#9cd6ed", "#807C7C",
           "#faea6b", "#f2791e")
           

c("#FF1414", "#FFFFFF", "#FFFFFF")


c("#9E8E4C", "#615656", "#FFFFFF")




"#4552a9", "#9cd6ed", "#615656",
"#faea6b","#9E8E4C", "#f2791e"



# 绘图：
# 先画一张：
ggplot(TIL_sum[which(TIL_sum$group2 == "PRAD.TMPRSS2-ERG Fusion"),]) +
  geom_area(aes(x, values, fill = group),
            position = "fill") +
  geom_vline(xintercept = 2, color = "white", size = 0.2)+
  scale_fill_manual(values = rev(c("#4552a9", "#9cd6ed","#615656",
                                            "#faea6b", "#f2791e")))+
                                              scale_x_continuous(expand = c(0, 0),
                                                                 breaks = c(1:3), 
                                                                 labels = c("WT", "gain", "Amp")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("")+
  ylab("")+
  ggtitle("PRAD.TMPRSS2-ERG Fusion")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust = 0.7),
        legend.position = "none")

ggsave("plot1.pdf", height = 6, width = 3)

fix(TIL_sum)

# 循环绘图：
p_list <- list()

for (i in 1:length(unique(TIL_sum$group2))) {
  if (i == 1) {
    p_list[[i]] <- ggplot(TIL_sum[which(TIL_sum$group2 == unique(TIL_sum$group2)[i]),]) +
      geom_area(aes(x, values, fill = group),
                position = "fill") +
      geom_vline(xintercept = 2, color = "white", size = 0.2)+
      scale_fill_manual(values = rev(c("#4552a9", "#9cd6ed",  "#615656",
                                                "#faea6b","#9E8E4C", "#f2791e")))+
                                                  scale_x_continuous(expand = c(0, 0),
                                                                     breaks = c(1:3), 
                                                                     labels = c("WT", "gain", "Amp")) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("")+
      ylab("")+
      ggtitle(unique(TIL_sum$group2)[i])+
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks.x = element_blank(), 
            axis.text.x = element_text(angle = 30, hjust = 0.7),
            legend.position = "none")
  } else {
    p_list[[i]] <- ggplot(TIL_sum[which(TIL_sum$group2 == unique(TIL_sum$group2)[i]),]) +
      geom_area(aes(x, values, fill = group),
                position = "fill") +
      geom_vline(xintercept = 2, color = "white", size = 0.2)+
      scale_fill_manual(values = rev(c("#4552a9", "#9cd6ed", "#615656",
                                                "#faea6b", "#f2791e")))+
                                                  scale_x_continuous(expand = c(0, 0),
                                                                     breaks = c(1:3), 
                                                                     labels = c("WT", "gain", "Amp")) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("")+
      ylab("")+
      ggtitle(unique(TIL_sum$group2)[i])+
      theme(plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(), 
            axis.text.y = element_blank(), 
            axis.text.x = element_text(angle = 30, hjust = 0.7),
            legend.position = "none")
  }
}

# 拼图：
library(cowplot)
plot_grid(plotlist = p_list)
plot_grid(plotlist = p_list, ncol = 10, rel_widths = c(1.15,1,1))

ggsave("plot2.pdf", height = 6, width = 8)


p_list$BRCA.Her2
p_list[[2]]




