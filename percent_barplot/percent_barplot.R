rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

# Load packages
library(ggplot2)



# Input table
data <- read.table("13_potential_therapeutic_targets.txt",header = T,sep = '\t',check.names = F)



# Example data (not used below)
data1 <- data.frame(
  Gene = paste0("Cancer", 1:21),
  Log2FC = c(runif(15, -0.5,0), runif(6, 0, 0.5)),
  "Log10qvalue" = c(runif(6, 1, 3),runif(15, 0, 1)))

# Preview
head(data)

# Color palettes (reference)
# browns: c("#f6fafd", "#FFCC66", "#FF9933", "#CC6600")
# purples: c("#f6fafd", "#CC99FF", "#9933FF", "#9933CC")
# blues:   c("#f6fafd", "#c8dfef", "#6fa6d1", "#2c49a2")


ggplot(data)+
  # Reference line
  geom_hline(yintercept = c( 0.2), 
             color = "black", linetype = "dotted")+
  # Bar plot
  geom_col(aes(x = reorder(cancer_type, rho), y = rho, fill = log10_fdr),
           color = "white")+
  # Zero line
  geom_hline(yintercept = 0, color = "black")+
  # Gradient fill
  scale_fill_gradientn(name="-Log10_\nq-value", # 修改图例标题
                       colours = c( "#FFD4D4","#EB7B7B","#D14141","#A80000" ),
                       breaks = 0:8,
                       labels = paste0(0:8, ".0"))+
  
  
#棕色c("#f6fafd","#FFCC66", "#FF9933","#CC6600")
  #紫色c("#f6fafd", "#CC99FF", "#9933FF", "#9933CC")
  #c("#f6fafd", "#c8dfef", "#6fa6d1", "#2c49a2")
  
  
####

  # Axis labels
  xlab("")+
  ylab("Spearman’s ρ")+
  # Reduce plot padding
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(limits = c(-0.1,0.4),breaks=seq(-0.1,1,0.4))+
  #scale_y_discrete(expand = c(0.008, 0.008))+
  # Annotation
  annotate(geom = "text", y = 0.39, x = 14.8, label = "")+
  # Title
  ggtitle("")+
  # Flip coordinates
  coord_flip()+
  # Theme
  theme_light()+
  theme(axis.text=element_text(family="sans",color="black"))+
  theme(panel.grid = element_blank(),  # 去掉网格线
        plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中、字体
        axis.ticks.y = element_blank(), # 去掉y轴刻度线
        axis.title.x = element_text(size = 10),  # x轴标题大小
        axis.text.y = element_text(size = 10),  # y轴刻度大小
        panel.border = element_rect(color = "black"),
        legend.position=c(0.27,0.25), legend.justification=c(1,0))  # 图例位置

# Save
ggsave("barplot.pdf", height = 5, width = 4)  




######################GSEA###################################################





