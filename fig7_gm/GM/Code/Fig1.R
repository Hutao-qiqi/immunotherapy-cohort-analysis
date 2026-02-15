###==Figure 1===###

rm (list=ls())

gc()

#====loading data===####


library(CytoTRACE)
library(Seurat)
library(dplyr)
library(ggpubr)
source('Code/plotCytoTRACE.v2.0.R') #modified plotCytoTRACE function implemented in CytoTRACE
source('Code/R_rainclouds.R')




####====Fig. 1A and 1B===####
#GSE115978
cohort <- "SKCM_GSE115978_aPD1" 
load(paste0("data/",cohort,'/scRNA.Rdata'))
load(paste0('data/',cohort,'/phenotype.Rdata'))

cell <- phenotype[phenotype$cellType == 'Malignant cells' ,'Cell'] 
response <- phenotype[phenotype$cellType == 'Malignant cells' ,'Response']
names(response) <- phenotype[phenotype$cellType == 'Malignant cells' ,'Cell']

scRNA <- scRNA[,colnames(scRNA) %in% cell]
scRNA <- as.data.frame(scRNA)
results <- CytoTRACE(scRNA)
  
Fig.1A <- plotCytoTRACE.v2.0(results,response)



boxdata <- data.frame(Cell = names(results$CytoTRACE),cytoTRACE=results$CytoTRACE)
boxdata <- left_join(boxdata,phenotype[,c('Cell','Response')],by='Cell')
boxdata$Response <- factor(boxdata$Response,levels = c("NR","TN"))


Fig.1B <- ggplot(boxdata, aes(x = Response, y = cytoTRACE, fill = Response)) +
  geom_flat_violin(aes(fill = Response),position = position_nudge(x = 0.1, y = 0),  trim = TRUE, alpha = .5, colour = NA)+
  geom_point(aes(x = .55, y = cytoTRACE, colour = Response),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = Response, y = cytoTRACE, fill = Response),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme_classic()+
  geom_signif(comparisons = list(c("NR","TN")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12))+
  # ylim(0,1)+
  ylab('CytoTRACE scores')

Fig.1B

####====Fig. 1C and 1D===####

cohort <- "BCC_GSE123813_aPD1"
load(paste0("data/",cohort,'/scRNA.Rdata')) 
load(paste0("data/",cohort,'/phenotype.Rdata'))

cell <- phenotype[phenotype$cellType == 'Malignant cells' ,'Cell'] 
response <- phenotype[phenotype$cellType == 'Malignant cells' ,'Response']
names(response) <- phenotype[phenotype$cellType == 'Malignant cells' ,'Cell']

scRNA <- scRNA[,colnames(scRNA) %in% cell]
scRNA <- as.data.frame(scRNA)
results <- CytoTRACE(scRNA)


Fig.1C <- plotCytoTRACE.v2.0(results,response)

boxdata <- data.frame(Cell = names(results$CytoTRACE),cytoTRACE=results$CytoTRACE)
boxdata <- left_join(boxdata,phenotype[,c('Cell','Response')],by='Cell')
boxdata$Response <- factor(boxdata$Response,levels = c('NR','R'))


Fig.1D <- ggplot(boxdata, aes(x = Response, y = cytoTRACE, fill = Response)) +
  geom_flat_violin(aes(fill = Response),position = position_nudge(x = 0.1, y = 0),  trim = TRUE, alpha = .5, colour = NA)+
  geom_point(aes(x = .55, y = cytoTRACE, colour = Response),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = Response, y = cytoTRACE, fill = Response),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme_classic()+
  geom_signif(comparisons = list(c("NR","R")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12))+
  ylab('CytoTRACE scores')

Fig.1D





