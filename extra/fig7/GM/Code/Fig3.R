####===Figure 3===####
rm (list=ls())
gc()




###=== GSVA ===###
library(GSVA)
library(GSEABase)
library(data.table)
source('Code/corplot.R')
source('Code/bnplot.R')
load('data/TCGA/bulkExpMatrix.Rdata')

# gmt <- getGmt('data/gmt/Stem.Sig.gmt')
# gsvaSig <- gsva(bulkExpMatrix,gmt,method='gsva',parallel.sz=30)
# gmt <- getGmt('data/gmt/h.all.v7.4.symbols.gmt')
# gsvaH <- gsva(bulkExpMatrix,gmt,method='gsva',parallel.sz=30)

# save(gsvaSig,file='results/GSVA/gsvaSig.Rdata')
# save(gsvaH,file='results/GSVA/gsvaH.Rdata')
load('results/GSVA/gsvaSig.Rdata')
load('results/GSVA/gsvaH.Rdata')

####===Fig. 3A===####
library(circlize)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(dplyr)

meta <- fread('data/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp')
Im <- fread('data/ImmuneGenes/Im.csv')
coef <- corplot(gsvaSig, bulkExpMatrix, mode = 2, return = 'coef',Im = Im, meta = meta)
IM75vsGSVAsig <- coef

annotation_row <- data.frame(Genes=rownames(IM75vsGSVAsig))
annotation_row <- left_join(annotation_row,Im,by='Genes')
rownames(annotation_row) <- annotation_row$Genes
annotation_row <- annotation_row[,-1]
annotation_row['PDCD1LG2','Function'] = 'Inhibitory'
annotation_row[annotation_row$Function == 'N/A','Function'] <- 'HLA'


annotation_row$order <- annotation_row$Pathway
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Atigen Presentation", '1')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Cell adhesion", '2')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Co-stimulator", '3')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Co-inhibitor", '4')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Other", '7')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"Receptor", '5')
annotation_row$order <- str_replace_all(annotation_row[,'order'],"ligand", '6')

annotation_row <- annotation_row[order(annotation_row$order),]

mat <- as.matrix(IM75vsGSVAsig[rownames(annotation_row),])

range(mat)

circlize_plot <- function(){
  
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                        just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  
  #####
  circos.clear()
  
  circos.par("gap.after" = c(1,1,1,1,1,1,10),
             'start.degree' = 90
  )
  
  path <- data.frame(Path = annotation_row$Pathway)
  rownames(path) <- rownames(annotation_row)
  
  col_path <- brewer.pal(7,"Set3")
  names(col_path) <- unique(annotation_row$Pathway)
  
  circos.heatmap(path, split = path$Path, col = col_path, track.height = 0.02,rownames.side = 'outside')
  
  col_mat = colorRamp2(c(-1, 0, 1),  c("#4783B4", "white", "#8D4A2A")) ##heatmap color###
  
  # col_mat = colorRamp2(c(-1, -0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),  c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B", "#67001F"))
  
  
  circos.heatmap(mat, split = annotation_row$Pathway, col = col_mat, track.height = 0.7)
  
  col_fun <- c("#BEBEBE","#E59E00","#323232")
  names(col_fun) <-  unique(annotation_row$Function)
  circos.heatmap(annotation_row$Function, col = col_fun, track.height = 0.01)
  
  circos.yaxis("left", at=1:30,labels = NULL, track.index =  3, sector.index = "Atigen Presentation",labels.cex = 0.7,tick=T,labels.niceFacing=F)
  
  lgd_path = Legend(title = "Track1-Pathway", at = names(col_path), 
                    legend_gp = gpar(fill = col_path),nrow=1)
  lgd_cor = Legend(title = "Track2-Corelation Spearman", col_fun = col_mat,direction = 'horizontal')
  
  lgd_fun = Legend(title = "Track3-Function", at = names(col_fun), 
                   legend_gp = gpar(fill = col_fun),nrow=1)
  
  upViewport()
  
  h = dev.size()[2]
  lgd_list = packLegend(lgd_path, lgd_cor, lgd_fun, max_height = unit(0.9*h, "inch"))
  draw(lgd_list, x = circle_size, just = "left")
  
  
}


circlize_plot()





####===Fig. 3B===####
library(MCPcounter)
# MCPcounter <- MCPcounter.estimate(bulkExpMatrix,featuresType = 'HUGO_symbols')
# save(MCPcounter,file='results/MCPcounter/MCPcounter.Rdata')
load('results/MCPcounter/MCPcounter.Rdata')
MCPcounter <- MCPcounter[1:7,] ## tumor infiltrated lymphocytes

coef <- corplot(gsvaSig, MCPcounter, mode = 1, return = 'coef',Im = Im, meta = meta)
p <- corplot(gsvaSig, MCPcounter, mode = 1, return = 'pValue',Im = Im, meta = meta)


Fig.3B <-  bnplot(t(coef), t(p), mode = 2, size = c(2,4), col = c("#4783B4", "white", "#8D4A2A"), limit = 1)

Fig.3B

####===Fig. 3C===####

coef <- corplot(gsvaSig,gsvaH,1,'coef',Im,meta)
p <- corplot(gsvaSig,gsvaH,1,"pValue",Im,meta)


mean <- apply(coef,1,mean)

mean <- mean[order(mean,decreasing = T)]

abs <- abs(mean)

abs <- abs[order(abs,decreasing = T)]

order <- mean[ names(mean) %in% names(abs)[1:10] ]

order <- names(order)

order <- rev(order)

Fig.3C <- bnplot(coef[order,],p[order,],mode=2,size=c(1,6),col=c('#3B7FC7',"yellow",'#D41E2A'),limit=1)
Fig.3C

####===Fig. 3D and 3E===####

library(ggpubr)
library(ggrepel)


load(file=paste0('data/TCGA/medianITH.Rdata'))
Fig.3D <- ggplot(medianITH,aes(x=mITH,y=mSig))+
  geom_point(data = medianITH,aes(x=mITH,y=mSig),color="#3B7FC7",size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_text_repel(aes(mITH, mSig, label =  CancerType))+
  #facet_wrap(~CancerType,scales = 'free')+
  # geom_smooth(method='lm',se=F,show.legend=F,linetype=1,color='darkred',size = 0.6)+
  geom_smooth(method='lm',se=T,show.legend=F,linetype=1,size = 0.6,color='#333333')+
  stat_cor(method = "spearman",show.legend = F)+
  xlab('median ITH')+
  ylab("median GSVA of Stem.Sig")

Fig.3D

load(file=paste0('data/TCGA/medianTMB.Rdata'))
Fig.3E <- ggplot(medianTMB,aes(x=log10mTMB,y=mSig))+
  geom_point(data = medianTMB,aes(x=log10mTMB,y=mSig),color="#3B7FC7",size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_text_repel(aes(log10mTMB, mSig, label =  CancerType))+
  #facet_wrap(~CancerType,scales = 'free')+
  # geom_smooth(method='lm',se=F,show.legend=F,linetype=1,color='darkred',size = 0.6)+
  geom_smooth(method='lm',se=T,show.legend=F,linetype=1,size = 0.6,color='#333333')+
  stat_cor(method = "spearman",show.legend = F)+
  xlab('median log10TMB')+
  ylab("median GSVA of Stem.Sig")
Fig.3E








