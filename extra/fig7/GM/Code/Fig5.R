####===Figure 5===####
rm (list=ls())
gc()

##===loading data===###
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggsci)



load('results/Machine Learning/AUC.Rdata')




###====Fig 5A====####

cohort <- rownames(AUC[['pancancer']]$Stem.Sig)
rename <- c('Validation Cohort',
            'Testing Cohort',
            'Hugo 2016 SKCM',
            'Van 2015 SKCM',
            'Kim 2018 GC',
            'Zhao 2019 GBM',
            'Synder 2017 UC'
)

# AUC <- lapply(AUC,as.data.frame)

auc <- AUC[['pancancer']]

ls_AUC <- lapply(cohort,function(c){
  cc <-lapply(auc,function(auc){auc[c,'AUC']}) %>% do.call(rbind,.)
  names(cc) <- names(auc) 
  cc <- as.data.frame(cc)
  colnames(cc) <- 'AUC'
  return(cc)
}) %>% `names<-`(rename)


circos_cormpare <- function(ls_AUC){
  bardata <- do.call(cbind, ls_AUC)
  bardata <- as.data.frame(bardata)
  colnames(bardata) <- names(ls_AUC)
  bardata <- as.data.frame(t(bardata))
  
  bardata$cohort <- rownames(bardata)
  
  bardata_m <- melt(bardata,value.name = 'y',variable.name = 'x')
  colnames(bardata_m)[1] <- 'sectors'
  bardata_m$y[is.na(bardata_m$y)] = 0
  
  df <- bardata_m
  
  
  Set2 <- brewer.pal(8,"Set2")
  Set1 <- brewer.pal(9,"Set1")
  col <- c(Set2,Set1)
  col[1:2] <- c("#CA6855",'#546B7E')
  # col <- col[1:7]
  
  sectors <- bardata_m$cohort
  x <- bardata_m$x
  y <- bardata_m$y
  x1 <- rep(1:(ncol(bardata)-1),nrow(bardata))
  
  x1 <- x1[order(x1)]
  df$x1 <- x1
  df$sectors <- factor(df$sectors,levels=unique(df$sectors))
  
  circos.clear()
  
  circos.par(gap.degree = c(rep(2,nrow(bardata)-1),15), cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.01, 0.01),start.degree=90)
  circos.initialize(df$sectors, xlim=c(0,ncol(bardata)))
  
  bgcol = c('#1A1A1A',rep('#999999',5))
  
  circos.track(df$sectors, y = df$y,track.height =0.03, bg.col=bgcol,bg.border='white',
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] + mm_y(3), 
                             # CELL_META$ycenter,
                             CELL_META$sector.index,
                             facing = 'inside',
                             niceFacing = T,
                             cex = 1,
                             font = 2)
               })
  
  
  
  circos.track(df$sectors,ylim=c(0,0.85),track.height=0.78,bg.border='white')
  for(i in unique(df$sectors)){
    circos.barplot(df[df$sectors==i,'y'],pos = df[df$sectors==i,'x1'],sector.index = i,col=col[1:7])
  }
  
  circos.yaxis(side='left', at=seq(0,1,0.2), labels=T, sector.index = unique(df$sectors)[1], track.index=2)
  
  # bgcol = c(rep(col[1],3),col[2],col[2],rep(col[1],4))
  
  # circos.track(df$sectors, ylim=c(0,1), bg.col = bgcol,track.height=0.05)
  
  draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
              rou1 = 0.15, col = "black", border = "#EEEEEE")
  
  
  text(0,0,"AUC",col = 'white',cex =1.5,font = 2 )
  
  
  sig <- col[1:(ncol(bardata)-1)]
  
  names(sig) <- unique(df$x)
  
  lgd_cancer = Legend(title = "Signature", at = names(sig),
                      legend_gp = gpar(fill = sig)) 
  
  draw(lgd_cancer, x = unit(1.3, "snpc"), just = "left")
  
}

circos_cormpare(ls_AUC[-1])




###====Fig 5B====####

#Heatmap 

hdata <- ls_AUC[-c(1)] %>% do.call(cbind,.) %>% `colnames<-`(names(ls_AUC)[-1])

hdata <- hdata[order(hdata[,1],decreasing = T) ,]


Heatmap(as.matrix(hdata), name = "AUC",
        # column_title = NULL,
        # row_title = NULL,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(c( 0.5, 0.80), c("#FFFFE3",  "#004628")),
        # show_row_names = T, show_column_names = F,
        rect_gp = gpar(col = "black", lwd = 2),
        width = nrow(hdata)*unit(6, "mm"),
        height = nrow(hdata)*unit(6, "mm"),
        # na_col = 'white',
        # column_names_side = c('top'),
        # row_split = c(rep('a',8),rep('b',8)),
        column_split = c('Combined',rep('Individual',5)),
        # top_annotation = column_ha
)




###====Fig 5C====####

SKCM <- AUC[['SKCM']] %>% do.call(rbind,.) %>% `colnames<-`('AUC') %>% as.data.frame
SKCM$signatures <- rownames(SKCM)
SKCM <- SKCM[order(SKCM$AUC,decreasing = T),]


Fig.5C <- ggplot(data=SKCM,aes(reorder(signatures,-AUC),fill=AUC)) + 
  geom_col(aes(y = AUC)) + 
  geom_text(data=SKCM, aes(y=AUC,label=round(AUC,2)),nudge_y = 0.02 ) + 
  theme_classic()+
  scale_fill_gradient(low = "#FFFFE3",high = "#004628")


Fig.5C









