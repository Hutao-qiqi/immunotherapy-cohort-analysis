bnplot <- function(coef,p,mode,col=NA,limit=1,size=c(1,6)){
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  coef <- as.data.frame(coef)
  rownames(coef) <- str_replace_all(rownames(coef),'HALLMARK','')
  rownames(coef) <- str_replace_all(rownames(coef),'_',' ')
  mdf <- coef
  
  
  order <- apply(mdf,1,mean)
  order <- order[order(order)]
  names(order)
  
  
  
  
  mdf$path <- rownames(mdf)
  mdf <-reshape2::melt(mdf,value.name = 'Spearman Correlation R',variable.name = 'cancer')
  
  mp <- as.data.frame(p)
  mp$path <- rownames(mp)
  mp <- reshape2::melt(mp,value.name = 'p value',variable.name = 'cancer')
  
  mdf$`p value` <- mp$`p value`
  
  mdf$FDR <- p.adjust(mdf$`p value`,method = 'BH')
  
  mdf$`-log10(q value)` <- -log10(mdf$FDR)
  
  if (mode == 1){
    mdf$path <- factor(mdf$path,levels = names(order))
  } else if (mode == 2) {
    mdf$path <- factor(mdf$path,levels = rownames(coef))
  }
  
  # for(i in 1:nrow(mdf)){
  #   # if(mdf$FDR[i] >= 0.05 | abs(mdf$`Spearman Correlation R`)[i] <0.4) {mdf$`-log10(q value)`[i] = NA}
  #   if(mdf$FDR[i] >= 0.05 ) {mdf$`-log10(q value)`[i] = NA}
  # }
  
  
  # mdf1 <- mdf[mdf$path %in% names(order)[c(1:5,length(order):(length(order)-19))],]
  
  mdf1 <- mdf
  
  median(mdf1$`Spearman Correlation R`)
  
  # b <- c(-1.0,-0.8,-0.4, 0, 0.4,0.8,1.0)
  
  if(is.na(col)){
    col = colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                             "#F4A582", "#D6604D" ,"#B2182B", "#67001F"))(200)
  }
  
  
  ggplot(mdf1, aes(x = cancer , y = path)) + 
    geom_point(aes(size = `-log10(q value)`, color = `Spearman Correlation R`)) +
    # scale_color_gradient2(low = "#3B7FC7", mid = "#F8EC27", high = "#D41E2A", space = "Lab" ) +
    # scale_colour_gradient2()+
    scale_colour_gradientn(colours =col,limits = c(-limit,limit))+
    # scale_colour_gradientn(colours = c("#4783B4", 'white', "#8D4A2A"), limits = c(-0.7,0.7))+
    scale_size(range = size)+  # Adjust the range of points size
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0))+
    scale_x_discrete(position = "top") +
    labs(x='',y='')
  # labs(x='',y='Top 25 Hallmark pathways')
  
}
