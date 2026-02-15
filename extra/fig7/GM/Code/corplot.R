#corplot sig matrix, bulk Matrix(mode 2) / immune Matrix(mode 1), Im
corplot <- function(sigMatrix,immune,mode,return,Im,meta){
  library(dplyr)
  id <- intersect(colnames(immune),colnames(sigMatrix))
  cor_sig_all <- as.data.frame(sigMatrix)
  cor_sig_all <- cor_sig_all[,id]
  if(mode == 1){
    library(corrplot)
    cor_Im_all <- as.data.frame(immune[,id])
    
    cancer <- data.frame(sample=colnames(cor_Im_all))
    cancer <- left_join(cancer,meta[,c(1,3)],by='sample')
    CA <- unique(cancer$`cancer type abbreviation`)
    CA <- CA[!is.na(CA)]
    #CA <- c(CA,'PanCancer')
    
    Im_sub <- rownames(cor_Im_all)
    sig_sub <- rownames(cor_sig_all)
    
    list_cor <- list()
    list_p <- list()
    
    for (c in CA){
      if(c == "PanCancer"){ 
        sample <- cancer[,'sample'] } 
      else {
        sample <- cancer[cancer$`cancer type abbreviation`== c & !is.na(cancer$`cancer type abbreviation`),'sample']
      }  
      
      cor_sig <- cor_sig_all[, colnames(cor_sig_all) %in% sample]
      cor_Im <- cor_Im_all[, colnames(cor_Im_all) %in% sample]
      cor_df <- data.frame()
      cor_p <- data.frame()
      for ( k in sig_sub){
        for(j in Im_sub){
          cor <- cor.test(as.numeric(cor_Im[j,]),as.numeric(cor_sig[k,]),method = "spearman",exact = FALSE)
          cor_df[j,k] <- cor$estimate
          cor_p[j,k] <- cor$p.value
        }
      }
      list_cor[[c]] <- cor_df
      list_p[[c]] <- cor_p
    }
    
    
    pheatmapData <- list_cor[[1]]
    for( i in 2:length(list_cor)){
      pheatmapData <- cbind(pheatmapData,list_cor[[i]])
    }
    
    pData <- list_p[[1]]
    for( i in 2:length(list_p)){
      pData <- cbind(pData,list_p[[i]])
    }
    
    
    colnames(pheatmapData) <- names(list_cor)
    colnames(pData) <- names(list_p)
    
    
    for(j in 1:ncol(pheatmapData)){
      pheatmapData[is.na(pheatmapData[,j]),j] <- 0
    }
    
    for(j in 1:ncol(pData)){
      pData[is.na(pData[,j]),j] <- 1
    }
    
    
    
    # p <- corrplot(as.matrix(pheatmapData),
    #               p.mat = as.matrix(pData),
    #               insig = 'blank',
    #               tl.col = 'black', 
    #               col = colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
    #                                        "#F4A582", "#D6604D" ,"#B2182B", "#67001F"))(200))
    
    if(return == "coef") {
      return(pheatmapData)
    } else if (return == "pValue") {
      return(pData)
    } else if ( return == 'p'){
      return(p)
    }
    
  } else if (mode == 2) {
    library(pheatmap)
    cor_Im_all <- as.data.frame(immune[Im$Genes,id])
    
    cancer <- data.frame(sample=colnames(cor_Im_all))
    cancer <- left_join(cancer,meta[,c(1,3)],by='sample')
    CA <- unique(cancer$`cancer type abbreviation`)
    CA <- CA[!is.na(CA)]
    # CA <- c(CA,'PanCancer')
    
    Im_sub <- rownames(cor_Im_all)
    sig_sub <- rownames(cor_sig_all)
    
    list_cor <- list()
    list_p <- list()
    
    for (c in CA){
      if(c == "PanCancer"){ 
        sample <- cancer[,'sample'] } 
      else {
        sample <- cancer[cancer$`cancer type abbreviation`== c & !is.na(cancer$`cancer type abbreviation`),'sample']
      }  
      
      cor_sig <- cor_sig_all[, colnames(cor_sig_all) %in% sample]
      cor_Im <- cor_Im_all[, colnames(cor_Im_all) %in% sample]
      cor_df <- data.frame()
      cor_p <- data.frame()
      for ( k in sig_sub){
        for(j in Im_sub){
          cor <- cor.test(as.numeric(cor_Im[j,]),as.numeric(cor_sig[k,]),method = "spearman",exact = FALSE)
          cor_df[j,k] <- cor$estimate
          cor_p[j,k] <- cor$p.value
        }
      }
      list_cor[[c]] <- cor_df
      list_p[[c]] <- cor_p
    }
    
    
    pheatmapData <- list_cor[[1]]
    for( i in 2:length(list_cor)){
      pheatmapData <- cbind(pheatmapData,list_cor[[i]])
    }
    
    pData <- list_p[[1]]
    for( i in 2:length(list_p)){
      pData <- cbind(pData,list_p[[i]])
    }
    
    
    colnames(pheatmapData) <- names(list_cor)
    colnames(pData) <- names(list_p)
    
    
    for(j in 1:ncol(pheatmapData)){
      pheatmapData[is.na(pheatmapData[,j]),j] <- 0
    }
    
    for(j in 1:ncol(pData)){
      pData[is.na(pData[,j]),j] <- 1
    }
    
    
    or <- Im[Im$Genes %in% rownames(pheatmapData),]
    or <- or[order(or$Pathway),]
    or <- or$Genes
    pheatmapData <- pheatmapData[or,]
    annotation_row <- data.frame(Genes=rownames(pheatmapData))
    annotation_row <- left_join(annotation_row,Im,by='Genes')
    rownames(annotation_row) <- annotation_row$Genes
    annotation_row <- annotation_row[,-1]
    
    bk <- c(seq(-0.8,-0.01,by=0.001),seq(0,0.8,by=0.001))
    
    # p <- pheatmap(pheatmapData,cluster_cols = F,cluster_rows = F,show_rownames = T,show_colnames = T, 
    #               color = c(colorRampPalette(colors = c("darkblue","white"))(length(bk)/2), 
    #                         colorRampPalette(colors = c('white',"darkred"))(length(bk)/2)),
    #               breaks =  bk,
    #               cellwidth = 10,
    #               cellheight = 10,
    #               annotation_row = annotation_row,
    #               # gaps_col = ncol(pheatmapData)-1,
    # )
    # 
    if(return == "coef") {
      return(pheatmapData)
    } else if (return == "pValue") {
      return(pData)
    } else if (return == 'p'){
      return(p)
    }
    
  }
  
}