
####===Development of Stem.Sig===####
library(CytoTRACE)
library(Seurat)
library(dplyr)
library(stringr)
library(cancerclass)
library(cutpointr)
library(OptimalCutpoints)
library(reshape2)
library(future)

rm(list=ls())
gc()


##################################################################################
####=================================1. load dir==============================####
##################################################################################


cells <- 'Malignant'
load(paste0('data/dir/dir_',cells,'.Rdata')) ## loading dir for scRNA-Seq cohorts with both available malignant and stromal/immune cells data. 
dir <- dir_Malignant 



##################################################################################
####===============2. run CytoTRACE for each scRNA-Seq cohort=================####
##################################################################################

#Plsease Download 34 scRNA-Seq datasets from http://tisch.comp-genomics.org/
##For demonstration, we provided scRNA-data and meta data of SRP114962 dataset. 

##Results of cytoTARCE analysis for each scRNA-Seq dataset were provided in 'results/cytoGenes/'


for(n in dir$cohort[3]){ #demonstration using dataset BRCA_SRP114962
  scRNA <- Read10X_h5(paste0('data/tisch_h5/',n,'_expression.h5'))
  phenotype <- read.delim(paste0('data/scRNAmeta/',n,'_CellMetainfo_table.tsv'))
  phenotype <- as.data.frame(phenotype)
  cellType <- phenotype[phenotype$Celltype..malignancy. == 'Malignant cells',]
  df <- scRNA[ ,colnames(scRNA) %in% cellType$Cell]
  df <- as.data.frame(df) 
  
  ####two samples from OV_GSE118828 need to be removed due to same expression level for all genes ####
  df <- df[, !colnames(df) == "PN1-P_dummy"]
  df <- df[, !colnames(df) == "HG3-M1_dummy"]
  ##################################################
  
  results <- CytoTRACE(df,ncores = 40)
  
  ###===Calculate Spearman correlation of cytoTRACE genes===####
  
  exprMatrix <- results[["exprMatrix"]]
  
  cytoTrace <- data.frame(cytoTrace = results$CytoTRACE)
  cytoTrace$cell <- rownames(cytoTrace)
  remove(results)
  
  cytoGenes <- data.frame(gene=rownames(exprMatrix),coef=NA,p=NA)
  
  for(i in 1:nrow(exprMatrix)){
    cor <- cor.test(exprMatrix[i,],cytoTrace$cytoTrace,method = 'pearson')
    cytoGenes$coef[i] <- cor$estimate
    cytoGenes$p[i] <- cor$p.value
  }
  
  cytoGenes$p.adjust <- p.adjust(cytoGenes$p,method = 'BH')
  
  # save(cytoGenes,file=paste0('results/cytoGenes/' , n , '_cytoGenes.Rdata')) 
}




##################################################################################
####==========3. Differential expression genes of malignant cells=============####
##################################################################################

## results were provided in 'results/DE/'

for(n in dir$cohort[3]){ #demonstration using dataset BRCA_SRP114962
  
  scRNA <- Seurat::Read10X_h5(paste0('data/tisch_h5/',n,'_expression.h5'))
  
  phenotype <- read.delim(paste0('data/scRNAmeta/',n,'_CellMetainfo_table.tsv'))
  phenotype <- as.data.frame(phenotype)
  phenotype <- phenotype[,c(1,4:6)]
  
  table(phenotype$Celltype..malignancy.)
  
  phenotype[!phenotype$Celltype..malignancy. == paste0(cells,' cells'),2] <- "control"
  phenotype[phenotype$Celltype..malignancy. == paste0(cells,' cells'),2] <- cells

  meta <- phenotype
  
  sce <- CreateSeuratObject(counts = scRNA, 
                            project = "test")
  
  sce@meta.data$group <- meta$Celltype..malignancy.
  
  #########future multicore##########
  future::plan("multisession", workers = 40)
  ##########################################
  
  DE <- FindMarkers(sce,ident.1 = cells, group.by = 'group',logfc.threshold = 0.25,min.pct = 0.1,base=exp(1)) 
  DE <- DE[DE$p_val_adj<1e-05,]
  
  # save(DE,file = paste0('results/DE/' , n , '_DE.Rdata')) 
}


##################################################################################
####===============================4. get Stem.Sig============================####
##################################################################################


getSig <- function(dir){

ls_Gn <- list()
allGenes <- c()  
  
for(n in dir$cohort){ 
  load(file=paste0('results/cytoGenes/' , n , '_cytoGenes.Rdata'))
  load(paste0('results/DE/' , n , '_DE.Rdata'))
  
  Gx <- cytoGenes[cytoGenes$coef>0 & cytoGenes$p.adjust < 1e-05,] 
  
  Gy <- rownames(DE[DE$avg_logFC>0,]) ##
  Gy <- rownames(DE) ##
  Gy <- Gy[!grepl("^RP[SL]", Gy,ignore.case = F)] ##(ribosome protein free)
  
  Gn <- Gx[Gx$gene %in% Gy,]
  
  allGenes <- c(allGenes,Gn$gene)
  ls_Gn[[n]] <- Gn
}


  allGenes <- unique(allGenes)
  allGenesDf <- data.frame(gene = allGenes)

  for(i in dir$cohort){
    allGenesDf <- left_join(allGenesDf,ls_Gn[[i]][,c('gene','coef')],by='gene')
  }
  
  allGenesDf <- allGenesDf[!is.na(allGenesDf$gene),] 
  rownames(allGenesDf) <- allGenesDf$gene
  allGenesDf <- allGenesDf[,-1]
  colnames(allGenesDf) <- dir$cohort
  genelist <- allGenesDf
  genelist$all_gmean <- compositions::geometricmeanRow(genelist[,1:length(dir$cohort)]) ##spearmanR geometric mean
  sig <- genelist[genelist$all_gmean>0.4,] ##filter genes with spearmanR geometric mean > 0.4
  sig <- sig[order(sig$all_gmean,decreasing = T),]
  sig <- rownames(sig)
  return(sig)
  
}


Stem.Sig <- getSig(dir)

# save(Stem.Sig,file='data/sig/Stem.Sig.Rdata')


##################################################################################
####=============================5. get Stem.SKCM.Sig=========================####
##################################################################################
dir_SKCM <- dir[grepl('SKCM',dir$cohort),] ## SKCM cohort
Stem.SKCM.Sig <- getSig(dir_SKCM) 

# save(Stem.SKCM.Sig,file='data/sig/Stem.SKCM.Sig.Rdata')



##################################################################################
####=============================6. get GCS.Sig===============================####
##################################################################################


getGCS <- function(dir){
  ls_stem <- list()
  allGenes_stem <- c()  
  for(n in dir$cohort){ 
    load(file=paste0('results/cytoGenes/' , n , '_cytoGenes.Rdata'))
    load(paste0('results/DE/' , n , '_DE.Rdata'))
    GCS <- cytoGenes[order(cytoGenes$coef,decreasing = T),'gene'][1:200]
    allGenes_stem <- c(allGenes_stem,GCS)
  }
  return(unique(allGenes_stem))
  
}

GCS.Sig <- getGCS(dir)






