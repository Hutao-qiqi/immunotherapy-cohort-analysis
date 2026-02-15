###===Figure S3===###

rm (list=ls())
gc()

####===loading data===####

library(GSVA)
library(corrplot)
library(data.table)
library(stringr)
library(ggplot2)
source('Code/corplot.R')
source('Code/bnplot.R')


load('data/gmt/stem_gmt.Rdata')
load('data/TCGA/bulkExpMatrix.Rdata')
load('results/GSVA/gsvaSig.Rdata')
load('results/GSVA/gsvaStem.Rdata')

mRNAsi <- read.csv("data/TCGA/mRNAsi.csv")
mRNAsi <- as.data.frame(mRNAsi)

# gsvaStem <- gsva(bulkExpMatrix,stem,method='gsva',parallel.sz=40)
# save(gsvaStem,file='results/GSVA/gsvaStem.Rdata')

stem <- rbind(gsvaStem,gsvaSig)
stem <- as.data.frame(stem[,colnames(stem) %in% mRNAsi$id])
rownames(stem)[12] <- 'Stem.Sig'

mRNAsi <- mRNAsi[mRNAsi$id %in% colnames(stem),]
mRNAsi <- mRNAsi[!duplicated(mRNAsi$id),]
rownames(mRNAsi) <- mRNAsi$id
mRNAsi <- mRNAsi[colnames(stem),]
mRNAsi <- as.data.frame(t(mRNAsi))
mRNAsi <- mRNAsi[6,]
Stem.Sig <- stem[12,]
Stem.other <- rbind(stem[1:11,],mRNAsi)



###===Figure S3A====####
stem.all <- rbind(Stem.other,Stem.Sig)
rn <- rownames(stem.all)
stem.all <- apply(stem.all,2,as.numeric)
stem.all <- as.data.frame(stem.all)
rownames(stem.all) <- rn
coef <- cor(t(stem.all))

col = colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7","#F4A582", "#D6604D" ,"#B2182B", "#67001F"))(200)


corrplot(corr = coef,order="AOE",type="upper",tl.pos="tp",col=col,tl.col='black')
corrplot(corr = coef,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n",col=col)




###===Figure S3B====####
meta <- fread('data/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp')
Im <- fread('data/ImmuneGenes/Im.csv')
coef <- corplot(Stem.Sig,Stem.other,1,'coef',Im,meta)
p <- corplot(Stem.Sig,Stem.other,1,'pValue',Im,meta)

bnplot(coef,p,2)







