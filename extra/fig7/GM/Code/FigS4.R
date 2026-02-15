###===Figure S4===###

rm (list=ls())
gc()
####===loading data===####
library(data.table)
library(stringr)
library(ggplot2)
source('Code/corplot.R')
source('Code/bnplot.R')
load('data/TCGA/bulkExpMatrix.Rdata')
load('results/GSVA/gsvaSig.Rdata')

Genes <- c('SKAP1','RBP5','PTGDS','LAT','EIF1AY','CETP','CD79B','CD1D','CCR6')
Genes <- rev(Genes)

TLS_df <- bulkExpMatrix[ Genes, colnames(gsvaSig)]

TLS_df <- rbind(TLS_df,apply(TLS_df,2,mean))

rownames(TLS_df)[10] <- 'TLS scores'

meta <- fread('data/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp')
Im <- fread('data/ImmuneGenes/Im.csv')

coef <- corplot(gsvaSig,TLS_df,1,'coef',TLS,meta)
p <- corplot(gsvaSig,TLS_df,1,'pValue',TLS,meta)

bnplot(coef,p,mode=2,limit=0.8)
