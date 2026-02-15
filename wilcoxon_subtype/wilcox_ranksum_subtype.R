
rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(tidyverse)
})



HNSC_HPV_Neg <- read.table("HPV-HNSC_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
HNSC_HPV_Pos <- read.table("HPV+HNSC_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
SKCM_RAS_Hotspot <- read.table("RAS_hotspot_mut_SKCM_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
SKCM_BRAF_Hotspot <- read.table("BRAF_hotspot_mut_SKCM_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_EGFR_mut <- read.table("EGFR_mut_LUAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_KRAS_hotspot <- read.table("KRAS_hotspot_mut_LUAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_STK11_mut <- read.table("STK11_mut_LUAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_CIN <- read.table("CIN_STAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_GS <- read.table("GS_STAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_MSI <- read.table("MSI_STAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
ESCA_CIN <- read.table("CIN_ESCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
ESCA_ESCC <- read.table("ESCC_ESCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
READ_CIN <- read.table("CIN_READ_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LGG_IDHmut_codel <- read.table("IDHmut-codel_LGG_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LGG_IDHmut_non_codel <- read.table("IDHmut-non-codel_LGG_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Differentiated <- read.table("differentiated_OV_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Immunoreactive <- read.table("Immunoreactive_OV_annoptation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Mesenchymal <- read.table("Mesenchymal_OV_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Proliferative <- read.table("Proliferative_OV_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Basal <- read.table("Basal_BRCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Her2 <- read.table("Her2_BRCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_LumA <- read.table("LumA_BRCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_LumB <- read.table("LumB_BRCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Normal <- read.table("Normal_BRCA_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_CN_HIGH <- read.table("CN_HIGH_UCEC_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_CN_LOW <- read.table("CN_low_UCEC_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_POLE <- read.table("POLE_UCEC_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)
PRAD_TMPRSS2_ERG_Fusion <- read.table("TMPRSS2-ERG-fusion_PRAD_annotation.txt",header=T,row.names = 1,sep = '\t',check.names = F)




b<-rbind(HNSC_HPV_Neg,HNSC_HPV_Pos,SKCM_RAS_Hotspot,
             SKCM_BRAF_Hotspot,LUAD_EGFR_mut,LUAD_KRAS_hotspot,
             LUAD_STK11_mut,STAD_CIN,STAD_GS,STAD_MSI,ESCA_CIN,
             ESCA_ESCC,READ_CIN,LGG_IDHmut_codel,LGG_IDHmut_non_codel,
             OV_Differentiated,OV_Immunoreactive,OV_Mesenchymal,
             OV_Proliferative,BRCA_Basal,BRCA_Her2,BRCA_LumA,BRCA_LumB,
             BRCA_Normal,UCEC_CN_HIGH,UCEC_CN_LOW,UCEC_POLE,
             PRAD_TMPRSS2_ERG_Fusion)

write.table(b, "18cancer_subtype_annotation.txt",sep="\t", quote=F,row.names = T,col.names = T)





patient_annotation <- read.table("18_cancer_subtype_annotation(immune_score).txt",header=T,row.names=1,sep = '\t',check.names = F)
immune_score <- read.table("BCR_Shnnon+Richness_18_cancer_type.txt",header=T,row.names=1,sep = '\t',check.names = F)


HNSC_HPV_Neg <- patient_annotation[patient_annotation$"Subtype_Selected"== "HNSC_HPV_Neg",]
HNSC_HPV_Pos <- patient_annotation[patient_annotation$"Subtype_Selected"== "HNSC_HPV_Pos",]
SKCM_RAS_Hotspot <- patient_annotation[patient_annotation$"Subtype_Selected"== "SKCM_RAS_Hotspot",]
SKCM_BRAF_Hotspot <- patient_annotation[patient_annotation$"Subtype_Selected"== "SKCM_BRAF_Hotspot",]
LUAD_EGFR_mut <- patient_annotation[patient_annotation$"Subtype_Selected"== "LUAD_EGFR_mut",]
LUAD_KRAS_hotspot <- patient_annotation[patient_annotation$"Subtype_Selected"== "LUAD_KRAS_hotspot",]
LUAD_STK11_mut <- patient_annotation[patient_annotation$"Subtype_Selected"== "LUAD_STK11_mut",]
STAD_CIN <- patient_annotation[patient_annotation$"Subtype_Selected"== "STAD_CIN",]
STAD_GS <- patient_annotation[patient_annotation$"Subtype_Selected"== "STAD_GS",]
STAD_MSI <- patient_annotation[patient_annotation$"Subtype_Selected"== "STAD_MSI",]
ESCA_CIN <- patient_annotation[patient_annotation$"Subtype_Selected"== "ESCA_CIN",]
ESCA_ESCC <- patient_annotation[patient_annotation$"Subtype_Selected"== "ESCA_ESCC",]
READ_CIN <- patient_annotation[patient_annotation$"Subtype_Selected"== "READ_CIN",]
LGG_IDHmut_codel <- patient_annotation[patient_annotation$"Subtype_Selected"== "LGG_IDHmut_codel",]
LGG_IDHmut_non_codel <- patient_annotation[patient_annotation$"Subtype_Selected"== "LGG_IDHmut_non_codel",]
OV_Differentiated <- patient_annotation[patient_annotation$"Subtype_Selected"== "OV_Differentiated",]
OV_Immunoreactive <- patient_annotation[patient_annotation$"Subtype_Selected"== "OV_Immunoreactive",]
OV_Mesenchymal <- patient_annotation[patient_annotation$"Subtype_Selected"== "OV_Mesenchymal",]
OV_Proliferative <- patient_annotation[patient_annotation$"Subtype_Selected"== "OV_Proliferative",]
BRCA_Basal <- patient_annotation[patient_annotation$"Subtype_Selected"== "BRCA_Basal",]
BRCA_Her2 <- patient_annotation[patient_annotation$"Subtype_Selected"== "BRCA_Her2",]
BRCA_LumA <- patient_annotation[patient_annotation$"Subtype_Selected"== "BRCA_LumA",]
BRCA_LumB <- patient_annotation[patient_annotation$"Subtype_Selected"== "BRCA_LumB",]
BRCA_Normal <- patient_annotation[patient_annotation$"Subtype_Selected"== "BRCA_Normal",]
UCEC_CN_HIGH <- patient_annotation[patient_annotation$"Subtype_Selected"== "UCEC_CN_HIGH",]
UCEC_CN_LOW <- patient_annotation[patient_annotation$"Subtype_Selected"== "UCEC_CN_LOW",]
UCEC_POLE <- patient_annotation[patient_annotation$"Subtype_Selected"== "UCEC_POLE",]
PRAD_TMPRSS2_ERG_Fusion <- patient_annotation[patient_annotation$"Subtype_Selected"== "PRAD_TMPRSS2_ERG_Fusion",]



### Full version
cancer_18types_annotation <- list(HNSC_HPV_Neg,HNSC_HPV_Pos,SKCM_RAS_Hotspot,
                                       SKCM_BRAF_Hotspot,LUAD_EGFR_mut,LUAD_KRAS_hotspot,
                                       LUAD_STK11_mut,STAD_CIN,STAD_GS,STAD_MSI,ESCA_CIN,
                                       ESCA_ESCC,READ_CIN,LGG_IDHmut_codel,LGG_IDHmut_non_codel,
                                       OV_Differentiated,OV_Immunoreactive,OV_Mesenchymal,
                                       OV_Proliferative,BRCA_Basal,BRCA_Her2,BRCA_LumA,BRCA_LumB,
                                       BRCA_Normal,UCEC_CN_HIGH,
                                       PRAD_TMPRSS2_ERG_Fusion)

### Reduced version
cancer_18types_annotation <- list(HNSC_HPV_Neg,SKCM_RAS_Hotspot,
                                  SKCM_BRAF_Hotspot,LUAD_EGFR_mut,LUAD_KRAS_hotspot,
                                  LUAD_STK11_mut,STAD_CIN,STAD_GS,STAD_MSI,ESCA_CIN,
                                  ESCA_ESCC,READ_CIN,LGG_IDHmut_non_codel,
                                  OV_Differentiated,OV_Immunoreactive,OV_Mesenchymal,
                                  OV_Proliferative,BRCA_Basal,BRCA_Her2,BRCA_LumA,BRCA_LumB,
                                  BRCA_Normal,UCEC_CN_HIGH,
                                  PRAD_TMPRSS2_ERG_Fusion)


### Full version
names(cancer_18types_annotation) <- c("HNSC_HPV_Neg","HNSC_HPV_Pos","SKCM_RAS_Hotspot",
                                      "SKCM_BRAF_Hotspot","LUAD_EGFR_mut","LUAD_KRAS_hotspot",
                                      "LUAD_STK11_mut","STAD_CIN","STAD_GS","STAD_MSI","ESCA_CIN",
                                      "ESCA_ESCC","READ_CIN","LGG_IDHmut_codel","LGG_IDHmut_non_codel",
                                      "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
                                      "OV_Proliferative","BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB",
                                      "BRCA_Normal","UCEC_CN_HIGH",
                                      "PRAD_TMPRSS2_ERG_Fusion")

### Reduced version
names(cancer_18types_annotation) <- c("HNSC_HPV_Neg","SKCM_RAS_Hotspot",
                                      "SKCM_BRAF_Hotspot","LUAD_EGFR_mut","LUAD_KRAS_hotspot",
                                      "LUAD_STK11_mut","STAD_CIN","STAD_GS","STAD_MSI","ESCA_CIN",
                                      "ESCA_ESCC","READ_CIN","LGG_IDHmut_non_codel",
                                      "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
                                      "OV_Proliferative","BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB",
                                      "BRCA_Normal","UCEC_CN_HIGH",
                                      "PRAD_TMPRSS2_ERG_Fusion")





cancer_18types_immune_score <- list()


for (i in seq_along(cancer_18types_annotation)) {
  # Current annotation table
  annotation <- cancer_18types_annotation[[i]]
  immune_score1 <- immune_score[rownames(annotation),]
  immune_score1 <-as.data.frame(immune_score1)
  immune_score1$"8q24"<-annotation$"8q24"
  immune_score1 <-na.omit(immune_score1)
  
  cancer_18types_immune_score[[i]] <- immune_score1
}

names(cancer_18types_immune_score) <- c("HNSC_HPV_Neg","SKCM_RAS_Hotspot",
                                        "SKCM_BRAF_Hotspot","LUAD_EGFR_mut","LUAD_KRAS_hotspot",
                                        "LUAD_STK11_mut","STAD_CIN","STAD_GS","STAD_MSI","ESCA_CIN",
                                        "ESCA_ESCC","READ_CIN","LGG_IDHmut_non_codel",
                                        "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
                                        "OV_Proliferative","BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB",
                                        "BRCA_Normal","UCEC_CN_HIGH",
                                        "PRAD_TMPRSS2_ERG_Fusion")






HNSC_HPV_Neg <- read.table("HNSC_HPV-_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
HNSC_HPV_Pos <- read.table("HNSC_HPV+_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
SKCM_RAS_Hotspot <- read.table("SKCM_RAS_Hotspot_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
SKCM_BRAF_Hotspot <- read.table("SKCM_BRAF_Hotspot_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_EGFR_mut <- read.table("LUAD_EGFR-mut_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_KRAS_hotspot <- read.table("LUAD_KRAS-hotspot_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LUAD_STK11_mut <- read.table("LUAD_STK11-mut_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_CIN <- read.table("STAD_CIN_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_GS <- read.table("STAD_GS_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
STAD_MSI <- read.table("STAD_MSI_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
ESCA_CIN <- read.table("ESCA_CIN_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
ESCA_ESCC <- read.table("ESCA_ESCC_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
READ_CIN <- read.table("READ_CIN_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LGG_IDHmut_codel <- read.table("LGG_IDHmut_codel_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
LGG_IDHmut_non_codel <- read.table("LGG_IDHmut_non_codel_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Differentiated <- read.table("OV_Differentiated_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Immunoreactive <- read.table("OV_Immunoreactive_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Mesenchymal <- read.table("OV_Mesenchymal_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
OV_Proliferative <- read.table("OV_Proliferative_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Basal <- read.table("BRCA_Basal_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Her2 <- read.table("BRCA_Her2_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_LumA <- read.table("BRCA_LumA_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_LumB <- read.table("BRCA_LumB_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
BRCA_Normal <- read.table("BRCA_Normal_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_CN_HIGH <- read.table("UCEC_CN_HIGH_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_CN_LOW <- read.table("UCEC_CN_LOW_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
UCEC_POLE <- read.table("UCEC_POLE_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)
PRAD_TMPRSS2_ERG_Fusion <- read.table("PRAD_TMPRSS2-ERG_Fusion_immune_score.txt",header=T,row.names = 1,sep = '\t',check.names = F)


#install.packages("robustbase")
#library("robustbase")

#getwd()
patient_annotation <- read.table("TMPRSS2-ERG-fusion_PRAD_annotation.txt",header=T,row.names=1,sep = '\t',check.names = F)
immune_score <- read.table("immune_infiltration_estimation_for_tcga.txt",header=T,row.names=1,sep = '\t',check.names = F)


#rownames(patient_annotation) <- patient_annotation[,1]
#patient_annotation <- patient_annotation[,-1]

#rownames(immune_score) <- immune_score[,1]
#immune_score <- immune_score[,-1]

pancancer_list <- list(HNSC_HPV_Neg,HNSC_HPV_Pos,SKCM_RAS_Hotspot,
                       SKCM_BRAF_Hotspot,LUAD_EGFR_mut,LUAD_KRAS_hotspot,
                       LUAD_STK11_mut,STAD_CIN,STAD_GS,STAD_MSI,ESCA_CIN,
                       ESCA_ESCC,READ_CIN,LGG_IDHmut_codel,LGG_IDHmut_non_codel,
                       OV_Differentiated,OV_Immunoreactive,OV_Mesenchymal,
                       OV_Proliferative,BRCA_Basal,BRCA_Her2,BRCA_LumA,BRCA_LumB,
                       BRCA_Normal,UCEC_CN_HIGH,UCEC_CN_LOW,UCEC_POLE,
                       PRAD_TMPRSS2_ERG_Fusion)
names(pancancer_list) <- c("HNSC_HPV_Neg","HNSC_HPV_Pos","SKCM_RAS_Hotspot",
                       "SKCM_BRAF_Hotspot","LUAD_EGFR_mut","LUAD_KRAS_hotspot",
                       "LUAD_STK11_mut","STAD_CIN","STAD_GS","STAD_MSI","ESCA_CIN",
                       "ESCA_ESCC","READ_CIN","LGG_IDHmut_codel","LGG_IDHmut_non_codel",
                       "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
                       "OV_Proliferative","BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB",
                       "BRCA_Normal","UCEC_CN_HIGH","UCEC_CN_LOW","UCEC_POLE",
                       "PRAD_TMPRSS2_ERG_Fusion")



immune_score <- immune_score[rownames(patient_annotation),]
immune_score <- as.data.frame(immune_score)
immune_score[1:5,1:5]

immune_score$"8q24"<-patient_annotation$"8q24"
immune_score <-na.omit(immune_score)
write.table(immune_score, "PRAD_TMPRSS2-ERG_Fusion_immune_score.txt",sep="\t", quote=F,row.names = T,col.names = T)
#write.table(data.frame(ID=rownames(immune_score),immune_score),"SKCM_immune_score.txt", sep="\t",quote=F,row.names=F,col.names = T)

wilcoxon_results <- list()

for (i in seq_along(cancer_18types_immune_score)) {
  # Current data frame
  df <- cancer_18types_immune_score[[i]]

  
  Amp <- grep("Amp",df$"8q24")
  WT <- grep("WT",df$"8q24")
  
  group <- c(rep(2,length(Amp)),rep(1,length(WT)))
  
  df<-df[,-ncol(df)]
  
  df<- t(df)
  df <- as.data.frame(df)
  
  # Data frame name (optional)
  #df_name <- deparse(substitute(pancancer_list[[i]]))
  
  pvalues <- sapply(1:nrow(df),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(df[i,])),group)
    p=wilcox.test(gene~group, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  log10_fdr=-log10(p.adjust(pvalues,method = "fdr"))
  
  # Calculate fold-change for each gene
  
  
  dataCon1=df[,c(which(group==1))]
  dataCon2=df[,c(which(group==2))]
  dataCon1<-as.data.frame(dataCon1)
  dataCon2<-as.data.frame(dataCon2)
  
  
  ### Robust weighted-mean idea (from robust volcano plot) to reduce outlier impact
  #log2foldChngCalc <- function (data, nSampG1, nSampG2) 
  #{
    #g1fold <- apply((data[, 1:nSampG1]), 1, weightedMean)
    #g2fold <- apply((data[, (nSampG1 + 1):(nSampG1 + nSampG2)]), 
                    #1, weightedMean)
    #foldChange <- log2(g1fold/g2fold)
    #return(foldChange)
  #}
  
  #Log2_foldChanges <- log2foldChngCalc(df,nSampG1 = ncol(dataCon2),nSampG2 = ncol(dataCon1))
  row_means_trimmed_2<- apply(dataCon2[, 1:ncol(dataCon2)], 1, function(x) mean(x))
  row_means_trimmed_1<- apply(dataCon1[, 1:ncol(dataCon1)], 1, function(x) mean(x))
  
  Log2_foldChanges=log2(row_means_trimmed_2/row_means_trimmed_1)
  
  #Log2_foldChanges=log2(rowMeans(dataCon2,trim=0.05)/rowMeans(dataCon1,trim=0.05))
  
  # Output results based on FDR threshold
  Log2_foldChanges<-as.data.frame(Log2_foldChanges)
  pvalues<-as.data.frame(pvalues)
  fdr<-as.data.frame(fdr)
  log10_fdr<-as.data.frame(log10_fdr)
  
  
  outRst<-data.frame(log2foldChange=Log2_foldChanges, pvalues=pvalues, fdr=fdr,log10_fdr=log10_fdr)
  rownames(outRst)=rownames(df)
  outRst=na.omit(outRst)
  wilcoxon_results[[i]] <- outRst
}

names(wilcoxon_results) <- c("HNSC_HPV_Neg","SKCM_RAS_Hotspot",
                             "SKCM_BRAF_Hotspot","LUAD_EGFR_mut","LUAD_KRAS_hotspot",
                             "LUAD_STK11_mut","STAD_CIN","STAD_GS","STAD_MSI","ESCA_CIN",
                             "ESCA_ESCC","READ_CIN","LGG_IDHmut_non_codel",
                             "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
                             "OV_Proliferative","BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB",
                             "BRCA_Normal","UCEC_CN_HIGH",
                             "PRAD_TMPRSS2_ERG_Fusion")

all_wilcoxon_results <- do.call(rbind,wilcoxon_results)
all_wilcoxon_results1 <- tibble::rownames_to_column(all_wilcoxon_results,"cancer_type+immune_score")
all_wilcoxon_results2 <- tidyr::separate(all_wilcoxon_results1,"cancer_type+immune_score",into = c("cancer_type","immune_score"),sep = "[.]")

write.table(all_wilcoxon_results2, "BCR_Richness+Shhonno_subtypes_all_wilcoxon_results.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


output_dir <- file.path("outputs", "immunomodulators_subtype_wilcoxon_result")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save each subtype result as a .txt file
for (i in seq_along(wilcoxon_results)) {
  df <- wilcoxon_results[[i]]

  file_name_all <- paste0("all_", names(wilcoxon_results)[i], ".txt")
  file_path_all <- file.path(output_dir, file_name_all)

  file_name_fdr <- paste0("0.05_", names(wilcoxon_results)[i], ".txt")
  file_path_fdr <- file.path(output_dir, file_name_fdr)

  write.table(df, file = file_path_all, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  fdrThres <- 0.05
  write.table(df[df$fdr <= fdrThres, ], file = file_path_fdr, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = TRUE)
}


















Amp <- grep("Amp",immune_score$"8q24")
WT <- grep("WT",immune_score$"8q24")

group <- c(rep(2,length(Amp)),rep(1,length(WT)))

immune_score<-immune_score[,-ncol(immune_score)]

immune_score<- t(immune_score)
immune_score <- as.data.frame(immune_score)




# Run the Wilcoxon rank-sum test for each gene

pvalues <- sapply(1:nrow(immune_score),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(immune_score[i,])),group)
  p=wilcox.test(gene~group, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")
log10_fdr=-log10(p.adjust(pvalues,method = "fdr"))

# Calculate fold-change for each gene


dataCon1=immune_score[,c(which(group==1))]
dataCon2=immune_score[,c(which(group==2))]
dataCon1<-as.matrix(dataCon1)
dataCon2<-as.matrix(dataCon2)
Log2_foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold
Log2_foldChanges<-as.data.frame(Log2_foldChanges)
pvalues<-as.data.frame(pvalues)
fdr<-as.data.frame(fdr)
log10_fdr<-as.data.frame(log10_fdr)



outRst<-data.frame(log2foldChange=Log2_foldChanges, pvalues=pvalues, fdr=fdr,log10_fdr=log10_fdr)
rownames(outRst)=rownames(immune_score)
outRst=na.omit(outRst)
write.table(outRst, "OV_WilcoxonTest1.txt",sep="\t", quote=F,row.names = T,col.names = T)
fdrThres=0.05
write.table(outRst[outRst$fdr<fdrThres,], "0.05_OV_WilcoxonTest1.txt",sep="\t", quote=F,row.names = T,col.names = T)





## Notes (not executed): general guidance about Wilcoxon and edgeR
# In cohort-level studies the sample size is typically large enough that strict
# parametric assumptions are often not required.
# As a non-parametric method, the Wilcoxon rank-sum test cannot adjust for
# confounders (e.g., sequencing depth). Consider proper normalization and batch
# correction before downstream testing.

if (FALSE) {
  library(edgeR)
  # Example: read count matrix and phenotype labels
  count <- read.table("TCGA-LUAD.htseq_counts.tsv", header = TRUE, row.names = 1)
  patient_ID <- read.table("normal_vs_EGFR_mut.txt", header = TRUE, row.names = 1)

  count1 <- t(count)
  count1 <- as.data.frame(count1)
  count1 <- count1[rownames(patient_ID), ]

  # Define groups (1 = control, 2 = case)
  group <- c(rep(1, 114), rep(2, 21))

  # Pre-processing: create DGEList
  dgelist <- DGEList(counts = count, group = group)

  keep <- filterByExpr(dgelist)

  # Manual filtering alternative
  # keep <- rowSums(cpm(dgelist) > 1) >= 2

  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

  # Normalization (TMM)
  dgelist_norm <- calcNormFactors(dgelist, method = "TMM")

  count_norm <- cpm(dgelist_norm)
  count_norm <- as.data.frame(count_norm)

  # Run the Wilcoxon rank-sum test for each gene
  pvalues <- sapply(1:nrow(count_norm), function(i) {
    data <- cbind.data.frame(gene = as.numeric(t(count_norm[i, ])), group)
    p <- wilcox.test(gene ~ group, data)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = "fdr")

  # Calculate fold-change for each gene
  dataCon1 <- count_norm[, c(which(group == 1))]
  dataCon2 <- count_norm[, c(which(group == 2))]
  foldChanges <- log2(rowMeans(dataCon2) / rowMeans(dataCon1))

  # Output results based on FDR threshold
  foldChanges <- as.data.frame(foldChanges)
  pvalues <- as.data.frame(pvalues)
  fdr <- as.data.frame(fdr)

  outRst <- data.frame(log2foldChange = foldChanges, pvalues = pvalues, FDR = fdr)
  rownames(outRst) <- rownames(count_norm)
  outRst <- na.omit(outRst)
  fdrThres <- 0.05
  write.table(outRst[outRst$FDR < fdrThres, ], "WilcoxonTest1.rst.txt",
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
}

