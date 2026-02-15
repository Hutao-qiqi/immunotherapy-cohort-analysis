rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(colourpicker)
  library(ggplot2)
  library(scales)
  library(ggprism)
})


###########################   Data import & preprocessing   ###########################
# This script expects input TXT files to be in the same folder as the script
# (or you can run it after setting the working directory to your input folder).

# Data import
#library(openxlsx)
# .txt格式导入

# Input table
df <- read.table("18_cancer_type+subtype_immune_pathway_wilcox_result.txt",header = T,sep = '\t',check.names = F)

# 
df$cancer_type<-as.character(df$cancer_type)
df$immune_score<-as.character(df$immune_score)
#ggplot()+geom_point(data=a,aes(x=s1,y=s2,size=s2))+
  #scale_x_discrete(limits=as.character(c(1,4,3,10,5,2,7,9,8,6)))



### Reduced version (subset of cancer types)
df$cancer_type=factor(df$cancer_type,levels =c(
  "SKCM_All","SKCM_BRAF_Hotspot","SKCM_RAS_Hotspot",
  "STAD_All","STAD_CIN","STAD_MSI",
  "HNSC_All","HNSC_HPV-","HNSC_HPV+",
  
  "LUAD_All","LUAD_KRAS_hotspot","PAAD_All","LUSC_All",
  
  "ESCA_All","ESCA_ESCC","ESCA_CIN","LGG_All","LIHC_All","UCEC_All",
  "OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal",
  "BRCA_LumA","BRCA_Basal","BRCA_Her2","BRCA_Normal","BRCA_LumB",
  "BLCA_All",
  
  "READ_All","READ_CIN","UVM_All","PRAD_TMPRSS2_ERG_Fusion","PRAD_All") )







#### Full version (all cancer types)

df$cancer_type= factor(df$cancer_type,levels =c(
  "STAD_All","STAD_CIN","STAD_MSI","STAD_GS","SKCM_All","SKCM_BRAF_Hotspot","SKCM_RAS_Hotspot",
  "LUAD_All","LUAD_KRAS_Hotspot","LUAD_EGFR_mut","LUAD_STK11_mut","HNSC_All","HNSC_HPV-","HNSC_HPV+","PAAD_All",
  "ESCA_All","ESCA_ESCC","ESCA_CIN","LUSC_All","LGG_All","LGG_IDHmut_non_codel","LGG_IDHmut_codel","LIHC_All","UCS_All","BRCA_All","BRCA_LumA","BRCA_Basal","BRCA_Her2","BRCA_Normal","BRCA_LumB","BLCA_All",
  "OV_All","OV_Differentiated","OV_Immunoreactive","OV_Mesenchymal","OV_Proliferative",
  "UCEC_All","UCEC_CN_LOW","UCEC_POLE","UCEC_CN_HIGH",
  
  
  "READ_All","READ_CIN","CHOL_All","UVM_All","PRAD_TMPRSS2_ERG_Fusion","PRAD_All") )



df$gene_symbol=factor(df$gene_symbol,levels =c("BIOCARTA_CTL_PATHWAY","PID_CD8_TCR_PATHWAY","PID_TCR_PATHWAY",
                                               "BIOCARTA_TCR_PATHWAY","KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY","PID_CD8_TCR_DOWNSTREAM_PATHWAY","BIOCARTA_BCR_PATHWAY",
                                               "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","BIOCARTA_IFNG_PATHWAY","PID_IFNG_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE","REACTOME_INTERFERON_GAMMA_SIGNALING",
                                               "BIOCARTA_IFNA_PATHWAY","HALLMARK_INTERFERON_ALPHA_RESPONSE","REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                                              "HALLMARK_INFLAMMATORY_RESPONSE","REACTOME_TOLL_LIKE_RECEPTOR_CASCADES","REACTOME_TOLL_LIKE_RECEPTOR_9_TLR9_CASCADE","REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE","KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","REACTOME_ADAPTIVE_IMMUNE_SYSTEM",
                                               "HALLMARK_COMPLEMENT","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING") )







####完整版###therapeutic_targets_顺序

df$cancer_type= factor(df$cancer_type,levels =c("PRAD_All","PRAD_TMPRSS2_ERG_Fusion","LIHC_All","SKCM_All","SKCM_BRAF_Hotspot","SKCM_RAS_Hotspot",
                                                "STAD_All","STAD_CIN","STAD_MSI","STAD_GS","HNSC_All","HNSC_HPV-","HNSC_HPV+","LGG_All","LGG_IDHmut_non_codel","LGG_IDHmut_codel","BLCA_All",
                                                "BRCA_All","BRCA_LumA","BRCA_Basal","BRCA_Her2","BRCA_Normal","BRCA_LumB",
                                                "OV_All","OV_Differentiated","OV_Immunoreactive","OV_Proliferative","OV_Mesenchymal",
                                                "PAAD_All","LUSC_All","ESCA_All","ESCA_ESCC","ESCA_CIN","UCEC_All","UCEC_CN_HIGH","LUAD_All","LUAD_KRAS_hotspot","LUAD_EGFR_mut","LUAD_STK11_mut",
                                                "UCS_All",
                                                "READ_All","READ_CIN","CHOL_All","UVM_All") )






df$gene_symbol=factor(df$gene_symbol,levels =c("CX3CL1","CXCL10","CXCL9",
                                               "CCL2","CCL5","MICB","MICA",
                                               "HLA-C","HLA-B","HLA-A","HLA-DPB1",
                                               "HLA-DPA1","HLA-DQB1","HLA-DQA1",
                                               "HLA-DRB5","HLA-DRB1","HLA-DRA","CD40",
                                               "CD27","GZMA","PRF1","ARG1","VEGFA","PVR",
                                               "KIR2DL3","KIR2DL1",
                                               "IL10","BTLA") )



df$immune_score=factor(df$immune_score,levels =c("AQP3","TP63","SNORA74A","PCSK9","LDLR",
                                               "RPL18","RPL36","SERPINB5","TYSND1","CNKSR3","PCAT1",
                                               "PHB2",
                                               "RPS2","OLA1","TTLL12") )




df$immune_cell=factor(df$immune_cell,levels =c("B cell_QUANTISEQ","B cell_MCPCOUNTER","B cell_EPIC","B cell plasma_CIBERSORT-ABS","B cell naive_CIBERSORT-ABS","B cell memory_CIBERSORT-ABS",
                                               "T cell CD8+_MCPCOUNTER","T cell CD8+_QUANTISEQ","T cell CD8+_CIBERSORT-ABS","T cell CD8+_EPIC",
                                               "T cell_MCPCOUNTER","T cell follicular helper_CIBERSORT-ABS","T cell regulatory (Tregs)_CIBERSORT-ABS","T cell regulatory (Tregs)_QUANTISEQ",
                                               "T cell CD4+ (non-regulatory)_QUANTISEQ","T cell CD4+_EPIC","T cell CD4+ memory resting_CIBERSORT-ABS",
                                               "Myeloid dendritic cell_MCPCOUNTER","Myeloid dendritic cell resting_CIBERSORT-ABS","Myeloid dendritic cell_QUANTISEQ","Myeloid dendritic cell activated_CIBERSORT-ABS",
                                               "NK cell_MCPCOUNTER","NK cell activated_CIBERSORT-ABS","NK cell resting_CIBERSORT-ABS","NK cell_EPIC","NK cell_QUANTISEQ",
                                               "Monocyte_MCPCOUNTER","Monocyte_CIBERSORT-ABS","Macrophage/Monocyte_MCPCOUNTER","Macrophage_EPIC",
                                               "Cancer associated fibroblast_MCPCOUNTER","Cancer associated fibroblast_EPIC") )



####颜色梯度为连续值的气泡热图
ggplot()+
  #geom_point(size=df$"Log10(q-value)",colour = "black",shape = 21, stroke = 2)+
  
  geom_point(data=df,aes(x=cancer_type,y=gene_symbol,size=`Log10(q-value)`,
                 fill=`Log2FC`),color="grey",shape=21,stroke = 0.01)+
  #scale_shape_manual(values=c(21,21))+
  #scale_x_discrete(limits=as.character(c(ESCA_All,OV_All,UCS_All,BRCA_All,LIHC_All,UVM_All,
                                         #STAD_All,PAAD_All,HNSC_All,LUAD_All,UCEC_All,LUSC_All,
   
  scale_fill_gradientn(colours = c("darkblue","blue","white","red","#CC0000"),
                       values = c(0,0.35,0.7,0.85,1),
                       breaks = c(0.1,0,-0.1,-0.2))+ 
  theme_bw()+
  theme(axis.text=element_text(family="sans",color="black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1))+
  
  #scale_fill_gradient2(low="blue",mid = "white",high="red",midpoint =0 )+
  #scale_colour_gradientn(colors = scales::viridis_pal()(10))+
  
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(tiltle="-Log10(q-value)",order=1),
         color=guide_legend(tiltle="Log2FC",order=2))

c("#FA8B0D", "#FFB875", "#FF9D42", "#FFA530", "#FABD61")
c("#FF6633","#FFCC66","#FFFF66","#FFFF99","#FFFFCC","#f6fafd")
###therapeutic_targets
#c("#f6fafd", "#c8dfef", "#6fa6d1", "#2c49a2")
####颜色梯度为离散值的气泡热图
ggplot()+
  #geom_point(size=df$"Log10(q-value)",colour = "black",shape = 21, stroke = 2)+
  geom_point(data=df[which(df$Neg_Pos=="Pos"),],aes(x=cancer_type,y=immune_score,size=rho,
                                                    color=FDR),shape=16)+
  geom_point(data=df[which(df$Neg_Pos=="Neg"),],aes(x=cancer_type,y=immune_score,size=rho,
                                                    fill=FDR),color="white",shape=21,stroke = 0.01)+
  scale_fill_manual(values=c("#2c49a2","#6699CC","#6fa6d1", "#c8dfef","#f6fafd"))+
  scale_color_manual(values=c("#FF6633","#FFA530","#FFCC66","#FFFF66","#FFFF99","#f6fafd"))+
  theme_bw()+
  theme(axis.text=element_text(family="sans",color="black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(tiltle="Spearman_ρ",order=1),
         fill=guide_legend(tiltle="Negtive_correlation_FDR_q-value",order=2),
         color=guide_legend(tiltle="Positive_correlation_FDR_q-value",order=3))

c("#f6fafd","#6fa6d1","#6699CC","#3333CC","#333399","#333366")
c("#f6fafd", "#7E75FF", "#6C5DF5", "#3424E3", "#1400CC", "#0D0080B8")
c("#f6fafd","#FFFF99","#FF99CC","#FF6699","#CC0066","#CC0000")
c("#f6fafd","#FFE0E0", "#FF8F8F", "#F23030", "#DE0404", "#B80000")
###immune_cells




####颜色梯度为离散值的气泡热图
ggplot()+
  #geom_point(size=df$"Log10(q-value)",colour = "black",shape = 21, stroke = 2)+
  
  geom_point(data=df,aes(x=cancer_type,y=immune_cell,size=`Log10(q-value)`,
                 fill=`Log2FC`),color="white",shape=21,stroke = 0.01)+
  geom_point(data=df[which(df$Neg_Pos=="Neg"),],aes(x=cancer_type,y=immune_cell,size=`Log10(q-value)`,
                 color=`Log2FC`),shape=16)+
  scale_fill_manual(values=c("#f6fafd","#FFE0E0", "#FF8F8F", "#F23030", "#DE0404", "#B80000"))+
  scale_color_manual(values=c("#f6fafd", "#7E75FF", "#6C5DF5", "#3424E3", "#1400CC", "#0D0080B8"))+
  theme_bw()+
  theme(axis.text=element_text(family="sans",color="black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(tiltle="-Log10(q-value)",order=1),
         fill=guide_legend(tiltle="Negtive Log2FC",order=2),
         color=guide_legend(tiltle="Positive Log2FC",order=3))
        
        
  
  



