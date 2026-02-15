####===Figure 2===####
rm (list=ls())
gc()


####===Fig. 2A===####
library(circlize)
library(RColorBrewer)
library(png)
library(graphics)
library(circlize)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(ggplot2)

dir <- data.table::fread('data/dir/dir.csv')
circos.clear()
Set3 <- brewer.pal(12,"Set3")
Set2 <- brewer.pal(8,"Set2")
col <- c(Set3,Set2)

sectors = dir$dataset
names(sectors) <- sectors

dataset <- data.frame(dataset = sectors)
rownames(dataset) <- sectors

col_cancer <- c(col[1], rep(col[2],2), col[3],col[4],rep(col[5],11),col[6],col[7],col[8],rep(col[9],2),col[10],col[11],rep(col[12],4),col[13],rep(col[14],2),
                rep(col[15],2),col[16],col[19])
names(col_cancer) <- sectors



image = 'data/png/venn.png'
image = as.raster(readPNG(image))

circos.clear()
circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.01, 0.01))
circos.heatmap(dataset, split = dataset$dataset, col = col_cancer, track.height = 0.02,rownames.side = 'outside')

circos.track(ylim = c(0, 1),track.height=0.15,bg.border = "#EEEEEE", panel.fun = function(x, y) {
  circos.raster(image, CELL_META$xcenter, CELL_META$ycenter, 
                width = "0.7cm",
                facing = "inside")
  
})

circos.track(ylim = c(0, 1), sectors = sectors,
             bg.col = "#8190A5", bg.border = "#EEEEEE" , track.height = 0.25)

circos.trackText(x = rep(0.5, 34), y = rep(0.5, 34),
                 labels = paste0(rep('G',34), 1:34),
                 cex = 0.5, sectors = sectors, col = "white", font = 2, facing = "clockwise",
                 niceFacing=T)

draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
            rou1 = 0.25, col = "#4D6381", border = "#EEEEEE")

text(0,0,"Stem.sig",col = 'white',cex = 2,font = 2)


cancer <- unique(col_cancer)
names(cancer) <- unique(dir$cancer)
lgd_cancer = Legend(title = "Cancer", at = names(cancer),
                 legend_gp = gpar(fill = cancer))
draw(lgd_cancer, x = unit(1.2, "snpc"), just = "left")



####===Fig. 2B===####
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

load('data/sig/Stem.Sig.Rdata')

sig <- Stem.Sig

sig_ID <- bitr(sig,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")

eReac <- enrichPathway(gene = sig_ID$ENTREZID,
                       organism = 'human',
                       pvalueCutoff = 0.05)

bardata <- eReac@result[1:20,]
bardata <- bardata[order(bardata$qvalue),]
bardata <- bardata[order(bardata$Description),]
bardata$Description <- factor(bardata$Description,levels=rev(bardata$Description))
barplot <- ggplot(bardata,aes(y=Description,x=Count))+geom_bar(stat = "identity",aes(fill=qvalue))+scale_fill_gradientn(colours =  c("#4D6381", "#BAC2CC"))+theme_bw()

eReacx <- setReadable(eReac, 'org.Hs.eg.db', 'ENTREZID')

cnetplot <-  cnetplot(eReacx, circular = TRUE, colorEdge = TRUE,categorySize="pvalue",showCategory = 20,layout = 'kk') 

barplot 
cnetplot







