library(pheatmap)
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)
library(biomaRt)
library(readr)
library(data.table)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(DESeq2)
library(gg3D)

# import TPMs
TPM <- read.delim('./data/35samples_TPM.txt')
rownames(TPM) <- TPM[,1]
TPM <- TPM[,-1]

ds<-TPM[rowSums(TPM) > 1,]
#QC on normalised counts
rld <- log2(ds+1)
#PCA
pca <- prcomp(t(rld), scale=T)
pcasum <- as.data.frame(apply(pca$x, 2, var))
colnames(pcasum)<-"Variance"
pcasum$Component <- rownames(pcasum)
pcasum <-within(pcasum, Variance_Fraction <- 100*(Variance/sum(Variance)))
pcasum <-within(pcasum, Variance_csum <- signif(cumsum(Variance_Fraction), digits = 2))
pcasum <-within(pcasum, Order <- as.numeric(substring(pcasum$Component, 3)))
pcasum$Component <- factor(pcasum$Component, levels = pcasum$Component[order(pcasum$Order)])

sampinfo<-samples
scores <- merge(sampinfo, as.data.frame(pca$x[,1:3]), by.x="row.names",by.y="row.names", all.x = FALSE, all.y = TRUE)
scores$condition<-factor(scores$sampleName)
scores$protein<-factor(scores$protein)
nb.cols <- length(unique(scores$sampleName))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

#PCA1vsPCA2
pc1.2 <- ggplot(scores, aes(x=PC1, y=PC2) )+geom_point(aes(color=sampleName, shape=protein))
pc1.2 <- pc1.2 + scale_x_continuous(paste('PC1 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC1"], digits = 2), ' %)', sep=""))+ scale_y_continuous(paste('PC2 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC2"], digits = 2), ' %)', sep="") ) + ggtitle("log2(TPM+1) PCA plot")
pc1.2 <- pc1.2 + theme(plot.title=element_text(size = 20,face = "bold"), axis.text.y = element_text(angle = 30, hjust = 0.7, size = 12), axis.text.x = element_text(angle = 30, hjust = 0.7, size = 12), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
pc1.2 <- pc1.2 + scale_color_manual(values = mycolors)
pc1.2 <- pc1.2 + geom_label_repel(aes(label = replicate), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') 
ggsave(file=paste('results/',"PCA1vsPCA2_log2TPM.pdf", sep=""), plot=pc1.2, width = 14, height = 7)
#PCA1vsPCA3
pc1.3 <- ggplot(scores, aes(x=PC1, y=PC3) )+geom_point(aes(color=sampleName, shape=protein))
pc1.3 <- pc1.3 + scale_x_continuous(paste('PC1 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC1"], digits = 2), ' %)', sep=""))+ scale_y_continuous(paste('PC3 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC3"], digits = 2), ' %)', sep="") ) + ggtitle("log2(TPM+1) PCA plot")
pc1.3 <- pc1.3 + theme(plot.title=element_text(size = 20,face = "bold"), axis.text.y = element_text(angle = 30, hjust = 0.7, size = 12), axis.text.x = element_text(angle = 30, hjust = 0.7, size = 12), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
pc1.3 <- pc1.3 + scale_color_manual(values = mycolors)
pc1.3 <- pc1.3 + geom_label_repel(aes(label = replicate), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') 
ggsave(file=paste('results/',"PCA1vsPCA3_log2TPM.pdf", sep=""), plot=pc1.3, width = 14, height = 7)
#PCA2vsPCA3
pc2.3 <- ggplot(scores, aes(x=PC2, y=PC3) )+geom_point(aes(color=sampleName, shape=protein))
pc2.3 <- pc2.3 + scale_x_continuous(paste('PC2 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC2"], digits = 2), ' %)', sep=""))+ scale_y_continuous(paste('PC3 (Variance ', signif(pcasum$Variance_Fraction[pcasum$Component=="PC3"], digits = 2), ' %)', sep="") ) + ggtitle("log2(TPM+1) PCA plot")
pc2.3 <- pc2.3 + theme(plot.title=element_text(size = 20,face = "bold"), axis.text.y = element_text(angle = 30, hjust = 0.7, size = 12), axis.text.x = element_text(angle = 30, hjust = 0.7, size = 12), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16))
pc2.3 <- pc2.3 + scale_color_manual(values = mycolors)
pc2.3 <- pc2.3 + geom_label_repel(aes(label = replicate), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50') 
ggsave(file=paste('results/',"PCA2vsPCA3_log2TPM.pdf", sep=""), plot=pc2.3, width = 14, height = 7)



#Correlation clustering
cormat<-cor(rld,  method="spearman")
sampleDists <- as.dist(1-cor(rld,  method="spearman"))
saminfo<-subset(sampinfo,select=c(sampleName,protein,producibility,cellType))
saminfo$sampleName<-as.character(saminfo$sampleName)
saminfo$protein <-as.character(saminfo$protein)
saminfo$producibility<-as.character(saminfo$producibility)
saminfo$cellType<-as.character(saminfo$cellType)

sample_colors = list(
  "sampleName" = c("CTRL293F"=mycolors[[1]],"CTRL293Free"=mycolors[[2]], "EPO7"=mycolors[[3]],"EPO8"=mycolors[[4]],"EPOB9"=mycolors[[5]],"EPOF21"=mycolors[[6]],"EPOI2"=mycolors[[7]],"EPOpoly"=mycolors[[8]],
                  "GFP1"=mycolors[[9]],"GFP25"=mycolors[[10]],"GFP26"=mycolors[[11]], "GFP27"=mycolors[[12]],"GFP28"=mycolors[[13]], "GFP29"=mycolors[[14]], "GFP3"=mycolors[[15]],"GFPpoly"=mycolors[[16]]),
  
  "protein" = c("EPO"="#f7f7f7","EPO_Control"="#cccccc", "GFP"="#969696","GFP_Control"="#525252"),
  
  "producibility"=c("Control"="#c994c7","Producer"="#dd1c77"),
  
  "cellType" = c("HEK293F"="#99d8c9","HEK293freestyle"="#2ca25f")
  )


pheatmap(cormat,cluster_rows=T, cluster_cols=T, 
         annotation_col=saminfo, 
         annotation_colors=sample_colors, 
         clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists, 
         show_rownames = F, show_colnames = T,annotation_legend = T,
         filename =paste0('results/',"CorrClust_log2TPM.pdf.pdf"),width = 12, height = 10 )