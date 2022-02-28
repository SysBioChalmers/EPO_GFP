library(DESeq2)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(reshape)

# import counts
counts <- read.delim("data/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
#counts$EPOF21_I_6 <- NULL
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# import samples table
nam <- colnames(counts)
samples <- data.frame(replicate = nam, sampleName = c(sapply(strsplit(nam, "_"),
    "[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), 
    c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), 
    growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), 
    c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), 
    c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), 
    c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))

#---------------------------------------------------------------
# import TPMs
TPM <- read.delim("data/35samples_TPM.txt")
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
colnames(TPM) <- gsub("abundance.","", colnames(TPM))
TPM <- TPM[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
samples$epoL2TPM <- log(as.vector(t(TPM["EPO_HPC4_splitGFP",])),2)
samples$epoL2TPM[which(!is.finite(samples$epoL2TPM))] <- 0
# check
all(colnames(TPM) == colnames(counts))
#---------------------------------------------------------------
# comparing EPO-producers with f-free control cell-line
samples <- samples[samples$protein == "EPO",]
counts <- counts[,rownames(samples)]
TPM <- TPM[,rownames(samples)]
samples$contrast = rep(c("cntrl","cntrl","cntrl","case","cntrl","cntrl"), c(2,2,2,3,2,2))
all(rownames(samples) == colnames(counts))
all(colnames(TPM) == colnames(counts))
#########################EPO7,EPO8,EPOB9,EPOF21,EPOI2,EPOpoly
#---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ contrast)
ds <- DESeq(dds)
norcounts <- counts(ds, normalized=TRUE)

#---------------------------------------------------------------
resultsNames(ds)
res <- results(ds, contrast=c("contrast", "case", "cntrl"))
res <- res[ ! is.na(res$padj), ]
#---------------------------------------------------------------
sig <- res
sig <- as.data.frame(sig[(abs(sig$log2FoldChange) > 1 & sig$padj < 0.05), ])
#---------------------------------------------------------------
detpm <- TPM[which(rownames(TPM) %in% rownames(sig)),]
detpm <- detpm[rowMeans(detpm) > 10,]
sig <- sig[which(rownames(sig) %in% rownames(detpm)),]
sig <- sig[order(sig$log2FoldChange),]
#--------------------------------------------------------------
denor <- norcounts[which(rownames(norcounts) %in% rownames(detpm)),]
denor <- denor[match(rownames(sig), rownames(denor)),]

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)
denor <- merge(x=denor, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
#detpm <- detpm[order(match(rownames(sig), detpm$Row.names)),]
denor <- denor[match(rownames(sig), denor$Row.names),]

denor$Row.names <- NULL

# de <- as.data.frame(rowMeans(detpm[,1:2]))
# de[,2] <- rowMeans(detpm[,3:4])
# de[,3] <- rowMeans(detpm[,5:6])
# de[,4] <- rowMeans(detpm[,7:9])
# de[,5] <- rowMeans(detpm[,10:11])
# de[,6] <- rowMeans(detpm[,12:13])
# de[,7] <- detpm[,14]
# colnames(de) <- c("EPO7","EPO8","EPOB9","EPOF21","EPOI2","EPOpoly","id")

df <- melt(denor, id = "external_gene_name")
colnames(df) <- c("Genes","cellType","norcounts")
df$EPO_producers <- as.character('Low EPO-producers')
df[grep('EPOF21',df$cellType),4] <- as.character('EPOF21')
df <- as.data.frame(df)
df$Genes <-  factor(unique(df$Genes), levels = unique(df$Genes))

ggplot(df, aes(x=Genes, y=log(norcounts+1,2))) +
  geom_boxplot(aes(colour=EPO_producers), outlier.shape=NA) +
  scale_colour_manual(values = c("Black", "forestgreen")) + 
  #geom_boxplot(data=df[df$EPO_producers=='low',] , aes(color = shape), colour = 'darkred') +
  #geom_boxplot(data=df[df$EPO_producers=='high',] , aes(color = shape), colour = 'Black') +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  panel.background = element_blank(),
  panel.border=element_rect(fill=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_blank(),
  axis.text.x=element_text(colour="black",size = 11,family="Arial",angle = 30, hjust = 1),
  axis.text.y=element_text(colour="black",size = 11,family="Arial"),
  axis.title=element_text(colour="black",size = 11),
  axis.ticks=element_line(colour="black"),
  plot.margin=unit(c(1,1,1,1),"line")
  ) +
ggsave(('results/box_plot.svg'), width=27, height=10, units="cm")

