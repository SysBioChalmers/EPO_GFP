# import needed packages
  library(biomaRt)
  library(dplyr)
  library(reshape)
  library(tidyverse)
  library(viridis)
  library(ggpubr)
  library(piano)
  library(DESeq2)
  library(msigdb)
  library(snowfall)
  library(snow)
  library(hrbrthemes)
  library(packcircles)

gsc = read.gmt('c5.all.v7.1.symbols.gmt')
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
# import counts
counts <- read.delim("data/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))

# import samples table
nam <- colnames(counts)
samples <- data.frame(replicate = nam, sampleName = c(sapply(strsplit(nam, "_"),
                                                             "[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), 
                                                                                                                                                                   c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), 
                      growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), 
                                   c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), 
                                                                                    c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))
#---------------------------------------------------------------
#comparing EPO-producers with GFP-producers 
samples <- samples[samples$producibility == "Producer",]
counts <- counts[,rownames(samples)]
all(rownames(samples) == colnames(counts))
#########################EPO7,EPO8,EPOB9,EPOF21,EPOI2,EPOpoly
#---------------------------------------------------------------
resultsNames(ds)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ protein)
ds <- DESeq(dds)
#norcounts <- counts(ds_epo7, normalized=TRUE)

#---------------------------------------------------------------
res <- results(ds, contrast=c("protein", "EPO", "GFP"))
res <- res[ ! is.na(res$padj), ]
sig <- res
sig <- res[ which(res$padj < 0.05), ]
#sig <- sig[ which(abs(sig$log2FoldChange) > 0.58), ]
#sig <- sig[order(sig$log2FoldChange),]
sig <- as.data.frame(sig)
sig <- merge(x=sig, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
epoRiboSig <- read.delim('results/epoRiboSig.txt', header = FALSE )
data <- sig[which(sig$external_gene_name %in% epoRiboSig$V1),]
packing <- circleProgressiveLayout(data$log2FoldChange, sizetype='radius')
data <- cbind(data, packing)
dat.gg <- circleLayoutVertices(packing, npoints=50)
dat.gg$value <- rep(data$log2FoldChange, each=51)
ggplot() + 
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=value), colour = "black", alpha = 0.6) +
  geom_text(data = data, aes(x, y, size=5, label = external_gene_name)) +
  theme_void() + 
  scale_size_continuous(range = c(1,4)) +
  scale_fill_distiller(palette = "BuPu", direction = 1 ) +
  coord_equal()
ggsave('results/bubbleplot_epoRibosig.svg',width = 10, height = 10, units="cm")


