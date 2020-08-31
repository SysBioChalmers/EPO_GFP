# PCA: all clones
# 
# sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils    
# [7] datasets  methods   base     

# other attached packages:
#  [1] ggplot2_3.3.0               DESeq2_1.28.1              
#  [3] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
#  [5] matrixStats_0.56.0          Biobase_2.48.0             
#  [7] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
#  [9] IRanges_2.22.2              S4Vectors_0.26.1           
# [11] BiocGenerics_0.34.0 
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
pathToEpoGfp <- '/path/to/local/epo_gfp_repo/'
  
library(DESeq2)
library(ggplot2)

# import counts
counts <- read.delim(paste(pathToEpoGfp,'data/35samples_counts.txt',sep=''))
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
#counts$EPOF21_I_6 <- NULL
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# import samples table
kallisto.dir <- list.files(path = paste(pathToEpoGfp,'data/kallisto_outputs/',sep=''))
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))
#---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~sampleName)
colnames(dds) <- colnames(counts)
ds <- DESeq(dds)
# Apply regularized-log transform to counts
rld <- rlogTransformation(ds)
# Principal component analysis
#png('/Users/rasools/Desktop/pca_plot.png')
pcaData <- plotPCA(rld, intgroup=c("protein","sampleName"), returnData=TRUE )
percentVar <- round(100 * attr(pcaData, "percentVar"))
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
ggplot(pcaData, aes(PC1, PC2, color=protein, label=samples$sampleName)) +
  geom_point(size=8) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme +
  #geom_text(aes(label=paste(samples$condition,"-",samples$growth,"-",samples$productiviy), size=0.5,check_overlap = TRUE,hjust=0, vjust=-2))
  geom_text(size=4,check_overlap = TRUE,hjust=0, vjust=-2) +
  ggsave(paste(pathToEpoGfp,'results/PCA_EPO.svg',sep=''), width=35, height=25, units="cm", dpi=300)
dev.off()
