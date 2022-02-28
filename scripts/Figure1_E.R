library(dplyr)
library(ggplot2)
library(reshape2)

# sessionInfo()
# R version 4.0.0 (2020-04-24)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

# locale:
# [1] C/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils    
# [7] datasets  methods   base     

# other attached packages:
#  [1] snowfall_1.84-6.1           snow_0.4-3                 
#  [3] biomaRt_2.44.0              piano_2.4.0                
#  [5] ggplot2_3.3.0               DESeq2_1.28.1              
#  [7] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
#  [9] matrixStats_0.56.0          Biobase_2.48.0             
# [11] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
# [13] IRanges_2.22.2              S4Vectors_0.26.1           
# [15] BiocGenerics_0.34.0  

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
samples$epoL2TPM <- log(as.vector(t(TPM["EPO_HPC4_splitGFP",])),2)
samples$epoL2TPM[which(!is.finite(samples$epoL2TPM))] <- 0
samples$gfpL2TPM <- log(as.vector(t(TPM["pD2529_GFP_extra",])),2)
samples$gfpL2TPM[which(!is.finite(samples$gfpL2TPM))] <- 0
# check
all(colnames(TPM) == colnames(counts))
#---------------------------------------------------------------
# comparing EPO-producers with F21 control cell-line
samples_epo <- samples[samples$cellType == "HEK293freestyle",]
#--------------------------------------------------------------- EPO scatter plot
data_epo <- as.data.frame(unique(samples_epo$sampleName))
data_epo$productivity <- c(0,3.67,4.05,3.35,13.9,2.87,2.35)
data_epo$log2epoTpm_lower <- c(0,14.50062,15.33458,14.38587,15.19843,13.88685,13.01370)
data_epo$log2epoTpm_higher <- c(0,14.86886,15.44983,14.51067,15.36507,14.19705,13.51256)
data_epo$TPM <- rowMeans(data_epo[,3:4])
data_epo$growth <- c(0.0257,0.0227,0.0200,0.0241,0.0240,0.0223,0.0257)
colnames(data_epo) <- c('cellType','productivity','log2epoTpm_lower','log2epoTpm_higher','TPM','growth')
data_epo[5,5] <- mean(samples[11:13,9])
data_epo <- data_epo[2:7,]

ggplot(data=data_epo, aes(y=productivity, x=TPM), size=7, shape=16) +
  geom_point() +
  #geom_smooth(method=lm) +
  geom_errorbarh(aes(xmin = log2epoTpm_lower, xmax = log2epoTpm_higher)) +
  xlab(paste0(bquote("EPO mRNA expression, TPM"))) +
  ylab(paste0("Productivity (pg/cell.day)")) + 
  geom_text(aes(label=cellType),hjust=1, vjust=-1) + 
  theme(
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 16),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_blank(),
  axis.text.x=element_text(colour="black",size = 18,family="Arial"),
  axis.text.y=element_text(colour="black",size = 18,family="Arial"),
  axis.title=element_text(colour="black",size = 20),
  axis.ticks=element_line(colour="black"),
  plot.margin=unit(c(3,1,1,1),"line"),
  axis.line=element_line(arrow = arrow()))
#ggsave(('./results/tpm_productivity_scatter_EPO.svg'), width=14, height=16, units="cm")






