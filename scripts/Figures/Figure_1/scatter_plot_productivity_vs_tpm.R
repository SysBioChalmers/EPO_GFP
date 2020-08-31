# boxplots: TPM vs. Productivity
# 
# sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

# Random number generation:
#  RNG:     Mersenne-Twister 
#  Normal:  Inversion 
#  Sample:  Rounding 
 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     

# other attached packages:
# [1] reshape2_1.4.4 dplyr_0.8.5    gdtools_0.2.2  ggplot2_3.3.0 
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
pathToEpoGfp <- '/path/to/local/epo_gfp_repo/'

library(dplyr)
library(ggplot2)
library(reshape2)

counts <- read.delim(paste(pathToEpoGfp,'data/35samples_counts.txt',sep=''))
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]

# import samples table
kallisto.dir <- list.files(path = paste(pathToEpoGfp,'data/kallisto_outputs/',sep=''))
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))
#---------------------------------------------------------------
# import TPMs
TPM <- read.delim(paste(pathToEpoGfp,'data/35samples_TPM.txt',sep=''))
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
colnames(TPM) <- gsub("abundance.","", colnames(TPM))
TPM <- TPM[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
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
  axis.line=element_line(arrow = arrow())) +
  ggsave(paste(pathToEpoGfp,'results/tpm_productivity_scatter_EPO.svg',sep=''), width=14, height=16, units="cm")
dev.off() 

