# Correlation analysis
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
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     

# other attached packages:
#  [1] piano_2.4.0       ggpubr_0.3.0      viridis_0.5.1    
#  [4] viridisLite_0.3.0 forcats_0.5.0     stringr_1.4.0    
#  [7] purrr_0.3.4       readr_1.3.1       tidyr_1.1.0      
# [10] tibble_3.0.1      ggplot2_3.3.0     tidyverse_1.3.0  
# [13] reshape_0.8.8     dplyr_0.8.5       biomaRt_2.44.0  
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***

library(biomaRt)
library(dplyr)
library(reshape)
library(tidyverse)
library(viridis)
library(ggpubr)
library(piano)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)

# import samples table
  sam <- read.delim(paste(getwd(),'/data/samples_for_correlation_analysis.txt',sep=''))
  rownames(sam) <- sam[,1]
  sam <- sam[,2:17]
  samples_epo <- sam[,grep('EPO',colnames(sam))]
  samples_epo['epoQP',] <- samples_epo['productivity',]
  samples_gfp <- sam[,grep('GFP',colnames(sam))]
  samples_gfp['gfpQP',] <- samples_gfp['productivity',]
  
# import TPMs
tpm <- read.delim(paste(getwd(),'/data/35samples_TPM.txt',sep=''))
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]
colnames(tpm) = rep(colnames(sam), times = c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2))

# preparation for correlation analysis
  y_epo <- tpm[,colnames(samples_epo)]
  y_gfp <- tpm[,colnames(samples_gfp)]

method = "spearman"
# correlation analysis
  corel <- matrix(nrow=nrow(y_gfp), ncol=4)
  for(i in 1:nrow(y_gfp)) {
    corel[i,1] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["estimate"]]
    corel[i,2] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["p.value"]]
    corel[i,3] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["estimate"]]
    corel[i,4] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["p.value"]]
  }
  corel <- as.data.frame(corel)
  colnames(corel) <- c("epo_coef",'epo_pval','gfp_coef',"gfp_pval")
  rownames(corel) <- rownames(y_gfp)

# finding significantly correlating genes
  df <- corel
  rownames(df) <- rownames(corel)
  th = 0.5
  ## filtering non significant genes based on pvalue lower than 0.05
    df$epo_coef[which(df$epo_pval > 0.05)] = 0
    df$gfp_coef[which(df$gfp_pval > 0.05)] = 0
  ## filtering genes with mean of tpm lower than 10
    df$epo_coef[which(rowMeans(tpm[,3:8]) < 10)] <- 0
    df$gfp_coef[which(rowMeans(tpm[,9:16]) < 10)] <- 0
  ## taging high correlating genes
    df$sig <- 'non'
    df$sig[which(df$epo_coef > th)] <- 'epo+'
    df$sig[which(df$epo_coef < -th)] <- 'epo-'
    df$sig[which(df$gfp_coef > th)] <- 'gfp+'
    df$sig[which(df$gfp_coef < -th)] <- 'gfp-'
    df$sig[which(df$epo_coef > th & df$gfp_coef > th)] <- 'epo+gfp+'
    df$sig[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'epo-gfp-'
    df$sig[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'epo+gfp-'
    df$sig[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'epo-gfp+'
    df$prot <- 'non'
    df$prot[which(df$epo_coef > th)] <- 'epo'
    df$prot[which(df$epo_coef < -th)] <- 'epo'
    df$prot[which(df$gfp_coef > th)] <- 'gfp'
    df$prot[which(df$gfp_coef < -th)] <- 'gfp'
    df$prot[which(df$epo_coef > th & df$gfp_coef > th)] <- 'both'
    df$prot[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'both'
    df$prot[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'both'
    df$prot[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'both'
  ## filtering non correlating genes
    sig <- df[-which(df$prot == 'non'),]
    sig <- merge(x=sig, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
    rownames(sig) <- sig$Row.names
    sig <- distinct(sig)
    sig <- sig[!is.na(sig$external_gene_name),]
    sig <- sig[!duplicated(sig$external_gene_name),]

# saving table
write.table(sig, paste(getwd(),'/results/',"35samples_corel_spearman_epo_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 







