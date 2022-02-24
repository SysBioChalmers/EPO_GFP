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
library(here)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

# import samples table
  sam <- read.delim('/Users/rasools/drive/projects/EPO_GFP/EPO_GFP_github/data/samples_for_correlation_analysis.txt')
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

# method = "pearson"
#   # correlation analysis
#     corel <- matrix(nrow=nrow(y_gfp), ncol=4)
#     for(i in 1:nrow(y_gfp)) {
#       corel[i,1] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["estimate"]]
#       corel[i,2] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["p.value"]]
#       corel[i,3] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["estimate"]]
#       corel[i,4] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["p.value"]]
#     }
#     corel <- as.data.frame(corel)
#     colnames(corel) <- c("epo_coef",'epo_pval','gfp_coef',"gfp_pval")
#     rownames(corel) <- rownames(y_gfp)
#   
#   # finding significantly correlating genes
#     df <- corel
#     rownames(df) <- rownames(corel)
#     th = 0.5
#     ## filtering non significant genes based on pvalue lower than 0.05
#       df$epo_coef[which(df$epo_pval > 0.05)] = 0
#       df$gfp_coef[which(df$gfp_pval > 0.05)] = 0
#     ## filtering genes with mean of tpm lower than 10
#       df$epo_coef[which(rowMeans(tpm[,3:8]) < 10)] <- 0
#       df$gfp_coef[which(rowMeans(tpm[,9:16]) < 10)] <- 0
#     ## tagging high correlating genes
#       df$sig <- 'non'
#       df$sig[which(df$epo_coef > th)] <- 'epo+'
#       df$sig[which(df$epo_coef < -th)] <- 'epo-'
#       df$sig[which(df$gfp_coef > th)] <- 'gfp+'
#       df$sig[which(df$gfp_coef < -th)] <- 'gfp-'
#       df$sig[which(df$epo_coef > th & df$gfp_coef > th)] <- 'epo+gfp+'
#       df$sig[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'epo-gfp-'
#       df$sig[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'epo+gfp-'
#       df$sig[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'epo-gfp+'
#       df$prot <- 'non'
#       df$prot[which(df$epo_coef > th)] <- 'epo'
#       df$prot[which(df$epo_coef < -th)] <- 'epo'
#       df$prot[which(df$gfp_coef > th)] <- 'gfp'
#       df$prot[which(df$gfp_coef < -th)] <- 'gfp'
#       df$prot[which(df$epo_coef > th & df$gfp_coef > th)] <- 'both'
#       df$prot[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'both'
#       df$prot[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'both'
#       df$prot[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'both'
#     ## filtering non correlating genes
#       sig <- df[-which(df$prot == 'non'),]
#       sig <- merge(x=sig, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
#       rownames(sig) <- sig$Row.names
#       sig <- distinct(sig)
#       sig <- sig[!is.na(sig$external_gene_name),]
#       sig <- sig[!duplicated(sig$external_gene_name),]
#   
#   # saving table
#   #epo positive
#   write.table(sig[sig$sig == 'epo+' | sig$sig == 'epo+gfp-' | sig$sig == 'epo+gfp+', 'external_gene_name'], 
#               paste(getwd(),'/results/',"35samples_corel_pear_epo_pos_thr_",th,".txt",sep=""), 
#               quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
#   #epo negative
#   write.table(sig[sig$sig == 'epo-' | sig$sig == 'epo-gfp-' | sig$sig == 'epo-gfp+', 'external_gene_name'], 
#               paste(getwd(),'/results/',"35samples_corel_pear_epo_neg_thr_",th,".txt",sep=""), 
#               quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
#   #gfp positive
#   write.table(sig[sig$sig == 'gfp+' | sig$sig == 'epo+gfp+' | sig$sig == 'epo-gfp+', 'external_gene_name'], 
#               paste(getwd(),'/results/',"35samples_corel_pear_gfp_pos_thr_",th,".txt",sep=""), 
#               quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
#   #gfp negative
#   write.table(sig[sig$sig == 'gfp-' | sig$sig == 'epo+gfp-' | sig$sig == 'epo-gfp-', 'external_gene_name'], 
#               paste(getwd(),'/results/',"35samples_corel_pear_gfp_neg_thr_",th,".txt",sep=""), 
#               quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
#   
#   write.table(sig, paste(getwd(),'/results/',"35samples_corel_pearson_thr_0.5",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t') 
# 

# spearman analysis
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
  ## tagging high correlating genes
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
  #epo all
  write.table(sig[sig$prot == 'epo' | sig$prot == 'both', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_epo_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 

  #epo all
  write.table(sig[sig$prot == 'epo' | sig$prot == 'both', c('Row.names','epo_coef')], 
              paste(getwd(),'/results/',"35samples_corel_coef_spear_epo_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
    
  #epo positive
  write.table(sig[sig$sig == 'epo+' | sig$sig == 'epo+gfp-' | sig$sig == 'epo+gfp+', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_epo_pos_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
  #epo negative
  write.table(sig[sig$sig == 'epo-' | sig$sig == 'epo-gfp-' | sig$sig == 'epo-gfp+', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_epo_neg_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
  #gfp all
  write.table(sig[sig$prot == 'gfp' | sig$prot == 'both', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_gfp_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
  
  #gfp all
  write.table(sig[sig$prot == 'gfp' | sig$prot == 'both', c('Row.names','gfp_coef')], 
              paste(getwd(),'/results/',"35samples_corel_coef_spear_gfp_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
  
  #gfp positive
  write.table(sig[sig$sig == 'gfp+' | sig$sig == 'epo+gfp+' | sig$sig == 'epo-gfp+', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_gfp_pos_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
  #gfp negative
  write.table(sig[sig$sig == 'gfp-' | sig$sig == 'epo+gfp-' | sig$sig == 'epo-gfp-', 'external_gene_name'], 
              paste(getwd(),'/results/',"35samples_corel_spear_gfp_neg_thr_",th,".txt",sep=""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 
  
  write.table(sig, paste(getwd(),'/results/',"35samples_corel_spearman_thr_0.5",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t') 
  
  #------------------------------------------------------------------------
  tpm <- merge(x=tpm, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
  #------
  gene_name = "UBE2J1"
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
  colnames(cordata) <- c('gene','epoQP')
  ggscatter(cordata, x = 'gene', y = 'epoQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', TPM', sep=''), ylab = 'EPO production, pg/cell.day')
  ggsave(paste('/Users/rasools/Desktop/',gene_name,'_epo.png', sep = ""),width = 10, height = 7, units="cm")
  
  
  gene_name = "EPOR"
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
  colnames(cordata) <- c('gene','gfpQP')
  ggscatter(cordata, x = 'gene', y = 'gfpQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', TPM', sep=''), ylab = 'GFP production, Relative to GFPpoly')
  ggsave(paste('/Users/rasools/Desktop/',gene_name,'_gfp.png',sep = ""),width=10,height=7, units="cm")
  





