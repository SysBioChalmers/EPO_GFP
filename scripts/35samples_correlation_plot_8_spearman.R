# import needed packages
library(biomaRt)
library(dplyr)
library(reshape)
library(tidyverse)
library(viridis)
library(ggpubr)
library(piano)

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

#------------------------------------------------------------------------
tpm <- merge(x=tpm, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
#------
gene_name = "SBDS"
gene = which(tpm$external_gene_name == gene_name)
cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
colnames(cordata) <- c('gene','epoQP')
ggscatter(cordata, x = 'gene', y = 'epoQP', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
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

#------------------------------------------------------------------------

sig_epo <- sig[sig$prot == "epo",]
for (i in 1:nrow(sig_epo)) {
  gene_name = sig_epo$external_gene_name[i]
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
  colnames(cordata) <- c('gene','epoQP')
  ggscatter(cordata, x = 'gene', y = 'epoQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', Log2 TPM', sep=''), ylab = 'EPO production, pg/cell.day')
  ggsave(paste('/Users/rasools/Desktop/corelA/epo/',gene_name,'_epo.png', sep = ""),width = 10, height = 10, units="cm")
}


sig_gfp <- sig[sig$prot == "gfp",]
for (i in 1:nrow(sig_gfp)) {
  gene_name = sig_gfp$external_gene_name[i]
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
  colnames(cordata) <- c('gene','gfpQP')
  ggscatter(cordata, x = 'gene', y = 'gfpQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', Log2 TPM', sep=''), ylab = 'GFP production, Relative to GFPpoly')
  ggsave(paste('/Users/rasools/Desktop/corelA/gfp/',gene_name,'_epo.png', sep = ""), width = 10, height = 10, units="cm")
}


sigER <- sig[which(sig$external_gene_name %in% erGsc[,1] & sig$prot == 'epo'),]

riboSig_gfp <- sig[which(sig$external_gene_name %in% riboGsc[,1] & sig$prot == 'gfp'),]


#------------------------------------------------------------------------

samples_epo <- sam[,c(2,grep('EPO',colnames(sam)))]
samples_epo['epoQP',] <- samples_epo['productivity',]
y_epo <- tpm[,colnames(samples_epo)]
gene_name = "EPOR"
gene = which(tpm$external_gene_name == gene_name)
cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
colnames(cordata) <- c('gene','epoQP')
ggscatter(cordata, x = 'epoQP', y = 'gene', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = 'EPO production, pg/cell.day', ylab = paste(gene_name, ' Expression, TPM', sep=''))
ggsave(paste('/Users/rasools/Desktop/',gene_name,'_epo.svg', sep = ""),width = 10, height = 7, units="cm")


samples_gfp <- sam[,c(1,grep('GFP',colnames(sam)))]
samples_gfp['gfpQP',] <- samples_gfp['productivity',]
y_gfp <- tpm[,colnames(samples_gfp)]
gene_name = "EPOR"
gene = which(tpm$external_gene_name == gene_name)
cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
colnames(cordata) <- c('gene','gfpQP')
ggscatter(cordata, x = 'gfpQP', y = 'gene', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = 'GFP production, Relative to GFPpoly', ylab = paste(gene_name, ' Expression, TPM', sep=''))
ggsave(paste('/Users/rasools/Desktop/',gene_name,'_gfp.svg', sep = ""),width = 10, height = 7, units="cm")

