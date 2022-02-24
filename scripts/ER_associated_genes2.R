pathToEpoGfp <- '/Users/rasools/drive/projects/EPO_GFP/EPO_GFP_github/'

library("ggplot2")
library("reshape2")
library("dplyr")
library("ggrepel")
library("tidyverse")
library("piano")
library("scales")
library("svglite")
library("cowplot")
library("ggpubr")
library("DESeq2")
# 
load(paste(pathToEpoGfp,'results/F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
# epo8 <- GSAsummaryTable(F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets)
# load(paste(pathToEpoGfp,'results/F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
# epo7 <- GSAsummaryTable(F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets)
# load(paste(pathToEpoGfp,'results/F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
# epob9 <- GSAsummaryTable(F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets)
# load(paste(pathToEpoGfp,'results/F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
# epoi2 <- GSAsummaryTable(F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets)
# load(paste(pathToEpoGfp,'results/F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
# epopoly <- GSAsummaryTable(F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets)

gsc <- loadGSC('/Users/rasools/drive/projects/EPO_GFP/Files/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.all.v7.0.symbols.gmt')
GO_ER_part_genes <- gsc[["gsc"]][["GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

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

#DEA
samples_epo <- samples[samples$protein == "EPO",]
counts_epo <- counts[, rownames(samples_epo)]
dds_epo <- DESeqDataSetFromMatrix(countData=counts_epo, colData=samples_epo, design = ~ sampleName)
ds_epo <- DESeq(dds_epo)
a <- c("EPO7","EPO8","EPOB9","EPOI2","EPOpoly")

for (i in 1:4){
  res_epo <- results(ds_epo, contrast=c('sampleName','EPOF21',a[i]))
  res_epo <- res_epo[!is.na(res_epo$padj), ]
  res_epo <- as.data.frame(res_epo)
  res_epo <- res_epo[ which(res_epo$padj < 0.05), ]
  res_epo <- res_epo[ which(abs(res_epo$log2FoldChange) > 1), ]
  res_epo <- merge(x=res_epo, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
  name_epo = paste('EPOF21_',a[i], sep="")
  sig_ER_part <- res_epo[match(GO_ER_part_genes, res_epo$external_gene_name), ]
  sig_ER_part <- sig_ER_part[complete.cases(sig_ER_part), ]
  sig_ER_part$comparison = name_epo
  sig_ER_part <- sig_ER_part[order(sig_ER_part$log2FoldChange),]
  # if (i==1) {
  #   sig_genes = sig_ER_part[which(sig_ER_part$padj < 0.05), c(1,8,9)]
  # } else {
  #   sig_genes = rbind(sig_genes, sig_ER_part[ which(sig_ER_part$padj < 0.05), c(1,8,9)])
  # }
  if (i==1) {
     top_genes = sig_ER_part[1:5, ]
     top_genes = rbind(top_genes, sig_ER_part[dim(sig_ER_part)[1]-4:dim(sig_ER_part)[1], ])
     
     
  rm(res_epo, name_epo, sig_ER_part)
  
}

sig_genes <- distinct(sig_genes[, c(1,2)])


write.table(sig_genes[,1], paste('results/', 'sum_epo_ER_sig_genes_p.05_l1',".txt", sep=""), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t') 


#------------------------------------------------------ mix_fig table
mix_fig <- matrix(nrow=length(sig_genes), ncol=4)
colnames(mix_fig) <- rep(c("EPO8","EPO7","EPOB9","EPOI2"),1)
rownames(mix_fig) <- unique(top_path$Name)

top_genes <- sig_epo8[1:5,]






















