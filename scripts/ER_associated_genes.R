
# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***

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
# load(paste(pathToEpoGfp,'results/F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
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
GO_ER_stress_genes <- c(gsc[["gsc"]][["GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]], 
                        gsc[["gsc"]][["GO_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN"]],
                        gsc[["gsc"]][["GO_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]])
#                       gsc[["gsc"]][["GO_SECRETION"]],
#                        gsc[["gsc"]][["GO_PROTEIN_CATABOLIC_PROCESS"]])


mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

# geneSets <- read.delim('results/heatmap_mix_f21_goSlimSecretion_2.txt')
# exploreGSAres(F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets,
#               browser = TRUE, geneAnnot = NULL)
# 
# F21_epo8_ERgenes <- read.delim('results/gene_table_F21_epo8_ER.txt')
# f21_8 <- read.delim('results/EPOF21_EPO8.txt')
# ER_genes_f21_vs_epo8 <- f21_8[match(F21_epo8_ERgenes$Gene.ID, f21_8$external_gene_name),]
# ER_genes_f21_vs_epo8 <- ER_genes_f21_vs_epo8[complete.cases(ER_genes_f21_vs_epo8), ]
# write.table(ER_genes_f21_vs_epo8, 'results/ER_genes_f21_vs_epo8.txt', sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# 
# 
# ER_genes_f21_vs_epo8_sig <- ER_genes_f21_vs_epo8[ER_genes_f21_vs_epo8$padj < 0.05, ]
# ER_genes_f21_vs_epo8_sig <- ER_genes_f21_vs_epo8_sig[abs(ER_genes_f21_vs_epo8_sig$log2FoldChange) > 0.58, ]
# ER_genes_f21_vs_epo8_sig_pos <- ER_genes_f21_vs_epo8_sig[ER_genes_f21_vs_epo8_sig$log2FoldChange > 0.58, c('Row.names','log2FoldChange','padj','external_gene_name')]
# write.table(ER_genes_f21_vs_epo8_sig_pos, 'results/ER_genes_f21_vs_epo8_sig_pos.txt', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# write.table(ER_genes_f21_vs_epo8_sig, 'results/ER_genes_f21_vs_epo8_sig.txt', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# #----------new
# library(DESeq2)
# library(biomaRt)
# library(ggplot2)
# library(dplyr)
# library(reshape)

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

# #---------------------------------------------------------------
# # import TPMs
# TPM <- read.delim("data/35samples_TPM.txt")
# rownames(TPM) <- TPM$Geneid
# TPM <- TPM[,-1]
# colnames(TPM) <- gsub("abundance.","", colnames(TPM))
# TPM <- TPM[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# samples$epoL2TPM <- log(as.vector(t(TPM["EPO_HPC4_splitGFP",])),2)
# samples$epoL2TPM[which(!is.finite(samples$epoL2TPM))] <- 0
# # check
# all(colnames(TPM) == colnames(counts))
#---------------------------------------------------------------
samples <- samples[samples$protein == "EPO",]
counts <- counts[,rownames(samples)]
# TPM <- TPM[,rownames(samples)]
samples$contrast = rep(c("cntrl","cntrl","cntrl","case","cntrl","cntrl"), c(2,2,2,3,2,2))
all(rownames(samples) == colnames(counts))
# all(colnames(TPM) == colnames(counts))
#########################EPO7,EPO8,EPOB9,EPOF21,EPOI2,EPOpoly
#---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ contrast)
ds <- DESeq(dds)
#---------------------------------------------------------------
resultsNames(ds)
res <- results(ds, contrast=c("contrast", "case", "cntrl"))
#---------------------------------------------------------------
sig <- res[ ! is.na(res$padj), ]
sig <-  as.data.frame(sig)
sig <- as.data.frame(sig[sig$padj < 0.05,])
#sig <- as.data.frame(sig[abs(sig$log2FoldChange) > 1,])
sig <- merge(x=sig, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
sig_ER_stress <- sig[match(GO_ER_stress_genes, sig$external_gene_name), ]
sig_ER_stress <- sig_ER_stress[complete.cases(sig_ER_stress), ]
sig_ER_stress <- sig_ER_stress[order(sig_ER_stress$log2FoldChange),]
sig_ER_stress <- distinct(sig_ER_stress)

ERresponse <- gsc[["gsc"]][["GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]
incorrect_prot <- gsc[["gsc"]][["GO_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN"]]
negERstress <- gsc[["gsc"]][["GO_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]
secretoin <- gsc[["gsc"]][["GO_SECRETION"]]
protcat <- gsc[["gsc"]][["GO_PROTEIN_CATABOLIC_PROCESS"]]


sig_ER_stress$GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS <- 'No'
sig_ER_stress$GO_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS[sig_ER_stress$external_gene_name %in% ERresponse] <- 'Yes'

sig_ER_stress$GO_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN <- 'No'
sig_ER_stress$GO_CELLULAR_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN[sig_ER_stress$external_gene_name %in% incorrect_prot] <- 'Yes'

sig_ER_stress$GO_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS <- 'No'
sig_ER_stress$GO_NEGATIVE_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS[sig_ER_stress$external_gene_name %in% negERstress] <- 'Yes'

# sig_ER_stress$GO_SECRETION <- 'No'
# sig_ER_stress$GO_SECRETION[sig_ER_stress$external_gene_name %in% secretoin] <- 'Yes'
# 
# sig_ER_stress$GO_PROTEIN_CATABOLIC_PROCESS <- 'No'
# sig_ER_stress$GO_PROTEIN_CATABOLIC_PROCESS[sig_ER_stress$external_gene_name %in% protcat] <- 'Yes'

write.table(sig_ER_stress, 'results/sig_ER_stress_0.05.txt', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# pos <- sig[sig$log2FoldChange > 0,]
# write.table(pos, 'results/f21_vs_epo_sig_pos_0.05.txt', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# 
# neg <- sig[sig$log2FoldChange < 0,]
# write.table(neg, 'results/f21_vs_epo_sig_neg_0.05.txt', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#--------------------------------------------------------------- heatmap of sig genes

norcounts <- counts(ds, normalized=TRUE)
norcounts_mean <- matrix(nrow=length(norcounts), ncol=5)
colnames(norcounts_mean) <- rep(c("EPOF21","EPO8","EPO7","EPOB9","EPOI2"),1)
norcounts_mean[,1] <- rowMeans(norcounts[,7:9])
norcounts_mean[,2] <- rowMeans(norcounts[,3:4])
norcounts_mean[,3] <- rowMeans(norcounts[,1:2])
norcounts_mean[,4] <- rowMeans(norcounts[,5:6])
norcounts_mean[,5] <- rowMeans(norcounts[,10:11])
#norcounts_mean <- log2(norcounts_mean+1)
norcounts_mean <- as.data.frame(norcounts_mean)
norcounts_mean$id <- rownames(norcounts)
norcounts_mean_sig_ER_stress <- norcounts_mean[match(sig_ER_stress$Row.names, norcounts_mean$id),]
rownames(norcounts_mean_sig_ER_stress) <- norcounts_mean_sig_ER_stress$id
norcounts_mean_sig_ER_stress$id <- c()
norcounts_mean_sig_ER_stress <- merge(x=norcounts_mean_sig_ER_stress, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
rownames(norcounts_mean_sig_ER_stress) <- norcounts_mean_sig_ER_stress$external_gene_name
norcounts_mean_sig_ER_stress$Row.names <- c()
norcounts_mean_sig_ER_stress$external_gene_name <- c()
write.table(norcounts_mean_sig_ER_stress, 'results/sig_ER_stress_0.05_norCounts.txt', sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)























