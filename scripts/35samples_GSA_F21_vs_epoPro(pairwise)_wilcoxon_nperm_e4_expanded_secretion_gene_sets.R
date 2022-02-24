
library(DESeq2)
# import counts
counts <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/alignment/tximport_outputs/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# import TPMs
TPM <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/alignment/tximport_outputs/35samples_TPM.txt")
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
colnames(TPM) <- gsub("abundance.","", colnames(TPM))
TPM <- TPM[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# check
all(colnames(TPM) == colnames(counts))
# import samples table
kallisto.dir <- list.files(path = "/Users/rasools/Box\ Sync/phd/ppf/293F_PSN/alignment/kallisto_outputs")
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))
samples$epoL2TPM <- log(as.vector(t(TPM["EPO_HPC4_splitGFP",])),2)
samples$epoL2TPM[which(!is.finite(samples$epoL2TPM))] <- 0
#---------------------------------------------------------------
# for comparing EPO-producers with f-free control cell-line
samples <- samples[samples$protein == "EPO",]
counts <- counts[,5:17]
all(rownames(samples) == colnames(counts))
#---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ sampleName)
ds <- DESeq(dds)
#---------------------------------------------------------------
library(piano)
library(biomaRt)
library(snowfall)
library(snow)
sec <- read.delim('/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/gene_sets/secretion_gene_sets.txt')
gsc <- loadGSC('/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/gene_sets/c5.all.v7.0.symbols.gmt')
gs <- gsc$gsc
genes=as.character(sec$Gene)
res <- runGSAhyper(genes=as.character(sec$Gene), gsc=gsc)
pval <- as.data.frame(res$pvalues)
pval$gs <- rownames(pval)
colnames(pval) <- c("pval","name")
siggs <- pval[pval$pval<0.05,]
sec_gsc <- gsc$gsc[which(names(gsc$gsc) %in% siggs$name)]
sec_gsc <- list(sec_gsc,gsc$addInfo[which(names(gsc$gsc) %in% siggs$name),])
names(sec_gsc) <- c("gsc","addInfo")
class(sec_gsc) <- "GSC"
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
#--------------------------------------------------------------- EPOf21 vs. EPO8
res_epo8 <- as.data.frame(results(ds, contrast=c("sampleName","EPOF21","EPO8")))
res_epo8 <- res_epo8[ ! is.na(res_epo8$padj), ]

geneLevelStats_epo8 <- as.data.frame(res_epo8[,c("log2FoldChange","padj")])
# Merge it with our gene-level statistics:
geneLevelStats_epo8 <- merge(x=geneLevelStats_epo8, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_epo8) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_epo8) <- geneLevelStats_epo8[,1]
geneLevelStats_epo8 <- geneLevelStats_epo8[order(geneLevelStats_epo8[,'gene'],geneLevelStats_epo8[,'padj'],decreasing = c(FALSE,FALSE)),]
geneLevelStats_epo8 <- geneLevelStats_epo8[!duplicated(geneLevelStats_epo8$gene),]
padj_epo8 <- geneLevelStats_epo8$padj
log2fc_epo8 <- geneLevelStats_epo8$log2fc
names(padj_epo8) <- names(log2fc_epo8) <- geneLevelStats_epo8$gene
F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets <- runGSA(geneLevelStats=padj_epo8, geneSetStat="wilcoxon", directions=log2fc_epo8, signifMethod="geneSampling", gsc=sec_gsc, ncpus=8, adjMethod="BH", nPerm=1e4)
save(F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets, file = "/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/results/F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData")
#--------------------------------------------------------------- EPOF21 vs. EPO7
res_epo7 <- as.data.frame(results(ds, contrast=c("sampleName","EPOF21","EPO7")))
res_epo7 <- res_epo7[ ! is.na(res_epo7$padj), ]

geneLevelStats_epo7 <- as.data.frame(res_epo7[,c("log2FoldChange","padj")])
# Merge it with our gene-level statistics:
geneLevelStats_epo7 <- merge(x=geneLevelStats_epo7, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_epo7) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_epo7) <- geneLevelStats_epo7[,1]
geneLevelStats_epo7 <- geneLevelStats_epo7[order(geneLevelStats_epo7[,'gene'],geneLevelStats_epo7[,'padj'],decreasing = c(FALSE,FALSE)),]
geneLevelStats_epo7 <- geneLevelStats_epo7[!duplicated(geneLevelStats_epo7$gene),]
padj_epo7 <- geneLevelStats_epo7$padj
log2fc_epo7 <- geneLevelStats_epo7$log2fc
names(padj_epo7) <- names(log2fc_epo7) <- geneLevelStats_epo7$gene
F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets <- runGSA(geneLevelStats=padj_epo7, geneSetStat="wilcoxon", directions=log2fc_epo7, signifMethod="geneSampling", gsc=sec_gsc, ncpus=8, adjMethod="BH", nPerm=1e4)
save(F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets, file = "/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/results/F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData")
#--------------------------------------------------------------- EPOF21 vs. EPOB9
res_epob9 <- as.data.frame(results(ds, contrast=c("sampleName","EPOF21","EPOB9")))
res_epob9 <- res_epob9[ ! is.na(res_epob9$padj), ]

geneLevelStats_epob9 <- as.data.frame(res_epob9[,c("log2FoldChange","padj")])
# Merge it with our gene-level statistics:
geneLevelStats_epob9 <- merge(x=geneLevelStats_epob9, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_epob9) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_epob9) <- geneLevelStats_epob9[,1]
geneLevelStats_epob9 <- geneLevelStats_epob9[order(geneLevelStats_epob9[,'gene'],geneLevelStats_epob9[,'padj'],decreasing = c(FALSE,FALSE)),]
geneLevelStats_epob9 <- geneLevelStats_epob9[!duplicated(geneLevelStats_epob9$gene),]
padj_epob9 <- geneLevelStats_epob9$padj
log2fc_epob9 <- geneLevelStats_epob9$log2fc
names(padj_epob9) <- names(log2fc_epob9) <- geneLevelStats_epob9$gene
F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets <- runGSA(geneLevelStats=padj_epob9, geneSetStat="wilcoxon", directions=log2fc_epob9, signifMethod="geneSampling", gsc=sec_gsc, ncpus=8, adjMethod="BH", nPerm=1e4)
save(F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets, file = "/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/results/F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData")
#--------------------------------------------------------------- EPOF21 vs. EPOI2
res_epoi2 <- as.data.frame(results(ds, contrast=c("sampleName","EPOF21","EPOI2")))
res_epoi2 <- res_epoi2[ ! is.na(res_epoi2$padj), ]

geneLevelStats_epoi2 <- as.data.frame(res_epoi2[,c("log2FoldChange","padj")])
# Merge it with our gene-level statistics:
geneLevelStats_epoi2 <- merge(x=geneLevelStats_epoi2, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_epoi2) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_epoi2) <- geneLevelStats_epoi2[,1]
geneLevelStats_epoi2 <- geneLevelStats_epoi2[order(geneLevelStats_epoi2[,'gene'],geneLevelStats_epoi2[,'padj'],decreasing = c(FALSE,FALSE)),]
geneLevelStats_epoi2 <- geneLevelStats_epoi2[!duplicated(geneLevelStats_epoi2$gene),]
padj_epoi2 <- geneLevelStats_epoi2$padj
log2fc_epoi2 <- geneLevelStats_epoi2$log2fc
names(padj_epoi2) <- names(log2fc_epoi2) <- geneLevelStats_epoi2$gene
F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets <- runGSA(geneLevelStats=padj_epoi2, geneSetStat="wilcoxon", directions=log2fc_epoi2, signifMethod="geneSampling", gsc=sec_gsc, ncpus=8, adjMethod="BH", nPerm=1e4)
save(F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets, file = "/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/results/F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData")
#--------------------------------------------------------------- EPOF21 vs. EPOpoly
res_epopoly <- as.data.frame(results(ds, contrast=c("sampleName","EPOF21","EPOpoly")))
res_epopoly <- res_epopoly[ ! is.na(res_epopoly$padj), ]

geneLevelStats_epopoly <- as.data.frame(res_epopoly[,c("log2FoldChange","padj")])
# Merge it with our gene-level statistics:
geneLevelStats_epopoly <- merge(x=geneLevelStats_epopoly, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
colnames(geneLevelStats_epopoly) <- c("ensembl","log2fc","padj","gene")
rownames(geneLevelStats_epopoly) <- geneLevelStats_epopoly[,1]
geneLevelStats_epopoly <- geneLevelStats_epopoly[order(geneLevelStats_epopoly[,'gene'],geneLevelStats_epopoly[,'padj'],decreasing = c(FALSE,FALSE)),]
geneLevelStats_epopoly <- geneLevelStats_epopoly[!duplicated(geneLevelStats_epopoly$gene),]
padj_epopoly <- geneLevelStats_epopoly$padj
log2fc_epopoly <- geneLevelStats_epopoly$log2fc
names(padj_epopoly) <- names(log2fc_epopoly) <- geneLevelStats_epopoly$gene
F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets <- runGSA(geneLevelStats=padj_epopoly, geneSetStat="wilcoxon", directions=log2fc_epopoly, signifMethod="geneSampling", gsc=sec_gsc, ncpus=8, adjMethod="BH", nPerm=1e4)
save(F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets, file = "/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/results/F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData")

