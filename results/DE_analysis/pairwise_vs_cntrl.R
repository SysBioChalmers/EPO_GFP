
library(DESeq2)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)

# import counts
counts <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/alignment/tximport_outputs/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
#counts$EPOF21_I_6 <- NULL
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# import samples table
kallisto.dir <- list.files(path = "/Users/rasools/Box\ Sync/phd/ppf/293F_PSN/alignment/kallisto_outputs")
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))

#---------------------------------------------------------------epo_plots
# #for comparing EPO-producers with HEK293freestyle control cell-line
samples_epo <- samples[samples$cellType == "HEK293freestyle",]
counts_epo <- counts[, rownames(samples_epo)]
dds_epo <- DESeqDataSetFromMatrix(countData=counts_epo, colData=samples_epo, design = ~ sampleName)
ds_epo <- DESeq(dds_epo)
a <- resultsNames(ds_epo)
for (i in 2:length(a)) {
	res_epo <- results(ds_epo, name=a[i])
	res_epo <- res_epo[ ! is.na(res_epo$padj), ]
	sig_epo <- res_epo[ which(res_epo$padj < 0.05), ]
	sig_epo <- sig_epo[ which(abs(sig_epo$log2FoldChange) > 1), ]
	sig_epo <- as.data.frame(sig_epo)
	sig_epo <- merge(x=sig_epo, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
	name_epo <- a[i]; name_epo <- gsub("sampleName_","",a[i]);
	write.table(sig_epo$external_gene_name, paste('/Users/rasools/Desktop/', name_epo,".txt", sep=""), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
#---------------------------------------------------------------gfp_plots
# #for comparing GFP-producers with HEK293F control cell-line
samples_gfp <- samples[samples$cellType == "HEK293F",]
counts_gfp <- counts[, rownames(samples_gfp)]
dds_gfp <- DESeqDataSetFromMatrix(countData=counts_gfp, colData=samples_gfp, design = ~ sampleName)
ds_gfp <- DESeq(dds_gfp)
b <- resultsNames(ds_gfp)
for (i in 2:length(b)) {
	res_gfp <- results(ds_gfp, name=b[i])
	res_gfp <- res_gfp[ ! is.na(res_gfp$padj), ]
	sig_gfp <- res_gfp[ which(res_gfp$padj < 0.05), ]
	sig_gfp <- sig_gfp[ which(abs(sig_gfp$log2FoldChange) > 1), ]
	sig_gfp <- as.data.frame(sig_gfp)
	sig_gfp <- merge(x=sig_gfp, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
	name_gfp <- b[i]; name_gfp <- gsub("sampleName_","",b[i]);
	write.table(sig_gfp$external_gene_name, paste('/Users/rasools/Desktop/', name_gfp,".txt", sep=""), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}





