
library(DESeq2)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(biomaRt)

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
          c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))

#---------------------------------------------------------------epo_plots
# #for comparing EPO-producers with HEK293freestyle control cell-line
samples_epo <- samples[samples$protein == "EPO",]
counts_epo <- counts[, rownames(samples_epo)]
dds_epo <- DESeqDataSetFromMatrix(countData=counts_epo, colData=samples_epo, design = ~ sampleName)
ds_epo <- DESeq(dds_epo)
a <- c("EPO7","EPO8","EPOB9","EPOI2","EPOpoly")
i = 1
res_epo <- results(ds_epo, contrast=c('sampleName','EPOF21',a[i]))
res_epo <- res_epo[!is.na(res_epo$padj), ]
res_epo <- as.data.frame(res_epo)
res_epo <- merge(x=res_epo, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
name_epo = paste('EPOF21_',a[i], sep="")

write.table(res_epo, paste('results/', name_epo,".txt", sep=""), 
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t') 


sig_epo <- res_epo[ which(res_epo$padj < 0.05), ]
sig_epo <- sig_epo[ which(abs(sig_epo$log2FoldChange) > 1), ]
sig_epo <- as.data.frame(sig_epo)
up <- dim(sig_epo[sig_epo$log2FoldChange > 1,])[1]
dn <- dim(sig_epo[sig_epo$log2FoldChange < -1,])[1]


name_epo = paste('EPOF21_',a[1], sep="")
svg(paste('/Users/rasools/Desktop/', name_epo,".svg", sep=""))
	EnhancedVolcano(res_epo,
		title = name_epo,
		lab = rownames(res_epo),
		x = 'log2FoldChange',
		y = 'pvalue',
		selectLab = c('EPO_HPC4_splitGFP'),
		xlab = bquote(~Log[2]~ 'fold change'),
		pCutoff = 0.05,
		FCcutoff = 1.0,
		pointSize = 1.0,
		labSize = 4.0,
		colAlpha = 1,
		legendPosition = 'right',
		legendLabSize = 14,
		legendIconSize = 5.0)
dev.off()






