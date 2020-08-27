
library(DESeq2)
# import counts
counts <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/alignment/tximport_outputs/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))
#counts$EPOF21_I_6 <- NULL
counts <- counts[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# import TPMs
TPM <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/alignment/tximport_outputs/35samples_TPM.txt")
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
colnames(TPM) <- gsub("abundance.","", colnames(TPM))
#TPM$EPOF21_I_6 <- NULL
TPM <- TPM[c(2,1,4,3,6,5,8,7,10,9,13,11,12,15,14,17,16,19,18,21,20,24,25,22,23,27,26,29,28,31,30,33,32,35,34)]
# check
all(colnames(TPM) == colnames(counts))
# import samples table
kallisto.dir <- list.files(path = "/Users/rasools/Box\ Sync/phd/ppf/293F_PSN/alignment/kallisto_outputs")
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))
rownames(samples) <- samples[,1]
all(rownames(samples) == colnames(counts))

# remove line related to EPO_F_21_6
# samples <- samples[-c(11),]
# add epoTPM column to samples table
samples$epoL2TPM <- log(as.vector(t(TPM["EPO_HPC4_splitGFP",])),2)
samples$epoL2TPM[which(!is.finite(samples$epoL2TPM))] <- 0
# check
#---------------------------------------------------------------
# #for comparing GFP-producers with HEK293F control cell-line
 samples <- samples[samples$cellType == "HEK293F",]
# samples$gfpTPM <- t(counts["pD2529_GFP_extra",])
 counts <- counts[, rownames(samples)]
#---------------------------------------------------------------
# for comparing EPO-producers with f-free control cell-line
 # samples <- samples[samples$cellType == "HEK293freestyle",]
 # counts <- counts[,3:17]
 # samples$contrast = rep(c("cntrl","case","case","case","case","case","case"), c(2,2,2,2,3,2,2))
 # all(rownames(samples) == colnames(counts))

 #########################,CTRL,EPO7,EPO8,EPOB9,EPOF21,EPOI2,EPOpoly

# samples$categ = rep(c("cntrl","low","ERs","low","high","low","low"), c(2,2,2,2,2,2,2))
#---------------------------------------------------------------
# for comparing EPO-producers
# samples <- samples[samples$protein == "EPO",]
# counts <- counts[,5:17]
# samples$contrast = rep(c("cntrl","out","cntrl","case","cntrl","out"), c(2,2,2,3,2,2))
##########################,EPO7,EPO8,EPOB9,EPOF21,EPOI2,EPOpoly
#---------------------------------------------------------------
# for comparing EPO-control with GFP-control samples
# samples <- samples[samples$producibility == "Control",]
# counts <- counts[,1:4]
#---------------------------------------------------------------
# #for comparing EPO-producers with GFP-producers
# samples <- samples[samples$producibility == "Producer",]
# counts <- counts[,5:34]
#---------------------------------------------------------------
# #for comparing EPO-producers
# samples <- samples[5:8,]
# counts <- counts[,5:8]
#---------------------------------------------------------------
# #for comparing EPO8 and EPO F21
# samples <- samples[7:12,]
# samples <- samples[-c(3:4),]

# counts <- counts[,7:12]
# counts <- counts[,-c(3:4)]
#---------------------------------------------------------------
# #for comparing EPO7 and EPO F21
# samples <- samples[5:12,]
# samples <- samples[-c(3:6),]

# counts <- counts[,5:12]
# counts <- counts[,-c(3:6)]
#---------------------------------------------------------------
# #for comparing EPO8 and EPO7
# samples <- samples[5:8,]
# counts <- counts[,5:8]
#---------------------------------------------------------------
# #for comparing all EPO producers against EPOF21
# samples <- samples[samples$protein == "EPO",]
# counts <- counts[,5:16]

# samples$protein[7] <- "GFP"
# samples$protein[8] <- "GFP"
#---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ producibility)
ds <- DESeq(dds)
#---------------------------------------------------------------
resultsNames(ds)
#res <- results(ds, name="epoL2TPM")
res <- results(ds, contrast=c("producibility", "Producer", "Control"))
res <- res[ ! is.na(res$padj), ]
#---------------------------------------------------------------
sig <- res
sig <- res[ which(res$padj < 0.05), ]
sig <- sig[ which(abs(sig$log2FoldChange) > 1), ]
sig <- sig[ order(sig$log2FoldChange), ]
sig <- as.data.frame(sig)
a <- as.data.frame(row.names(sig))
rownames(a) <- a$`row.names(sig)`
#---------------------------------------------------------------
library(biomaRt)
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = rownames(sig),
                 mart = ensembl )
symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_gene_id, paste, collapse="; ")
sig$symbol <- symbols[ rownames(sig) ]
#---------------------------------------------------------------
write.table(sig, "/Users/rasools/Desktop/f21_pairwise/F21_vs_EPOproducers_excl_EPOpoly&EPO8", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# write.table(sig$symbol, "/Users/rasools/Desktop/all_vs_F21_sym", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(sig, "/Users/rasools/Desktop/F21_vs_producers", sep="\t", quote = FALSE)
#---------------------------------------------------------------
# write.table(a, "/Users/rasools/Box Sync/phd/ppf/293F_PSN/DEA/geneLists_ID/proVScntrl_ID", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(sig$symbol, "/Users/rasools/Box Sync/phd/ppf/293F_PSN/DEA/geneLists/proVScntrl_sym", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(sig, "/Users/rasools/Box Sync/phd/ppf/293F_PSN/DEA/fullResults/proVScntrl",col.names = NA)
#---------------------------------------------------------------
d <- plotCounts(ds, gene="EPO_HPC4_splitGFP", intgroup="protein", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=protein, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))  
  #+ scale_y_log10(breaks=c(25,100,400))
#---------------------------------------------------------------ReportingTools
library("ReportingTools")
desReport <- HTMLReport(shortName = 'EPO_producers',
    title = 'EPO_producers',
    reportDirectory = "./reports")
publish(ds,desReport, pvalueCutoff=0.1,
	annotation.db="org.Mm.eg.db", factor = samples$sampleName,
	reportDir="./reports")
finish(desReport)
#---------------------------------------------------------------regionReport
library("regionReport")
report <- DESeq2Report(ds, project = "EPO_producers", c('sampleName','protein'), outdir = "DESeq2Report-example")
#---------------------------------------------------------------Heatmap of the count matrix
library("pheatmap")
ntd <- normTransform(ds)
select <- order(rowMeans(counts(ds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(ds)[,c("sampleName","protein")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
#---------------------------------------------------------------
vsd <- vst(ds, blind=FALSE)
df <- as.data.frame(colData(ds)[,c("sampleName","protein")])
select <- order(rowMeans(counts(ds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
#---------------------------------------------------------------
rld <- rlog(ds, blind=FALSE)
df <- as.data.frame(colData(ds)[,c("sampleName","protein")])
select <- order(rowMeans(counts(ds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
#---------------------------------------------------------------Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sampleName, vsd$protein, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#---------------------------------------------------------------Principal component plot of the samples
library("ggplot2")
vsd <- vst(ds, blind=FALSE)
rld <- rlog(ds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("sampleName", "protein"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sampleName, shape=protein)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#---------------------------------------------------------------
dds$group <- dds$growth*dds$productivity
design(dds) <- ~ group
ds <- DESeq(dds)
resultsNames(ds)
res <- results(ds, contrast=c("group", 0.081000, 0.333600))









