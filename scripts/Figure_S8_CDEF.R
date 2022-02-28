
library(DESeq2)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(biomaRt)
library(reshape)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)

# import counts
counts <- read.delim("data/35samples_counts.txt")
rownames(counts) <- counts$Geneid
counts <- counts[,-1]
counts <- ceiling(counts)
colnames(counts) <- gsub("counts.","", colnames(counts))

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
a <- c("EPO7","EPO8","EPOB9","EPOI2")
EPOR_genes = read.delim("heatmap_EPOR_related_genes/string_protein_annotations.tsv")
l2fc <- matrix(, nrow = 9, ncol = 4)
padj <- matrix(, nrow = 9, ncol = 4)
rownames(l2fc) = rownames(padj) = EPOR_genes$genes
colnames(l2fc) = colnames(padj) = a

for (i in 1:length(a)) {
	res_epo <- results(ds_epo, contrast=c("sampleName","EPOF21",a[i]))
	res_epo <- res_epo[ ! is.na(res_epo$padj), ]
	#sig_epo <- res_epo[ which(res_epo$padj < 0.05), ]
	#sig_epo <- sig_epo[ which(abs(sig_epo$log2FoldChange) > 1), ]
	sig_epo <- as.data.frame(res_epo)
	sig_epo <- merge(x=sig_epo, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
	#name_epo = paste(a[i],'_vs_EPOF21',sep="")
	l2fc[,i] = sig_epo[match(EPOR_genes$genes,sig_epo$external_gene_name),"log2FoldChange"]
	padj[,i] = sig_epo[match(EPOR_genes$genes,sig_epo$external_gene_name),"padj"]
}

l2fc <- melt(l2fc)
padj <- melt(padj)
df <- l2fc
df$padj <- -log10(padj$value)

ggplot(df ,aes(x = X2, y = X1, size = padj)) + 
	geom_point(aes(color = value)) +
	scale_color_gradient2(midpoint = 0, low = "blue", mid = "gray",
                            high = "red", space = "Lab" ) +
	theme(
	legend.title = element_text(size = 6),
	legend.text = element_text(size = 6),
	panel.background = element_blank(),
	panel.border=element_rect(fill=NA),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	strip.background=element_blank(),
	axis.text.x=element_text(colour="black",size = 7,family="Arial",angle = 30, hjust = 1),
	axis.text.y=element_text(colour="black",size = 7,family="Arial"),
	axis.title=element_text(colour="black",size = 8),
	axis.ticks=element_line(colour="black"),
	plot.margin=unit(c(1,1,1,1),"line")
	) +
ggsave(('results/box_plot_EPOR_TPM.svg'), width=6.5, height=6.5, units="cm")





