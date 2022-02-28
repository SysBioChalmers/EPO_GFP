library(DESeq2)
library(EnhancedVolcano)
library(magrittr)
library(biomaRt)

# sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] magrittr_2.0.2              EnhancedVolcano_1.10.0      ggrepel_0.9.1              
# [4] DESeq2_1.32.0               SummarizedExperiment_1.22.0 Biobase_2.52.0             
# [7] MatrixGenerics_1.4.3        matrixStats_0.61.0          GenomicRanges_1.44.0       
# [10] GenomeInfoDb_1.28.4         IRanges_2.26.0              S4Vectors_0.30.2           
# [13] BiocGenerics_0.38.0         piano_2.8.0                 ggpubr_0.4.0               
# [16] viridis_0.6.2               viridisLite_0.4.0           forcats_0.5.1              
# [19] stringr_1.4.0               purrr_0.3.4                 readr_2.1.2                
# [22] tidyr_1.2.0                 tibble_3.1.6                ggplot2_3.3.5              
# [25] tidyverse_1.3.1             reshape_0.8.8               dplyr_1.0.7                
# [28] biomaRt_2.48.3             

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

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
a <- c("EPO7","EPO8","EPOB9","EPOI2","EPOpoly")
for (i in 1:length(a)) {
  res_epo <- results(ds_epo, contrast=c('sampleName','EPOF21',a[i]))
  res_epo <- res_epo[!is.na(res_epo$padj), ]
  res_epo <- as.data.frame(res_epo)
  res_epo <- merge(x=res_epo, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
  name_epo = paste('EPOF21_',a[i], sep="")
  write.table(res_epo, paste('./results/DE_results/', name_epo,".txt", sep=""), 
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep='\t') 
  # svg(paste('/results/DE_results/', name_epo,".svg", sep=""))
  # EnhancedVolcano(res_epo,
  #                 title = name_epo,
  #                 lab = rownames(res_epo),
  #                 x = 'log2FoldChange',
  #                 y = 'pvalue',
  #                 selectLab = c('EPO_HPC4_splitGFP'),
  #                 xlab = bquote(~Log[2]~ 'fold change'),
  #                 pCutoff = 0.05,
  #                 FCcutoff = 1.0,
  #                 pointSize = 1.0,
  #                 labSize = 4.0,
  #                 colAlpha = 1,
  #                 legendPosition = 'right',
  #                 legendLabSize = 14,
  #                 legendIconSize = 5.0)
  # dev.off()
}


