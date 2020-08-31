# Network of correlating genes with EPO and GFP
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
# [1] parallel  stats4    stats     graphics  grDevices utils    
# [7] datasets  methods   base     

# other attached packages:
#  [1] biomaRt_2.44.0              DESeq2_1.28.1              
#  [3] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
#  [5] matrixStats_0.56.0          Biobase_2.48.0             
#  [7] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
#  [9] IRanges_2.22.2              S4Vectors_0.26.1           
# [11] BiocGenerics_0.34.0  
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
pathToEpoGfp <- '/path/to/local/epo_gfp_repo/'

library(DESeq2)
library(biomaRt)
#Load averages of TPM values between replicates
TPM <- read.delim(paste(pathToEpoGfp,'data/avtpm.txt',sep=''))
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
TPM <- merge(x=TPM, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)

#------epo
#subnetwork_epo.sif has generated in networkanalyst.ca from EPO-correlated genes in table S4 and downloaded as .sif
epo_net <- read.delim(paste(pathToEpoGfp,'results/subnetwork_epo.sif',sep=''), header=F)
df_epo <- unique(c(epo_net$V1, epo_net$V3))
tpm_epoNet <- TPM[match(df_epo, TPM$external_gene_name), c(4:8,18)];
tpm_epoNet <- na.omit(tpm_epoNet)
#finding genes with tpm lower than 10
tpm_epoNet_10 <- tpm_epoNet[rowMeans(tpm_epoNet[,1:5]) < 10,]
#removing lines with these genes from primary network
new_epo_net <- epo_net[!(epo_net$V1 %in% tpm_epoNet_10$external_gene_name | epo_net$V3 %in% tpm_epoNet_10$external_gene_name),c(1,3)]
new_epo_net <-  new_epo_net[complete.cases(new_epo_net),]
write.table(new_epo_net, file = paste(pathToEpoGfp,'results/new_epo_net.txt',sep=''), quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
#------gfp
#subnetwork_gfp.sif has generated in networkanalyst.ca from GFP-correlated genes in table S4 and downloaded as .sif
epo_net <- read.delim(paste(pathToEpoGfp,'results/subnetwork_gfp.sif',sep=''), header=F)
df_gfp <- unique(c(gfp_net$V1, gfp_net$V3))
tpm_gfpNet <- TPM[match(df_gfp, TPM$external_gene_name), c(10:16,18)];
tpm_gfpNet <- na.omit(tpm_gfpNet)
#finding genes with tpm lower than 10
tpm_gfpNet_10 <- tpm_gfpNet[rowMeans(tpm_gfpNet[,1:7]) < 10,]
#removing lines with these genes from primary network
new_gfp_net <- gfp_net[!(gfp_net$V1 %in% tpm_gfpNet_10$external_gene_name | gfp_net$V3 %in% tpm_gfpNet_10$external_gene_name),c(1,3)]
new_gfp_net <-  new_gfp_net[complete.cases(new_gfp_net),]
write.table(new_gfp_net, file = paste(pathToEpoGfp,'results/new_gfp_net.txt',sep=''), quote = FALSE, row.names = FALSE, col.names = FALSE, sep='\t')
