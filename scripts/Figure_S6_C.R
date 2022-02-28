library(UpSetR)
library(biomaRt)

EPO <- scan('results/DE_analysis/EPO7_vs_EPOF21.txt', what = list(EPO7=''))
EPO <- append(EPO, scan('results/DE_analysis/EPO8_vs_EPOF21.txt', what = list(EPO8='')))
EPO <- append(EPO, scan('results/DE_analysis/EPOB9_vs_EPOF21.txt', what = list(EPOB9='')))
EPO <- append(EPO, scan('results/DE_analysis/EPOI2_vs_EPOF21.txt', what = list(EPOI2='')))

pdf('results/DE_analysis/f21_vs_epo_upset.pdf', width=5, height=4)
upset(fromList(EPO), order.by = "freq", nsets = 4 , queries = list(list(query = intersects, params = list("EPO7", "EPO8", "EPOB9", "EPOI2"), color = "red", active = T)))
dev.off()



#-----------------------------------------------------

a <- as.data.frame(Reduce(intersect, EPO))
a <- a[complete.cases(a), ]

write.table(a$symbol, "results/DE_analysis/common_DE_epo", quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t')
TPM <- read.delim("data/avCounts.txt")
TPM_epo <- TPM[,2:8]
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = rownames(TPM_epo),
                 mart = ensembl )
symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_gene_id, paste, collapse="; ")
TPM_epo$symbol <- symbols[ rownames(TPM_epo) ]

TPM_epo_sig <- TPM_epo[which(TPM_epo$symbol %in% a),]
























