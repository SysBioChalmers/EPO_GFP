library(UpSetR)
library(biomaRt)

EPO <- scan('results/DE_analysis/EPO7_vs_CTRL293Free.txt', what = list(EPO7=''))
EPO <- append(EPO, scan('results/DE_analysis/EPO8_vs_CTRL293Free.txt', what = list(EPO8='')))
EPO <- append(EPO, scan('results/DE_analysis/EPOB9_vs_CTRL293Free.txt', what = list(EPOB9='')))
EPO <- append(EPO, scan('results/DE_analysis/EPOF21_vs_CTRL293Free.txt', what = list(EPOF21='')))
EPO <- append(EPO, scan('results/DE_analysis/EPOI2_vs_CTRL293Free.txt', what = list(EPOI2='')))

pdf('results/epo_cnt_upset.pdf', width=6.5, height=4)
upset(fromList(EPO), order.by = "freq", nsets = 5 ,queries = list(list(query = intersects, params = list("EPO7", "EPO8", "EPOB9", "EPOF21", "EPOI2"), color = "red", active = T)))
dev.off()

#--------------------------------------------------

GFP <- scan('results/DE_analysis/GFP1_vs_CTRL293F.txt', what = list(GFP1=''))
GFP <- append(GFP, scan('results/DE_analysis/GFP3_vs_CTRL293F.txt', what = list(GFP3='')))
GFP <- append(GFP, scan('results/DE_analysis/GFP25_vs_CTRL293F.txt', what = list(GFP25='')))
GFP <- append(GFP, scan('results/DE_analysis/GFP26_vs_CTRL293F.txt', what = list(GFP26='')))
GFP <- append(GFP, scan('results/DE_analysis/GFP27_vs_CTRL293F.txt', what = list(GFP27='')))
GFP <- append(GFP, scan('results/DE_analysis/GFP28_vs_CTRL293F.txt', what = list(GFP28='')))
GFP <- append(GFP, scan('results/DE_analysis/GFP29_vs_CTRL293F.txt', what = list(GFP29='')))

pdf('results/DE_analysis/gfp_cnt_upset.pdf', width=5.5, height=3.5)
upset(fromList(GFP), order.by = "freq", nsets = 7 ,queries = list(list(query = intersects, params = list("GFP1", "GFP3", "GFP25", "GFP26", "GFP27", "GFP28", "GFP29"), color = "green4", active = T)))
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
























