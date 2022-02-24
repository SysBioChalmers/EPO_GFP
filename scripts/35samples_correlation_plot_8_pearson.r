# import needed packages
  library(biomaRt)
  library(dplyr)
  library(reshape)
  library(tidyverse)
  library(viridis)
  library(ggpubr)
  library(piano)

  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
  ensembl2name <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=mart)
  #gsc <- loadGSC('/Users/rasools/Box Sync/phd/ppf/293F_PSN/GSA/gene_sets/c5.all.v7.1.symbols.gmt')
  #riboGsc <- as.data.frame(gsc[["gsc"]][["GO_RIBOSOME"]])

# import tpms
  tpm <- read.delim("/Users/rasools/drive/projects/EPO_GFP/EPO_GFP_github/data/avtpm.txt")

# import sam table
  sam <- read.delim("/Users/rasools/Box Sync/phd/ppf/293F_PSN/sam.txt")
  rownames(sam) <- sam[,1]
  sam <- sam[,2:17]
  samples_epo <- sam[,grep('EPO',colnames(sam))]
  samples_epo['epoQP',] <- samples_epo['productivity',]
  samples_gfp <- sam[,grep('GFP',colnames(sam))]
  samples_gfp['gfpQP',] <- samples_gfp['productivity',]

# prepreation for correlation analysis
  y_epo <- tpm[,colnames(samples_epo)]
  y_gfp <- tpm[,colnames(samples_gfp)]

method = "pearson"
# correlation analysis
  corel <- matrix(nrow=nrow(y_gfp), ncol=4)
  for(i in 1:nrow(y_gfp)) {
    corel[i,1] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["estimate"]]
    corel[i,2] <- cor.test(as.numeric(samples_epo['epoQP',]),as.numeric(y_epo[i,]), method=method)[["p.value"]]
    corel[i,3] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["estimate"]]
    corel[i,4] <- cor.test(as.numeric(samples_gfp['gfpQP',]),as.numeric(y_gfp[i,]), method=method)[["p.value"]]
  }
  corel <- as.data.frame(corel)
  colnames(corel) <- c("epo_coef",'epo_pval','gfp_coef',"gfp_pval")
  rownames(corel) <- rownames(y_gfp)

# finding significantly correlating genes
  df <- corel
  rownames(df) <- rownames(corel)
  th = 0.5
  ## filtering non significant genes based on pvalue lower than 0.05
    df$epo_coef[which(df$epo_pval > 0.05)] = 0
    df$gfp_coef[which(df$gfp_pval > 0.05)] = 0
  ## filtering genes with mean of tpm lower than 10
    df$epo_coef[which(rowMeans(tpm[,3:8]) < 10)] <- 0
    df$gfp_coef[which(rowMeans(tpm[,9:16]) < 10)] <- 0
  ## taging high correlating genes
    df$sig <- 'non'
    df$sig[which(df$epo_coef > th)] <- 'epo+'
    df$sig[which(df$epo_coef < -th)] <- 'epo-'
    df$sig[which(df$gfp_coef > th)] <- 'gfp+'
    df$sig[which(df$gfp_coef < -th)] <- 'gfp-'
    df$sig[which(df$epo_coef > th & df$gfp_coef > th)] <- 'epo+gfp+'
    df$sig[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'epo-gfp-'
    df$sig[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'epo+gfp-'
    df$sig[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'epo-gfp+'
    df$prot <- 'non'
    df$prot[which(df$epo_coef > th)] <- 'epo'
    df$prot[which(df$epo_coef < -th)] <- 'epo'
    df$prot[which(df$gfp_coef > th)] <- 'gfp'
    df$prot[which(df$gfp_coef < -th)] <- 'gfp'
    df$prot[which(df$epo_coef > th & df$gfp_coef > th)] <- 'both'
    df$prot[which(df$epo_coef < -th & df$gfp_coef < -th)] <- 'both'
    df$prot[which(df$epo_coef > th & df$gfp_coef < -th)] <- 'both'
    df$prot[which(df$epo_coef < -th & df$gfp_coef > th)] <- 'both'
  ## filtering non correlating genes
    sig <- df[-which(df$prot == 'non'),]
    sig <- merge(x=sig, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
    rownames(sig) <- sig$Row.names
    sig <- distinct(sig)
    sig <- sig[!is.na(sig$external_gene_name),]
    sig <- sig[!duplicated(sig$external_gene_name),]

# saving tables
  write.table(sig[, c(1)], paste("/Users/rasools/Desktop/35samples_corel_pear_epo_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries
  #write.table(sig[sig$prot == 'epo'|sig$prot == 'both', c(1)], paste("/Users/rasools/Desktop/35samples_corel_pear_epo_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries
  #write.table(sig[sig$prot == 'gfp'|sig$prot == 'both', c(1)], paste("/Users/rasools/Desktop/35samples_corel_pear_gfp_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries
  #write.table(sig[sig$sig == 'epo+', 1], paste("/Users/rasools/Desktop/35samples_corel_pear_epoUp_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries
  #write.table(sig[sig$sig == 'epo-', 1], paste("/Users/rasools/Desktop/35samples_corel_pear_epoDn_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') # 47 entries
  #write.table(sig[sig$sig == 'gfp+', 1], paste("/Users/rasools/Desktop/35samples_corel_pear_gfpUp_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries
  #write.table(sig[sig$sig == 'gfp-', 1], paste("/Users/rasools/Desktop/35samples_corel_pear_gfpDn_thr_",th,".txt",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, , sep='\t') #  entries

# circle plot
  meltedsig <- melt(data=sig, id.vars=c('external_gene_name','sig'), measure.vars=c('epo_coef', 'gfp_coef'))
  meltedsig$value[abs(meltedsig$value) < th] <- 0
  ## plot adjustment
  empty_bar <- 2
  nObsType <- nlevels(as.factor(meltedsig$variable))
  to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(meltedsig$sig))*nObsType, ncol(meltedsig)))
  colnames(to_add) <- colnames(meltedsig)
  to_add$sig <- rep(levels(as.factor(meltedsig$sig)), each=empty_bar*nObsType)
  meltedsig <- rbind(meltedsig, to_add)
  meltedsig <- meltedsig %>% arrange(sig, external_gene_name)
  meltedsig$id <- rep( seq(1, nrow(meltedsig)/nObsType) , each=nObsType)
  ## Get the name and the y position of each label
  label_data <- meltedsig %>% group_by(id, external_gene_name) %>% dplyr::summarize(tot=max(value))
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  ## prepare a data frame for base lines
  base_data <- meltedsig %>% 
    group_by(sig) %>% 
    dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  ## prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1,]
  p <- ggplot(meltedsig) +        
  ## Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=variable), stat="identity", alpha=0.9, position='dodge') +
  scale_fill_viridis(discrete=TRUE) +
  ## Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = -1, xend = start, yend = -1), colour = "Black", alpha=1, size=1 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "Black", alpha=1, size=1 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "Black", alpha=1, size=1 , inherit.aes = FALSE ) +  
  ## Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(meltedsig$id),5), y = c(-1, -0.5, 0, 0.5, 1), label = c("-1", '-0.5', "0", "0.5","1") , color="Black", size=18 , angle=0, fontface="bold", hjust=1) +  
  ylim(-1,1) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +  
  ## Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=1, label=external_gene_name, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=7.5, angle= label_data$angle, inherit.aes = FALSE ) +
  ## Add base line information
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  scale_fill_manual(values = c("red", "seagreen3"))
  ggsave(p, file=paste("/Users/rasools/Desktop/35samples_corel_pear_thr_",th,".png",sep=""), width=45, height=45, limitsize = FALSE)

#------------------------------------------------------------------------
tpm <- merge(x=tpm, y=ensembl2name, by.x=0, by.y=1, all.x=TRUE)
#------
gene_name = "ATF6B"
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
  colnames(cordata) <- c('gene','epoQP')
  ggscatter(cordata, x = 'gene', y = 'epoQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', TPM', sep=''), ylab = 'EPO production, pg/cell.day')
  ggsave(paste('/Users/rasools/Desktop/',gene_name,'_epo.png', sep = ""),width = 10, height = 7, units="cm")


gene_name = "EPOR"
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
  colnames(cordata) <- c('gene','gfpQP')
  ggscatter(cordata, x = 'gene', y = 'gfpQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', TPM', sep=''), ylab = 'GFP production, Relative to GFPpoly')
  ggsave(paste('/Users/rasools/Desktop/',gene_name,'_gfp.png',sep = ""),width=10,height=7, units="cm")

#------------------------------------------------------------------------

sig_epo <- sig[sig$prot == "epo",]
for (i in 1:nrow(sig_epo)) {
  gene_name = sig_epo$external_gene_name[i]
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
  colnames(cordata) <- c('gene','epoQP')
  ggscatter(cordata, x = 'gene', y = 'epoQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', Log2 TPM', sep=''), ylab = 'EPO production, pg/cell.day')
  ggsave(paste('/Users/rasools/Desktop/corelA/epo/',gene_name,'_epo.png', sep = ""),width = 10, height = 10, units="cm")
}


sig_gfp <- sig[sig$prot == "gfp",]
for (i in 1:nrow(sig_gfp)) {
  gene_name = sig_gfp$external_gene_name[i]
  gene = which(tpm$external_gene_name == gene_name)
  cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
  colnames(cordata) <- c('gene','gfpQP')
  ggscatter(cordata, x = 'gene', y = 'gfpQP', 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = paste(gene_name, ', Log2 TPM', sep=''), ylab = 'GFP production, Relative to GFPpoly')
  ggsave(paste('/Users/rasools/Desktop/corelA/gfp/',gene_name,'_epo.png', sep = ""), width = 10, height = 10, units="cm")
}


sigER <- sig[which(sig$external_gene_name %in% erGsc[,1] & sig$prot == 'epo'),]

riboSig_gfp <- sig[which(sig$external_gene_name %in% riboGsc[,1] & sig$prot == 'gfp'),]


#------------------------------------------------------------------------

samples_epo <- sam[,c(2,grep('EPO',colnames(sam)))]
samples_epo['epoQP',] <- samples_epo['productivity',]
y_epo <- tpm[,colnames(samples_epo)]
gene_name = "EPOR"
gene = which(tpm$external_gene_name == gene_name)
cordata <- as.data.frame(t(rbind(y_epo[gene,],samples_epo['epoQP',])))
colnames(cordata) <- c('gene','epoQP')
ggscatter(cordata, x = 'epoQP', y = 'gene', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = 'EPO production, pg/cell.day', ylab = paste(gene_name, ' Expression, TPM', sep=''))
ggsave(paste('/Users/rasools/Desktop/',gene_name,'_epo.svg', sep = ""),width = 10, height = 7, units="cm")


samples_gfp <- sam[,c(1,grep('GFP',colnames(sam)))]
samples_gfp['gfpQP',] <- samples_gfp['productivity',]
y_gfp <- tpm[,colnames(samples_gfp)]
gene_name = "EPOR"
gene = which(tpm$external_gene_name == gene_name)
cordata <- as.data.frame(t(rbind(y_gfp[gene,],samples_gfp['gfpQP',])))
colnames(cordata) <- c('gene','gfpQP')
ggscatter(cordata, x = 'gfpQP', y = 'gene', 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = 'GFP production, Relative to GFPpoly', ylab = paste(gene_name, ' Expression, TPM', sep=''))
ggsave(paste('/Users/rasools/Desktop/',gene_name,'_gfp.svg', sep = ""),width = 10, height = 7, units="cm")







