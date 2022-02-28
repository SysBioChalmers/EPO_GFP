
# import needed packages
  library(biomaRt)
  library(dplyr)
  library(reshape)
  library(tidyverse)
  library(viridis)
  library(ggpubr)
  library(piano)

# import tpms
  tpm <- read.delim("data/avtpm.txt")

samples <- as.data.frame(colnames(tpm))
samples[,2] <- rep(c("GFP","EPO","GFP"), c(1,7,8))
samples[,3] <- t(tpm['ENSG00000187266',])
colnames(samples) <- c("sampleName","EPOR","TPM")
df <- melt(data=samples, id.vars=c('EPOR','sampleName'), measure.vars=c('TPM'))
colnames(df) <- c("EPOR","sampleName","variable","TPM")

ggplot(df, aes(x=EPOR, y=TPM, label=sampleName)) +
  geom_boxplot(aes(), outlier.shape=NA) +
  geom_point() +
  geom_text(size=2.5, hjust=-0.2, vjust=0) +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  panel.background = element_blank(),
  panel.border=element_rect(fill=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_blank(),
  axis.text.x=element_text(colour="black",size = 11,family="Arial",angle = 30, hjust = 1),
  axis.text.y=element_text(colour="black",size = 11,family="Arial"),
  axis.title=element_text(colour="black",size = 11),
  axis.ticks=element_line(colour="black"),
  plot.margin=unit(c(1,1,1,1),"line")
  ) +
ggsave(('results/box_plot_EPOR_TPM.svg'), width=9, height=15, units="cm")



# samples <- as.data.frame(colnames(tpm[,3:16]))
# samples[,2] <- rep(c("EPO","GFP"), c(6,8))
# a = t(tpm['ENSG00000187266',1:2])
# samples[,3] <- t(tpm['ENSG00000187266',3:16])
# samples[1:6,3] <- samples[1:6,3]/a[2]
# samples[7:14,3] <- samples[7:14,3]/a[1]

# colnames(samples) <- c("sampleName","EPOR","TPM")
# df <- melt(data=samples, id.vars=c('EPOR','sampleName'), measure.vars=c('TPM'))
# colnames(df) <- c("EPOR","sampleName","variable","L2FC")

# ggplot(df, aes(x=EPOR, y=L2FC, label=sampleName)) +
#   geom_boxplot(aes(), outlier.shape=NA) +
#   geom_point() +
#   geom_text(size=2.5, hjust=-0.2, vjust=0) +
#   theme(
#   legend.title = element_text(size = 8),
#   legend.text = element_text(size = 8),
#   panel.background = element_blank(),
#   panel.border=element_rect(fill=NA),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(),
#   strip.background=element_blank(),
#   axis.text.x=element_text(colour="black",size = 11,family="Arial",angle = 30, hjust = 1),
#   axis.text.y=element_text(colour="black",size = 11,family="Arial"),
#   axis.title=element_text(colour="black",size = 11),
#   axis.ticks=element_line(colour="black"),
#   plot.margin=unit(c(1,1,1,1),"line")
#   ) +
# ggsave(('results/box_plot_EPOR_l2fc.svg'), width=9, height=15, units="cm")