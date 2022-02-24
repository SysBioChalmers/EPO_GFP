# TFAM
# 
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***

library(biomaRt)
library(dplyr)
library(reshape)
library(tidyverse)
library(viridis)
library(ggpubr)
library(hrbrthemes)

# import samples table
sam <- read.delim(paste(getwd(),'/data/samples_for_correlation_analysis.txt',sep=''))
rownames(sam) <- sam[,1]
sam <- sam[,2:17]
samples_epo <- sam[,grep('EPO',colnames(sam))]
samples_epo['epoQP',] <- samples_epo['productivity',]
samples_gfp <- sam[,grep('GFP',colnames(sam))]
samples_gfp['gfpQP',] <- samples_gfp['productivity',]

# import TPMs
tpm <- read.delim(paste(getwd(),'/data/35samples_TPM.txt',sep=''))
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]
colnames(tpm) = rep(colnames(sam), times = c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2))

tfam = 'ENSG00000108064'

tpm_tfam = tpm[tfam,]
nam = c('CTRL293Free','EPOpoly','EPOI2','EPOB9','EPO7','EPO8','EPOF21','CTRL293F','GFPpoly','GFP3','GFP27','GFP29','GFP28','GFP25','GFP26','GFP1')

data = data.frame(name = nam, tpm_tfam=0)
data$tpm_tfam[1] = mean(as.numeric(tpm_tfam[3:4]))
data$tpm_tfam[2] = mean(as.numeric(tpm_tfam[16:17]))
data$tpm_tfam[3] = mean(as.numeric(tpm_tfam[14:15]))
data$tpm_tfam[4] = mean(as.numeric(tpm_tfam[9:10]))
data$tpm_tfam[5] = mean(as.numeric(tpm_tfam[5:6]))
data$tpm_tfam[6] = mean(as.numeric(tpm_tfam[7:8]))
data$tpm_tfam[7] = mean(as.numeric(tpm_tfam[11:13]))
data$tpm_tfam[8] = mean(as.numeric(tpm_tfam[1:2]))
data$tpm_tfam[9] = mean(as.numeric(tpm_tfam[34:35]))
data$tpm_tfam[10] = mean(as.numeric(tpm_tfam[32:33]))
data$tpm_tfam[11] = mean(as.numeric(tpm_tfam[26:27]))
data$tpm_tfam[12] = mean(as.numeric(tpm_tfam[30:31]))
data$tpm_tfam[13] = mean(as.numeric(tpm_tfam[28:29]))
data$tpm_tfam[14] = mean(as.numeric(tpm_tfam[20:21]))
data$tpm_tfam[15] = mean(as.numeric(tpm_tfam[22:25]))
data$tpm_tfam[16] = mean(as.numeric(tpm_tfam[18:19]))

data$sd[1] = sd(as.numeric(tpm_tfam[3:4]))
data$sd[2] = sd(as.numeric(tpm_tfam[16:17]))
data$sd[3] = sd(as.numeric(tpm_tfam[14:15]))
data$sd[4] = sd(as.numeric(tpm_tfam[9:10]))
data$sd[5] = sd(as.numeric(tpm_tfam[5:6]))
data$sd[6] = sd(as.numeric(tpm_tfam[7:8]))
data$sd[7] = sd(as.numeric(tpm_tfam[11:13]))
data$sd[8] = sd(as.numeric(tpm_tfam[1:2]))
data$sd[9] = sd(as.numeric(tpm_tfam[34:35]))
data$sd[10] = sd(as.numeric(tpm_tfam[32:33]))
data$sd[11] = sd(as.numeric(tpm_tfam[26:27]))
data$sd[12] = sd(as.numeric(tpm_tfam[30:31]))
data$sd[13] = sd(as.numeric(tpm_tfam[28:29]))
data$sd[14] = sd(as.numeric(tpm_tfam[20:21]))
data$sd[15] = sd(as.numeric(tpm_tfam[22:25]))
data$sd[16] = sd(as.numeric(tpm_tfam[18:19]))

data$name <-  factor(data$name , levels=data$name)

ggplot(data) +
  geom_bar( aes(x=name, y=tpm_tfam), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=tpm_tfam-sd, ymax=tpm_tfam+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) + ylab('TPM') +  xlab('tfam expression') +
  ggsave(paste(getwd(),'/results/tfam_tpm.svg',sep=''), bg = "transparent", width=20, height=10, units="cm")

data$prod[1] = sam[2,2]
data$prod[2] = sam[2,8]
data$prod[3] = sam[2,7]
data$prod[4] = sam[2,5]
data$prod[5] = sam[2,3]
data$prod[6] = sam[2,4]
data$prod[7] = sam[2,6]
data$prod[8] = sam[2,1]
data$prod[9] = sam[2,16]
data$prod[10] = sam[2,15]
data$prod[11] = sam[2,12]
data$prod[12] = sam[2,14]
data$prod[13] = sam[2,13]
data$prod[14] = sam[2,10]
data$prod[15] = sam[2,11]
data$prod[16] = sam[2,9]

data_epo = data[1:7,]
ggplot(data_epo, aes(x=tpm_tfam, y=prod)) +
  geom_point() + ggtitle('EPO producers') + ylab('productivity') + xlab('tpfm TPM') +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()

data_gfp = data[8:16,]
ggplot(data_gfp, aes(x=tpm_tfam, y=prod)) +
  geom_point() + ggtitle('GFP producers') + ylab('productivity') + xlab('tpfm TPM') +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_ipsum()


