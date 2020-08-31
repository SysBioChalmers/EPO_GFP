# boxplots: growth and productivity
# 
# sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

# Random number generation:
#  RNG:     Mersenne-Twister 
#  Normal:  Inversion 
#  Sample:  Rounding 
 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     

# other attached packages:
# [1] reshape2_1.4.4 dplyr_0.8.5    gdtools_0.2.2  ggplot2_3.3.0 
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
pathToEpoGfp <- '/path/to/local/epo_gfp_repo/'

library(dplyr)
library(ggplot2)
library(reshape2)

kallisto.dir <- list.files(path = paste(pathToEpoGfp,'data/kallisto_outputs/',sep=''))
nam <- substring(sapply(strsplit(kallisto.dir, "_lib"), "[", 1),10)
samples <- data.frame(replicate = c(gsub("lib.+_.+_","",substring(kallisto.dir,10))), sampleName = c(sapply(strsplit(nam, "_"),"[",1)), protein = rep(c("GFP_Control","EPO_Control","EPO","GFP"), c(2,2,13,18)), producibility = rep(c("Control","Producer"), c(4,31)), cellType = rep(c("HEK293F","HEK293freestyle","HEK293F"), c(2,15,18)), growth = rep(c(0.0254,0.0257,0.0227,0.02,0.0241,0.024,0.0223,0.0257,0.0254,0.0246,0.0218,0.0226,0.0228,0.0254,0.0243,0.0226), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), epoQP = rep(c(NaN,0,3.67,4.05,3.35,13.9,2.87,2.35,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)), gfpQP = rep(c(0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,4.31,3.86,4.29,2.27,3.15,2.48,1.16,1), c(2,2,2,2,2,3,2,2,2,2,4,2,2,2,2,2)))

samples <- distinct(samples[5:35,c('sampleName','protein','growth','epoQP','gfpQP')])
samples$productivity <- samples$epoQP
samples$productivity[7:14] <- samples$gfpQP[7:14]
samples$epoQP <- NULL
samples$gfpQP <- NULL

df <- melt(data=samples, id.vars=c('protein','sampleName'), measure.vars=c('growth', 'productivity'))

ggplot(df[which(df$variable=='growth'),], aes(x=protein, y=value, label=sampleName)) +
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
ggsave(paste(pathToEpoGfp,'results/box_plot_growth.svg',sep=''), width=6, height=15, units="cm")

ggplot(df[which(df$variable=='productivity'),], aes(x=protein, y=value, label=sampleName)) +
  geom_boxplot(aes()) +
  stat_boxplot(geom ='errorbar') + 
  geom_point() +
  geom_text(size=2.5,hjust=-0.2, vjust=0) +
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
ggsave(paste(pathToEpoGfp,'results/box_plot_productivity_absolute.svg',sep=''), width=5, height=15, units="cm")

#generating separate plots for productivity
df_gfp <- df[df$protein=='GFP',]
ggplot(df_gfp[which(df_gfp$variable=='productivity'),], aes(x=protein, y=value, label=sampleName)) +
  geom_boxplot(aes()) +
  stat_boxplot(geom ='errorbar') + 
  geom_point() +
  geom_text(size=2.5,hjust=-0.2, vjust=0) +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  panel.background = element_blank(),
  panel.border=element_rect(fill=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_blank(),
  axis.text.x=element_text(colour="black",size = 11),
  axis.text.y=element_text(colour="black",size = 11),
  axis.title=element_text(colour="black",size = 11),
  axis.ticks=element_line(colour="black"),
  plot.margin=unit(c(1,1,1,1),"line")
  ) +
ggsave(paste(pathToEpoGfp,'results/box_plot_productivity_gfp.svg',sep=''), width=5, height=15, units="cm")

df_epo <- df[df$protein=='EPO',]
ggplot(df_epo[which(df_epo$variable=='productivity'),], aes(x=protein, y=value, label=sampleName)) +
  geom_boxplot(aes()) +
  stat_boxplot(geom ='errorbar') + 
  geom_point() +
  geom_text(size=2.5,hjust=-0.2, vjust=0) +
  theme(
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  panel.background = element_blank(),
  panel.border=element_rect(fill=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_blank(),
  axis.text.x=element_text(colour="black",size = 11),
  axis.text.y=element_text(colour="black",size = 11),
  axis.title=element_text(colour="black",size = 11),
  axis.ticks=element_line(colour="black"),
  plot.margin=unit(c(1,1,1,1),"line")
  ) +
ggsave(paste(pathToEpoGfp,'results/box_plot_productivity_epo.svg',sep=''), width=5, height=15, units="cm")




