# IPAheatmap: pairwise comparison of producers vs. control
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
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     

# other attached packages:
#  [1] RColorBrewer_1.1-2 ggpubr_0.3.0       cowplot_1.0.0     
#  [4] svglite_1.2.3      scales_1.1.1       piano_2.4.0       
#  [7] forcats_0.5.0      stringr_1.4.0      purrr_0.3.4       
# [10] readr_1.3.1        tidyr_1.1.0        tibble_3.0.1      
# [13] tidyverse_1.3.0    ggrepel_0.8.2      dplyr_0.8.5       
# [16] reshape2_1.4.4     ggplot2_3.3.0   
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
pathToEpoGfp <- '/path/to/local/epo_gfp_repo/'

library("ggplot2")
library("reshape2")
library("dplyr")
library("ggrepel")
library("tidyverse")
library("piano")
library("scales")
library("svglite")
library("cowplot")
library("ggpubr")
library('RColorBrewer')

epo <- read.delim(paste(pathToEpoGfp,'results/epo_vs_cnt_pairwise_IPA_l2fc_0.58_padj_0.05.txt',sep=''))
gfp <- read.delim(paste(pathToEpoGfp,'results/gfp_vs_cnt_pairwise_IPA_l2fc_0.58_padj_0.05.txt',sep=''))

all <- merge(epo, gfp, by.x=1, by.y=1, all.x=TRUE, all.y=TRUE) #merge epo & gfp
all <- all[ ,c(1,5,3,2,4,6,7,9,8,11,12,10,13)] #make order of columns the same order as productivity
all[is.na(all)] <- 0
mix_fig <- data.frame()
for (k in colnames(all)) mix_fig[[k]] <- as.character()

for (i in 2:dim(all)[2]) {
  all <- all[order(-all[,i], all[,1]),]
  mix_fig <- rbind(mix_fig, all[all[,i]>1.3,])
}
mix_fig <- unique(mix_fig)

#------------------------------------------------------ df table
df <- mix_fig
cons <- mix_fig
cons[cons < 1.3] <- 0
cons[cons > 1.3] <- 1
cons$sum <- rowSums(cons[,2:13])

df <- mix_fig[rowSums(cons[,2:6])>2|rowSums(cons[,7:13])>3,]
df[,2:13] <- 10^(-df[,2:13])
df <- df[c(7,2,1,4,5,3,6,8,9,10),]
melted_df <- melt(df, id = "Canonical.Pathways")
melted_df$Canonical.Pathways <-  factor(melted_df$Canonical.Pathways, levels=df$Canonical.Pathways)

ggplot(data = melted_df, aes(x = variable, y = Canonical.Pathways)) + 
  geom_tile(aes(fill = value), colour = "black") + 
  coord_flip() +
  theme(
    aspect.ratio = 1/2,
    axis.text.y = element_text(angle = 0, hjust = 1, size=8), 
    axis.text.x = element_text(angle = 30, hjust = 0, size=10),
    plot.margin = unit(c(0,5,0,0),"cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_text(size=13), 
    legend.text=element_text(size=12),
    ) + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(labels = wrap_format(70), position = "right") +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 9, name = "YlOrRd")),na.value = "lightgoldenrodyellow",limits=c(0,0.1)) +
ggsave(paste(pathToEpoGfp,'results/heatmap_pro_vs_cnt_pairwise.svg',sep=''), bg = "transparent", width=20, height=10, units="cm")

