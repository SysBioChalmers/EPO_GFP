# scatterPlot: TPM vs. Productivity
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
#  [1] ggpubr_0.3.0                cowplot_1.0.0              
#  [3] svglite_1.2.3               scales_1.1.1               
#  [5] piano_2.4.0                 forcats_0.5.0              
#  [7] stringr_1.4.0               purrr_0.3.4                
#  [9] readr_1.3.1                 tidyr_1.1.0                
# [11] tibble_3.0.1                tidyverse_1.3.0            
# [13] ggrepel_0.8.2               dplyr_0.8.5                
# [15] reshape2_1.4.4              ggplot2_3.3.0              
# [17] biomaRt_2.44.0              DESeq2_1.28.1              
# [19] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
# [21] matrixStats_0.56.0          Biobase_2.48.0             
# [23] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
# [25] IRanges_2.22.2              S4Vectors_0.26.1           
# [27] BiocGenerics_0.34.0 
# 
# Rasool Saghaleyni   2020-08-25
####################################################################

# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***

pathToEpoGfp <- '/Users/rasools/drive/projects/EPO_GFP/EPO_GFP_github/'

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

load(paste(pathToEpoGfp,'results/F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
epo8 <- GSAsummaryTable(F21_vs_epo8_wilcoxon_nperm_e4_extended_secretion_gene_sets)
load(paste(pathToEpoGfp,'results/F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
epo7 <- GSAsummaryTable(F21_vs_epo7_wilcoxon_nperm_e4_extended_secretion_gene_sets)
load(paste(pathToEpoGfp,'results/F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
epob9 <- GSAsummaryTable(F21_vs_epob9_wilcoxon_nperm_e4_extended_secretion_gene_sets)
load(paste(pathToEpoGfp,'results/F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
epoi2 <- GSAsummaryTable(F21_vs_epoi2_wilcoxon_nperm_e4_extended_secretion_gene_sets)
load(paste(pathToEpoGfp,'results/F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets.RData',sep=''))
epopoly <- GSAsummaryTable(F21_vs_epopoly_wilcoxon_nperm_e4_extended_secretion_gene_sets)

#epo8
sig_epo8 <- epo8[,c("Name", "p adj (dist.dir.dn)", "p adj (non-dir.)", "p adj (dist.dir.up)")]
colnames(sig_epo8) <- c("Name", "down", "non", "up")
co=0.05
sig_epo8 <- sig_epo8[sig_epo8$non<co,]
sig_epo8 <- sig_epo8[sig_epo8$down<co|sig_epo8$up<co,]
sig_epo8 <- as.data.frame(sig_epo8)
#epo7
sig_epo7 <- epo7[,c("Name", "p adj (dist.dir.dn)", "p adj (non-dir.)", "p adj (dist.dir.up)")]
colnames(sig_epo7) <- c("Name", "down", "non", "up")
co=0.05
sig_epo7 <- sig_epo7[sig_epo7$non<co,]
sig_epo7 <- sig_epo7[sig_epo7$down<co|sig_epo7$up<co,]
sig_epo7 <- as.data.frame(sig_epo7)
#epob9
sig_epob9 <- epob9[,c("Name", "p adj (dist.dir.dn)", "p adj (non-dir.)", "p adj (dist.dir.up)")]
colnames(sig_epob9) <- c("Name", "down", "non", "up")
co=0.05
sig_epob9 <- sig_epob9[sig_epob9$non<co,]
sig_epob9 <- sig_epob9[sig_epob9$down<co|sig_epob9$up<co,]
sig_epob9 <- as.data.frame(sig_epob9)
#epoi2
sig_epoi2 <- epoi2[,c("Name", "p adj (dist.dir.dn)", "p adj (non-dir.)", "p adj (dist.dir.up)")]
colnames(sig_epoi2) <- c("Name", "down", "non", "up")
co=0.05
sig_epoi2 <- sig_epoi2[sig_epoi2$non<co,]
sig_epoi2 <- sig_epoi2[sig_epoi2$down<co|sig_epoi2$up<co,]
sig_epoi2 <- as.data.frame(sig_epoi2)
#epopoly
sig_epopoly <- epopoly[,c("Name", "p adj (dist.dir.dn)", "p adj (non-dir.)", "p adj (dist.dir.up)")]
colnames(sig_epopoly) <- c("Name", "down", "non", "up")
co=0.05
sig_epopoly <- sig_epopoly[sig_epopoly$non<co,]
sig_epopoly <- sig_epopoly[sig_epopoly$down<co|sig_epopoly$up<co,]
sig_epopoly <- as.data.frame(sig_epopoly)

#epo8
sig_epo8 <- sig_epo8[order(sig_epo8$non),]
sig_epo8 <- sig_epo8[order(sig_epo8$up),]
top_path <- sig_epo8[1:5,]
sig_epo8 <- sig_epo8[order(sig_epo8$non),]
sig_epo8 <- sig_epo8[order(sig_epo8$down),]
top_path <- rbind(top_path,sig_epo8[1:5,])
#epo7
sig_epo7 <- sig_epo7[order(sig_epo7$non),]
sig_epo7 <- sig_epo7[order(sig_epo7$up),]
top_path <- rbind(top_path,sig_epo7[1:5,])
sig_epo7 <- sig_epo7[order(sig_epo7$non),]
sig_epo7 <- sig_epo7[order(sig_epo7$down),]
top_path <- rbind(top_path,sig_epo7[1:5,])
#epob9
sig_epob9 <- sig_epob9[order(sig_epob9$non),]
sig_epob9 <- sig_epob9[order(sig_epob9$up),]
top_path <- rbind(top_path,sig_epob9[1:5,])
#epoi2
sig_epoi2 <- sig_epoi2[order(sig_epoi2$non),]
sig_epoi2 <- sig_epoi2[order(sig_epoi2$up),]
top_path <- rbind(top_path,sig_epoi2[1:5,])
#epopoly
# sig_epopoly <- sig_epopoly[order(sig_epopoly$non),]
# sig_epopoly <- sig_epopoly[order(sig_epopoly$up),]
# top_path <- rbind(top_path,sig_epopoly[1:10,])
# sig_epopoly <- sig_epopoly[order(sig_epopoly$non),]
# sig_epopoly <- sig_epopoly[order(sig_epopoly$down),]
# top_path <- rbind(top_path,sig_epopoly[1:10,])
#------------------------------------------------------ mix_fig table
mix_fig <- matrix(nrow=length(unique(top_path$Name)), ncol=5)
colnames(mix_fig) <- rep(c("EPO8","EPO7","EPOB9","EPOI2","EPOpoly"),1)
rownames(mix_fig) <- unique(top_path$Name)
#------------------------------------------------------
epo8_topPath <- epo8[epo8$Name %in% rownames(mix_fig),]
epo8_topPath <- epo8_topPath[,c("Name", "p adj (dist.dir.dn)", "p adj (dist.dir.up)","p adj (non-dir.)")]
colnames(epo8_topPath) <- c("Name","dn","up","non")
epo8_topPath <- epo8_topPath[match(rownames(mix_fig), epo8_topPath$Name),]

for (i in 1:length(rownames(mix_fig))) {
    a <- which.min(epo8_topPath[i,2:3])
    if (a==1){
    mix_fig[i,1] <- log(epo8_topPath$dn[i],10)
    }
    if (a==2){
    mix_fig[i,1] <- -log(epo8_topPath$up[i],10)
    }
    if (epo8_topPath$non[i]>co){
    mix_fig[i,1] <- 0
    }
}
#--------
epo7_topPath <- epo7[epo7$Name %in% rownames(mix_fig),]
epo7_topPath <- epo7_topPath[,c("Name", "p adj (dist.dir.dn)", "p adj (dist.dir.up)","p adj (non-dir.)")]
colnames(epo7_topPath) <- c("Name","dn","up","non")
epo7_topPath <- epo7_topPath[match(rownames(mix_fig), epo7_topPath$Name),]

for (i in 1:length(rownames(mix_fig))) {
    a <- which.min(epo7_topPath[i,2:3])
    if (a==1){
    mix_fig[i,2] <- log(epo7_topPath$dn[i],10)
    }
    if (a==2){
    mix_fig[i,2] <- -log(epo7_topPath$up[i],10)
    }
    if (epo7_topPath$non[i]>co){
    mix_fig[i,2] <- 0
    }
}
#--------
epob9_topPath <- epob9[epob9$Name %in% rownames(mix_fig),]
epob9_topPath <- epob9_topPath[,c("Name", "p adj (dist.dir.dn)", "p adj (dist.dir.up)","p adj (non-dir.)")]
colnames(epob9_topPath) <- c("Name","dn","up","non")
epob9_topPath <- epob9_topPath[match(rownames(mix_fig), epob9_topPath$Name),]

for (i in 1:length(rownames(mix_fig))) {
    a <- which.min(epob9_topPath[i,2:3])
    if (a==1){
    mix_fig[i,3] <- log(epob9_topPath$dn[i],10)
    }
    if (a==2){
    mix_fig[i,3] <- -log(epob9_topPath$up[i],10)
    }
    if (epob9_topPath$non[i]>co){
    mix_fig[i,3] <- 0
    }
}
#--------
epoi2_topPath <- epoi2[epoi2$Name %in% rownames(mix_fig),]
epoi2_topPath <- epoi2_topPath[,c("Name", "p adj (dist.dir.dn)", "p adj (dist.dir.up)","p adj (non-dir.)")]
colnames(epoi2_topPath) <- c("Name","dn","up","non")
epoi2_topPath <- epoi2_topPath[match(rownames(mix_fig), epoi2_topPath$Name),]

for (i in 1:length(rownames(mix_fig))) {
    a <- which.min(epoi2_topPath[i,2:3])
    if (a==1){
    mix_fig[i,4] <- log(epoi2_topPath$dn[i],10)
    }
    if (a==2){
    mix_fig[i,4] <- -log(epoi2_topPath$up[i],10)
    }
    if (epoi2_topPath$non[i]>co){
    mix_fig[i,4] <- 0
    }
}
#--------
# epopoly_topPath <- epopoly[epopoly$Name %in% rownames(mix_fig),]
# epopoly_topPath <- epopoly_topPath[,c("Name", "p adj (dist.dir.dn)", "p adj (dist.dir.up)","p adj (non-dir.)")]
# colnames(epopoly_topPath) <- c("Name","dn","up","non")
# epopoly_topPath <- epopoly_topPath[match(rownames(mix_fig), epopoly_topPath$Name),]

# for (i in 1:length(rownames(mix_fig))) {
#     a <- which.min(epopoly_topPath[i,2:3])
#     if (a==1){
#     mix_fig[i,5] <- log(epopoly_topPath$dn[i],10)
#     }
#     if (a==2){
#     mix_fig[i,5] <- -log(epopoly_topPath$up[i],10)
#     }
#     if (epopoly_topPath$non[i]>co){
#     mix_fig[i,5] <- 0
#     }
# }

mix_fig <- as.data.frame(mix_fig)
mix_fig$EPOpoly <- NULL

write.table(mix_fig, 'results/heatmap_mix_f21_goSlimSecretion_2.txt', 
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep='\t') 

df <- mix_fig
df$Name <- rownames(df)
df$mean <- rowMeans(df[,1:length(df)-1])
df <- df[order(df$mean, decreasing = FALSE),]
df$mean <- NULL
melted_df <- melt(df, id = "Name")
melted_df$Name <- gsub('GO_', '', melted_df$Name)
melted_df$Name <- gsub('_', ' ', melted_df$Name)
melted_df$Name <-  factor(melted_df$Name, levels = unique(melted_df$Name))

ggplot(data = melted_df, aes(x = variable, y = Name)) + 
  geom_tile(aes(fill = value), colour = "black") + 
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0, size=4), 
    axis.text.y = element_text(angle = 0, hjust = 1, size=6),
    plot.margin = unit(c(0,0,0.6,0.2),"cm"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_text(size=13), 
    legend.text=element_text(size=12),
    aspect.ratio=1/6,
    ) + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(labels = wrap_format(40), position = "right") +
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0) +
ggsave(paste(pathToEpoGfp,'results/heatmap_mix_f21_goSlimSecretion_2.svg', sep=''), bg = "transparent", width=17, height=8, units="cm")
