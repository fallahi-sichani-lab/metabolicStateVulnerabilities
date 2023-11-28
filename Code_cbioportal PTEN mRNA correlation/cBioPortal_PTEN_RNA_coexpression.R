rm(list=ls()) #clear all
cat("\014") #clc

#libraries needed
library(ggplot2)
library(ggpubr)
library(purrr)
library(factoextra)
library(dplyr)
library(jcolors)
library(randomcoloR)
library(NbClust)
library(gridExtra)
library(tidyverse)
library(data.table)

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019")

#importing gene datasets
PTEN_coexpression <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/PTEN_RNA_co-expression.tsv",sep="\t")


PTEN_coexpression$corr_direction <- ifelse(PTEN_coexpression$Spearman.s.Correlation < 0, "Neg.", "Pos.")

#scatter plot of significant positive/negative correlations with PTEN mRNA levels

p <- ggplot() + geom_point(data = PTEN_coexpression[PTEN_coexpression$q.Value > 0.05,],aes(x=Spearman.s.Correlation, y=-log10(q.Value)), color= "grey80", shape = 16, size = 4, alpha = 0.2) + 
  geom_point(data=PTEN_coexpression[PTEN_coexpression$q.Value <= 0.05,], aes(x=Spearman.s.Correlation, y=-log10(q.Value), color=corr_direction), alpha = 0.7,shape = 16, size = 4) + theme_classic() + labs(x = "Spearman Correlation (rho)", y = "-log10 (q-value)", title = "PTEN mRNA co-expression correlation") +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
        strip.text=element_text(size=12), text = element_text(size =12))  +
  scale_color_manual(values=c("skyblue3", "red3"))  +
  geom_hline(yintercept = -log10(0.05), colour="black", linetype = "dashed", linewidth = 0.7) + 
  geom_vline(xintercept = 0, colour="black", linetype = "solid", linewidth = 0.7)

ggsave(file.path(outDir, 'PTEN_mRNA_coexpression_scatter_plot.pdf'), plot = p,device = 'pdf',width = 5, height = 6,dpi=300)
