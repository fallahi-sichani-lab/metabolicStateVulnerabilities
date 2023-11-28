rm(list=ls()) #clear all
cat("\014") #clc

#importing libraries needed
library(viridis)
library(stringr)
library(plyr)
library(dplyr)
library(gtools)
library(matrixStats)
library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)
library(tidyr)
library(devtools)
library(gridExtra)
library(ggrepel)
library(randomcoloR)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/Enrichr analysis")

#loading enrichr results
GO_BP_PTEN_negative <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/Enrichr analysis/Negative_PTEN_coexpression_GO_Biological_Process_2023_table.txt", sep = '\t', header = TRUE)
GO_CC_PTEN_negative<- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/Enrichr analysis/Negative_PTEN_coexpression_GO_Cellular_Component_2023_table.txt", sep = '\t', header = TRUE)


GO_BP_PTEN_positive <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/Enrichr analysis/Positive_PTEN_coexpression_GO_Biological_Process_2023_table.txt", sep = '\t', header = TRUE)
GO_CC_PTEN_positive <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/cBioPortal/GLASS Consortium, Nature 2019/Enrichr analysis/Positive_PTEN_coexpression_GO_Cellular_Component_2023_table.txt", sep = '\t', header = TRUE)


#keeping GO terms with FDR <= 0.05
GO_BP_PTEN_negative <- GO_BP_PTEN_negative[GO_BP_PTEN_negative$Adjusted.P.value<= 0.05, ]
GO_CC_PTEN_negative <- GO_CC_PTEN_negative[GO_CC_PTEN_negative$Adjusted.P.value<= 0.05, ]

GO_BP_PTEN_positive <- GO_BP_PTEN_positive[GO_BP_PTEN_positive$Adjusted.P.value<= 0.05, ]
GO_CC_PTEN_positive <- GO_CC_PTEN_positive[GO_CC_PTEN_positive$Adjusted.P.value<= 0.05, ]


#barplot of GO biological processes enriched in genes negatively correlated with PTEN
p <- ggplot(GO_BP_PTEN_negative, aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="black", width = 0.7, alpha = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth =0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(FDR)") +
  labs(title = "GO Biological Process (FDR <= 0.05)", size = 14) + 
  scale_x_continuous(limits = c(0,5), expand = c(0, 0)) 

ggsave(file.path(outDir, 'PTEN_negative_coexpression_GOBP_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 5,dpi=300) 

#barplot of GO cellular component enriched in genes negatively correlated with PTEN
p <- ggplot(GO_CC_PTEN_negative, aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="black", width = 0.7, alpha = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(FDR)") + 
  labs(title = "GO Cellular Component (FDR <= 0.05)", size = 14) +
  scale_x_continuous(limits = c(0,15), expand = c(0, 0)) 
ggsave(file.path(outDir, 'PTEN_negative_coexpression_GCC_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 5,dpi=300)


#barplot of GO biological processes enriched in genes positively correlated with PTEN
p <- ggplot(na.omit(GO_BP_PTEN_positive), aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="black", width = 0.7, alpha = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth =0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(FDR)") +
  labs(title = "GO Biological Process (FDR <= 0.05)", size = 14) + 
  scale_x_continuous(limits = c(0,5), expand = c(0, 0)) 

ggsave(file.path(outDir, 'PTEN_positive_coexpression_GOBP_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 8,dpi=300) 

#barplot of GO cellular component enriched in genes positively correlated with PTEN
p <- ggplot(GO_CC_PTEN_positive, aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="black", width = 0.7, alpha = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(FDR)") + 
  labs(title = "GO Cellular Component (FDR <= 0.05)", size = 14) +
  scale_x_continuous(limits = c(0,10), expand = c(0, 0)) 
ggsave(file.path(outDir, 'PTEN_positive_coexpression_GCC_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 5,dpi=300)

