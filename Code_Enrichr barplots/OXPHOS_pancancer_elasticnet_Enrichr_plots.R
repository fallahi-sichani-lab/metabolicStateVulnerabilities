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
library(ggplot2)
library(gridExtra)
library(fgsea)
library(ggrepel)
library(randomcoloR)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/Enrichr analysis")

#loading enrichr results
GO_CC_OXPHOS_high <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/Enrichr analysis/OXPHOS_high_pancancer_GO_Cellular_Component_2023_table.txt", sep = '\t', header = TRUE)
GO_CC_OXPHOS_low <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/Enrichr analysis/OXPHOS_low_pancancer_GO_Cellular_Component_2023_table.txt", sep = '\t', header = TRUE)

Reactome_OXPHOS_high <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/Enrichr analysis/OXPHOS_high_pancancer_Reactome_2022_table.txt", sep = '\t', header = TRUE)
Reactome_OXPHOS_low <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/Enrichr analysis/OXPHOS_low_pancancer_Reactome_2022_table.txt", sep = '\t', header = TRUE)


#keeping GO Cellular Component terms with p-values <= 0.05
GO_CC_OXPHOS_high <- GO_CC_OXPHOS_high[GO_CC_OXPHOS_high$P.value<= 0.05, ]
GO_CC_OXPHOS_low <- GO_CC_OXPHOS_low[GO_CC_OXPHOS_low$P.value <= 0.05, ]


#GO CC OXPHOS high barplot
GO_CC_OXPHOS_high <- GO_CC_OXPHOS_high[order(GO_CC_OXPHOS_high$P.value),] #ordering coefficient values
GO_CC_OXPHOS_high$Term <- factor(GO_CC_OXPHOS_high$Term, levels = GO_CC_OXPHOS_high$Term )


p <- ggplot(GO_CC_OXPHOS_high[c(1:20),], aes(x=-log10(P.value), y=reorder(Term, -P.value))) + geom_bar(stat="identity", fill="#C90D2E", width = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth =0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(p-value)") + 
  labs(title = "OXPHOS high enriched GO CC terms (p <= 0.05, top 20 terms)", size = 14) + xlim(0,9) + 
  scale_x_continuous(limits = c(0,9), expand = c(0, 0)) 
ggsave(file.path(outDir, 'OXPHOS_high_GO_CC_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 5,dpi=300) 



#GO CC OXPHOS low barplot
GO_CC_OXPHOS_low <- GO_CC_OXPHOS_low[order(GO_CC_OXPHOS_low$P.value),] #ordering coefficient values
GO_CC_OXPHOS_low$Term <- factor(GO_CC_OXPHOS_low$Term, levels = GO_CC_OXPHOS_low$Term )

p <- ggplot(GO_CC_OXPHOS_low[c(1:20),], aes(x=-log10(P.value), y=reorder(Term, -P.value))) + geom_bar(stat="identity", fill="#426ED9", width = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(p-value)") + 
  labs(title = "OXPHOS low enriched GO CC terms (p <= 0.05), top 20 terms", size = 14) + xlim(0,9) + 
  scale_x_continuous(limits = c(0,9), expand = c(0, 0)) 
ggsave(file.path(outDir, 'OXPHOS_low_GO_CC_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 10, height = 5,dpi=300)



#keeping Reactome terms
Reactome_OXPHOS_high <- Reactome_OXPHOS_high[Reactome_OXPHOS_high$Adjusted.P.value <= 0.05, ]
Reactome_OXPHOS_low <- Reactome_OXPHOS_low[Reactome_OXPHOS_low$Adjusted.P.value <= 0.05, ]


#Reactome OXPHOS high barplot
Reactome_OXPHOS_high <- Reactome_OXPHOS_high[order(Reactome_OXPHOS_high$Adjusted.P.value),] #ordering coefficient values
Reactome_OXPHOS_high$Term <- factor(Reactome_OXPHOS_high$Term, levels = Reactome_OXPHOS_high$Term )


p <- ggplot(Reactome_OXPHOS_high, aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="#C90D2E", width = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth =0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(adj. p-value)") + 
  labs(title = "OXPHOS high enriched Reactome terms (adj. p-value <= 0.05)", size = 14) + xlim(0,3) + 
  scale_x_continuous(limits = c(0,3), expand = c(0, 0)) 
ggsave(file.path(outDir, 'OXPHOS_high_Reactome_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 15, height = 8,dpi=300) 


#Reactome OXPHOS low barplot
Reactome_OXPHOS_low <- Reactome_OXPHOS_low[order(Reactome_OXPHOS_low$Adjusted.P.value),] #ordering coefficient values
Reactome_OXPHOS_low$Term <- factor(Reactome_OXPHOS_low$Term, levels = Reactome_OXPHOS_low$Term )

p <- ggplot(Reactome_OXPHOS_low, aes(x=-log10(Adjusted.P.value), y=reorder(Term, -Adjusted.P.value))) + geom_bar(stat="identity", fill="#426ED9", width = 0.7)+theme_classic()  + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(linewidth=0.2,color="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("") + xlab("-log10(adj. p-value)") + 
  labs(title = "OXPHOS low enriched Reactome terms (adj. p-value <= 0.05)", size = 14) + xlim(0,3) + 
  scale_x_continuous(limits = c(0,3), expand = c(0, 0)) 

ggsave(file.path(outDir, 'OXPHOS_low_Reactome_enriched_terms_barplot.pdf'), plot = p,device = 'pdf',width = 12, height = 3,dpi=300)





