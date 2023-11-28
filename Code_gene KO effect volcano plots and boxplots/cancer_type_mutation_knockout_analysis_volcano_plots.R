rm(list=ls()) #clear all
cat("\014") #clc

#importing libraries needed
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
library(ggrepel)

outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/Comparing mutation and cancer type analysis")

#loading mutation permutation results 
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/permutation_mutation_analysis_results.RData")

#load cancer type permutation analysis
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/pancancer_cancer_type_permutation_results.RData")


#plotting OXPHOS low vulnerability KO effects across cancer types and mutations
OXPHOS_low_candidates <- sort(unique(cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Median_sub <= -0.15 & 
                                                                                     cancer_type_pancancer_permutation_results.df$Cancer_Type == "Pancancer",]$Gene))

plot_list <- list()
plot_list <- lapply(OXPHOS_low_candidates, function(g) { 
  #getting all cancer type results
  df_cancer_type <- cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Gene == g, c(1:4)]
  df_cancer_type$Median_sub <- (df_cancer_type$Median_sub)*-1
  df_cancer_type$Cancer_or_Mutation <- c("Cancer")
  
  #getting all mutation results
  df_mutation_type <- permutation_each_mutation_results.df[permutation_each_mutation_results.df$Gene_KO == g, c(1:4)]
  pancancer_result <- cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Cancer_Type == "Pancancer" & cancer_type_pancancer_permutation_results.df$Gene == g, c(1:4)]
  colnames(pancancer_result) <-colnames(df_mutation_type)
  df_mutation_type <- rbind(df_mutation_type, pancancer_result)
  df_mutation_type$Cancer_or_Mutation <- c("Mutation")

  colnames(df_mutation_type) <- colnames(df_cancer_type)
  
  
  df_merged <- rbind(df_mutation_type, df_cancer_type)
  
  #labeling pancancer result for plot
  df_merged[df_merged$Cancer_Type== "Pancancer", ]$Cancer_or_Mutation <- "True_Pancancer"
  
  #pancancer KO effect threshold for plot
  pancancer_threshold = df_cancer_type[df_cancer_type$Cancer_Type == "Pancancer", ]$Median_sub
  
  
  
  plot_list[[g]] <- ggplot(data=df_merged, aes(x=Median_sub, y=-log10(adj_pvalue), fill = Cancer_or_Mutation)) + geom_point(color ="black",shape = 21, size = 6) + theme_classic() + labs(x = "Difference in median gene dependency scores", y = "-log10 (p-value)", title = paste(g, "KO effect")) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12),
          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
          strip.text=element_text(size=12), text = element_text(size =12),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both", type = "closed")))  +
    geom_hline(yintercept = -log10(0.01), colour="red", linetype = "dashed", linewidth = 0.7) + geom_vline(xintercept = c(-pancancer_threshold, pancancer_threshold), colour="black", linetype = "dashed", linewidth = 0.4) +
    geom_text(aes(label=ifelse(-log10(adj_pvalue)>= -log10(0.01) & abs(Median_sub) >= abs(pancancer_threshold) ,as.character(Cancer_Type),'')),hjust=0.2,vjust=-1, angle = 30, size = 4) + xlim(-1,1) + 
    ylim(0,5) + scale_fill_manual(values=c("#E69F00", "#56B4E9", "black")) +
    geom_point(data = df_merged[abs(df_merged$Median_sub) < abs(pancancer_threshold) |df_merged$adj_pvalue > 0.01,],aes(x=Median_sub, y=-log10(adj_pvalue)), color ="black",shape = 21, size = 6, fill = "grey") +
    geom_segment(aes(x=0, y=0, xend=0, yend=5), arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"), lwd=0.3) + 
    scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,5), expand = c(0, 0)) 
  
  
  }
)

p <- ggarrange(plotlist = plot_list)

ggsave(file.path(outDir, 'cancer_type_versus_mutation_OXPHOS_low_KO_effect.pdf'), plot = p,device = 'pdf',width =20, height = 20,dpi=300)



#plotting OXPHOS high vulnerability KO effects across cancer types and mutations
OXPHOS_high_candidates <- sort(unique(cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Median_sub >= 0.15 & 
                                                                                   cancer_type_pancancer_permutation_results.df$Cancer_Type == "Pancancer",]$Gene))
plot_list <- list()
plot_list <- lapply(OXPHOS_high_candidates, function(g) {

    
    #getting all cancer type results
    df_cancer_type <- cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Gene == g, c(1:4)]
    df_cancer_type$Cancer_or_Mutation <- c("Cancer")
    
    #getting all mutation results
    df_mutation_type <- permutation_each_mutation_results.df[permutation_each_mutation_results.df$Gene_KO == g, c(1:4)]
    pancancer_result <- df_cancer_type[df_cancer_type$Cancer_Type == "Pancancer", (1:4)]
    colnames(pancancer_result) <-colnames(df_mutation_type)
    df_mutation_type <- rbind(df_mutation_type, pancancer_result)
    df_mutation_type$Cancer_or_Mutation <- c("Mutation")
    df_mutation_type$Median_sub <- (df_mutation_type$Median_sub)*-1
    colnames(df_mutation_type) <- colnames(df_cancer_type)
    
    
    df_merged <- rbind(df_mutation_type, df_cancer_type)
    
    #labeling pancancer result for plot
    df_merged[df_merged$Cancer_Type== "Pancancer", ]$Cancer_or_Mutation <- "True_Pancancer"
    
    #pancancer KO effet for plot
    pancancer_threshold = df_cancer_type[df_cancer_type$Cancer_Type == "Pancancer", ]$Median_sub
    

  
  plot_list[[g]] <- ggplot(data=df_merged, aes(x=Median_sub, y=-log10(adj_pvalue), fill = Cancer_or_Mutation)) + geom_point(color ="black",shape = 21, size = 6) + theme_classic() + labs(x = "Difference in median gene dependency scores", y = "-log10 (p-value)", title = paste(g, "KO effect")) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12),
          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
          strip.text=element_text(size=12), text = element_text(size =12),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both", type = "closed")))  +
    geom_hline(yintercept = -log10(0.01), colour="red", linetype = "dashed", linewidth = 0.7) + geom_vline(xintercept = c(-pancancer_threshold, pancancer_threshold), colour="black", linetype = "dashed", linewidth = 0.4) +
    geom_text(aes(label=ifelse(-log10(adj_pvalue)>= -log10(0.01) & abs(Median_sub) >= abs(pancancer_threshold) ,as.character(Cancer_Type),'')),hjust=0.2,vjust=-1, angle = 30, size = 4) + xlim(-1,1) + 
    ylim(0,5) + scale_fill_manual(values=c("#E69F00", "#56B4E9", "black")) +
    geom_point(data = df_merged[abs(df_merged$Median_sub) < abs(pancancer_threshold) |df_merged$adj_pvalue > 0.01,],aes(x=Median_sub, y=-log10(adj_pvalue)), color ="black",shape = 21, size = 6, fill = "grey") +
    geom_segment(aes(x=0, y=0, xend=0, yend=5), arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last", type = "closed"), lwd=0.3) + 
    scale_x_continuous(limits = c(-1,1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,5), expand = c(0, 0)) 
  }

)

p <- ggarrange(plotlist = plot_list)

ggsave(file.path(outDir, 'cancer_type_versus_mutation_OXPHOS_high_KO_effect.pdf'), plot = p,device = 'pdf',width =20, height = 15,dpi=300)













 