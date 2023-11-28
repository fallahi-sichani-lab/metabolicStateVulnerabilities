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

#loading KO results across cancer types and mutations
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/driver_mutations_OXPHOS_variable_cell_lines.RData")

#loading mutation permutation results 
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/permutation_mutation_analysis_results.RData")

#load cancer type permutation analysis
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/pancancer_cancer_type_permutation_results.RData")
pancancer_sig_vulnerabilities <- cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Cancer_Type == "Pancancer" & 
                                                                                cancer_type_pancancer_permutation_results.df$adj_pvalue <= 0.01 &
                                                                                abs(cancer_type_pancancer_permutation_results.df$Median_sub) >= 0.15,]


#boxplots of OXPHOS high vulnerabilities
OXPHOS_high_vulnerabilites <- c("DNM1L", "EP300")

plot_list <- list()
plot_list <- lapply(OXPHOS_high_vulnerabilites, function(g) { 
  
  #getting pancancer difference in median results
  pancancer_median_diff <- pancancer_sig_vulnerabilities[pancancer_sig_vulnerabilities$Gene == g, ]$Median_sub
  
  #getting mutations and cancer types' differences that are greater than or equal to pancancer difference
  sig_cancer_types <-cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Gene == g &
                                                                     cancer_type_pancancer_permutation_results.df$adj_pvalue <= 0.01 &
                                                                     cancer_type_pancancer_permutation_results.df$Median_sub >= pancancer_median_diff,c(1,3,4:6)]
  
  sig_mutations <- permutation_each_mutation_results.df[permutation_each_mutation_results.df$Gene_KO == g &
                                                                 permutation_each_mutation_results.df$adj_pvalue <= 0.01 &
                                                                 permutation_each_mutation_results.df$Median_sub >= pancancer_median_diff,-c(2)]
  
  colnames(sig_mutations) <- colnames(sig_cancer_types)
  
  #combining significant cancer types and mutations for labels
  sig_mutations_cancer_types <- rbind(sig_cancer_types, sig_mutations)
  
  
  #extracting pancancer cell lines
  pancancer_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("Cancer_Type", "OXPHOS_state", g)]
  pancancer_KO_results$Cancer_Type <- c("Pancancer")
  
  #extracting individual cancer type cell lines
  cancer_type_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0 &
                                                                           mutation_OXPHOS_variable_cell_lines_final.df$Cancer_Type %in% sig_cancer_types$Cancer_Type[-(3)],c("Cancer_Type","OXPHOS_state", g)]
  
  #combining pancancer and cancer type KO results
  cancer_type_pancancer_KO_results <- rbind(cancer_type_KO_results, pancancer_KO_results)
  colnames(cancer_type_pancancer_KO_results)[3] <- c("Gene_KO")
  
  
  #extracting mutated cell lines
  mutation_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("StrippedCellLineName","OXPHOS_state", g, sig_mutations$Cancer_Type)]
  
  mutation_KO_results_for_plotting <- data.frame()
  for (i in 1:length(sig_mutations$Cancer_Type)) {
    
    mutation <- sig_mutations$Cancer_Type[i]
    df <- mutation_KO_results[, c(1:3, 3+i)]
    colnames(df)[4] <- "Mutation"
    mutated_in_gene <- na.omit(df[df$Mutation == 1, ])
    mutated_in_gene$Mutation_status <- mutation
    
    mutation_KO_results_for_plotting <- rbind(mutation_KO_results_for_plotting,mutated_in_gene[,-c(1,4)])

  }
  
  #rearranging dataframe so we can bind it with cancer/pancancer results
  mutation_KO_results_for_plotting <- mutation_KO_results_for_plotting[,c(3,1,2)]
  colnames(mutation_KO_results_for_plotting) <- colnames(cancer_type_pancancer_KO_results)
  
  #combining mutation, cancer type,and pancancer results
  mutation_cancer_type_pancancer_KO_results <- rbind(cancer_type_pancancer_KO_results,mutation_KO_results_for_plotting)
  mutation_cancer_type_pancancer_KO_results <- merge(mutation_cancer_type_pancancer_KO_results,sig_mutations_cancer_types[,c(1:5)], by = "Cancer_Type")
  
  #converting table from wide to long for plotting purposes, and using boxplots to plot Gene_dependency scores in OXPHOS high vs. OXPHOS low cell lines
  mutation_cancer_type_pancancer_KO_results$OXPHOS_state <- factor(mutation_cancer_type_pancancer_KO_results$OXPHOS_state, levels=c("OXPHOS_low", "OXPHOS_high"))
  mutation_cancer_type_pancancer_KO_results$Cancer_Type <- factor(mutation_cancer_type_pancancer_KO_results$Cancer_Type, levels=c("Pancancer", sig_mutations$Cancer_Type, sig_cancer_types[sig_cancer_types$Cancer_Type != "Pancancer", ]$Cancer_Type))
  
  
  #generating plots
  plot_list[[g]] <- ggplot(mutation_cancer_type_pancancer_KO_results, aes(x=Cancer_Type, y= Gene_KO, fill=OXPHOS_state)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + labs(title = paste(g, "KO effects across cancers and mutations")) + xlab("") + ylab("Gene dependency") +
    theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(),
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=12),
                            axis.line=element_line(linewidth=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",linewidth=0.5),
                            axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12)) + ylim(0,1.2) +
    geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 1, alpha = 0.6) + scale_fill_manual(values=c("skyblue1", "tomato")) + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
              aes(label = paste("Diff.  = ", round(Median_sub, digits = 2)), x = Cancer_Type, y = 1.15), size = 3, color = "red") +
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
             aes(label = paste("p = ", format(adj_pvalue, scientific = TRUE, digits = 3)), x = Cancer_Type, y = 1.2), size = 3, color = "red") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("low = ", num_OXPHOS_low_cell_lines), x = Cancer_Type, y = 1.08), size = 3, color = "black") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("high = ", num_OXPHOS_high_cell_lines), x = Cancer_Type, y = 1.12), size = 3, color = "black")
  }
)

p <- ggarrange(plotlist = plot_list, ncol = 1)

ggsave(file.path(outDir, 'OXPHOS_high_vulnerabilities_boxplots.pdf'), plot = p,device = 'pdf',width =15, height = 15,dpi=300)


#boxplots of ETC/mitochondrial related targets and the KO effect comparing pancancer vs. PTEN
ETC_vulnerabilites <- c("COQ4", "MRPL58", "NDUFB8", "NDUFC2", "NDUFS1")

plot_list <- list()
plot_list <- lapply(ETC_vulnerabilites, function(g) { 
  
  
  #getting mutations and cancer types' differences that are greater than or equal to pancancer difference
  sig_cancer_types <-cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Gene == g &
                                                                    cancer_type_pancancer_permutation_results.df$Cancer_Type == "Pancancer",c(1,3,4:6)]
  
  sig_mutations <- permutation_each_mutation_results.df[permutation_each_mutation_results.df$Gene_KO == g &
                                                          permutation_each_mutation_results.df$Mutation == "PTEN",-c(2)]
  
  colnames(sig_mutations) <- colnames(sig_cancer_types)
  
  #combining significant cancer types and mutations for labels
  sig_mutations_cancer_types <- rbind(sig_cancer_types, sig_mutations)
  
  
  #extracting pancancer cell lines
  pancancer_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("Cancer_Type", "OXPHOS_state", g)]
  pancancer_KO_results$Cancer_Type <- c("Pancancer")
  colnames(pancancer_KO_results)[3] <- c("Gene_KO")
  
  
  #extracting mutated cell lines
  mutation_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("StrippedCellLineName","OXPHOS_state", g, sig_mutations$Cancer_Type)]
  
  mutation_KO_results_for_plotting <- data.frame()
  for (i in 1:length(sig_mutations$Cancer_Type)) {
    
    mutation <- sig_mutations$Cancer_Type[i]
    df <- mutation_KO_results[, c(1:3, 3+i)]
    colnames(df)[4] <- "Mutation"
    mutated_in_gene <- na.omit(df[df$Mutation == 1, ])
    mutated_in_gene$Mutation_status <- mutation
    
    mutation_KO_results_for_plotting <- rbind(mutation_KO_results_for_plotting,mutated_in_gene[,-c(1,4)])
    
  }
  
  #rearranging dataframe so we can bind it with cancer/pancancer results
  mutation_KO_results_for_plotting <- mutation_KO_results_for_plotting[,c(3,1,2)]
  colnames(mutation_KO_results_for_plotting) <- colnames(pancancer_KO_results)
  
  #combining mutation, cancer type,and pancancer results
  mutation_cancer_type_pancancer_KO_results <- rbind(pancancer_KO_results,mutation_KO_results_for_plotting)
  mutation_cancer_type_pancancer_KO_results <- merge(mutation_cancer_type_pancancer_KO_results,sig_mutations_cancer_types[,c(1:5)], by = "Cancer_Type")
  
  #converting table from wide to long for plotting purposes, and using boxplots to plot Gene_dependency scores in OXPHOS high vs. OXPHOS low cell lines
  mutation_cancer_type_pancancer_KO_results$OXPHOS_state <- factor(mutation_cancer_type_pancancer_KO_results$OXPHOS_state, levels=c("OXPHOS_low", "OXPHOS_high"))
  mutation_cancer_type_pancancer_KO_results$Cancer_Type <- factor(mutation_cancer_type_pancancer_KO_results$Cancer_Type, levels=c("Pancancer", sig_mutations$Cancer_Type))
  
  plot_list[[g]] <- ggplot(mutation_cancer_type_pancancer_KO_results, aes(x=Cancer_Type, y= Gene_KO, fill=OXPHOS_state)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + labs(title = paste(g, "KO effect")) + xlab("") + ylab("Gene dependency") +
    theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(),
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=12),
                            axis.line=element_line(linewidth=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",linewidth=0.5),
                            axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12)) + ylim(0,1.2) +
    geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 1, alpha = 0.6) + scale_fill_manual(values=c("skyblue1", "tomato")) + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
              aes(label = paste("Diff.  = ", round(Median_sub, digits = 2)), x = Cancer_Type, y = 1.15), size = 3, color = "red") +
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
              aes(label = paste("p = ", format(adj_pvalue, scientific = TRUE, digits = 3)), x = Cancer_Type, y = 1.2), size = 3, color = "red") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("low = ", num_OXPHOS_low_cell_lines), x = Cancer_Type, y = 1.08), size = 3, color = "black") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("high = ", num_OXPHOS_high_cell_lines), x = Cancer_Type, y = 1.12), size = 3, color = "black")
  
  
}
)

p <- ggarrange(plotlist = plot_list)

ggsave(file.path(outDir, 'ETC_pancancer_PTEN_boxplots.pdf'), plot = p,device = 'pdf',width =15, height = 12,dpi=300)



#boxplots of OXPHOS low vulnerabilities
OXPHOS_low_vulnerabilites <- c("EFR3A","RAB6A", "SNAP23")

plot_list <- list()
plot_list <- lapply(OXPHOS_low_vulnerabilites, function(g) { 
  
  #getting pancancer difference in median results
  pancancer_median_diff <- pancancer_sig_vulnerabilities[pancancer_sig_vulnerabilities$Gene == g, ]$Median_sub
  
  #getting mutations and cancer types' differences that are greater than or equal to pancancer difference
  sig_cancer_types <-cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$Gene == g &
                                                                    cancer_type_pancancer_permutation_results.df$adj_pvalue <= 0.01 &
                                                                    cancer_type_pancancer_permutation_results.df$Median_sub <= pancancer_median_diff,c(1,3,4:6)]
  
  sig_mutations <- permutation_each_mutation_results.df[permutation_each_mutation_results.df$Gene_KO == g &
                                                          permutation_each_mutation_results.df$adj_pvalue <= 0.01 &
                                                          permutation_each_mutation_results.df$Median_sub <= pancancer_median_diff,-c(2)]
  
  colnames(sig_mutations) <- colnames(sig_cancer_types)
  
  #combining significant cancer types and mutations for labels
  sig_mutations_cancer_types <- rbind(sig_cancer_types, sig_mutations)
  
  
  #extracting pancancer cell lines
  pancancer_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("Cancer_Type", "OXPHOS_state", g)]
  pancancer_KO_results$Cancer_Type <- c("Pancancer")
  
  #extracting individual cancer type cell lines
  cancer_type_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0 &
                                                                           mutation_OXPHOS_variable_cell_lines_final.df$Cancer_Type %in% sig_cancer_types$Cancer_Type[-(3)],c("Cancer_Type","OXPHOS_state", g)]
  
  #combining pancancer and cancer type KO results
  cancer_type_pancancer_KO_results <- rbind(cancer_type_KO_results, pancancer_KO_results)
  colnames(cancer_type_pancancer_KO_results)[3] <- c("Gene_KO")
  
  
  #extracting mutated cell lines
  mutation_KO_results <- mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0,c("StrippedCellLineName","OXPHOS_state", g, sig_mutations$Cancer_Type)]
  
  mutation_KO_results_for_plotting <- data.frame()
  for (i in 1:length(sig_mutations$Cancer_Type)) {
    
    mutation <- sig_mutations$Cancer_Type[i]
    df <- mutation_KO_results[, c(1:3, 3+i)]
    colnames(df)[4] <- "Mutation"
    mutated_in_gene <- na.omit(df[df$Mutation == 1, ])
    mutated_in_gene$Mutation_status <- mutation
    
    mutation_KO_results_for_plotting <- rbind(mutation_KO_results_for_plotting,mutated_in_gene[,-c(1,4)])
    
  }
  
  #rearranging dataframe so we can bind it with cancer/pancancer results
  mutation_KO_results_for_plotting <- mutation_KO_results_for_plotting[,c(3,1,2)]
  colnames(mutation_KO_results_for_plotting) <- colnames(cancer_type_pancancer_KO_results)
  
  #combining mutation, cancer type,and pancancer results
  mutation_cancer_type_pancancer_KO_results <- rbind(cancer_type_pancancer_KO_results,mutation_KO_results_for_plotting)
  mutation_cancer_type_pancancer_KO_results <- merge(mutation_cancer_type_pancancer_KO_results,sig_mutations_cancer_types[,c(1:5)], by = "Cancer_Type")
  
  #converting table from wide to long for plotting purposes, and using boxplots to plot Gene_dependency scores in OXPHOS high vs. OXPHOS low cell lines
  mutation_cancer_type_pancancer_KO_results$OXPHOS_state <- factor(mutation_cancer_type_pancancer_KO_results$OXPHOS_state, levels=c("OXPHOS_low", "OXPHOS_high"))
  mutation_cancer_type_pancancer_KO_results$Cancer_Type <- factor(mutation_cancer_type_pancancer_KO_results$Cancer_Type, levels=c("Pancancer", sig_mutations$Cancer_Type, sig_cancer_types[sig_cancer_types$Cancer_Type != "Pancancer", ]$Cancer_Type))
  
  
  #generating plots
  plot_list[[g]] <- ggplot(mutation_cancer_type_pancancer_KO_results, aes(x=Cancer_Type, y= Gene_KO, fill=OXPHOS_state)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + labs(title = paste(g, "KO effects across cancers and mutations")) + xlab("") + ylab("Gene dependency") +
    theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(),
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=12),
                            axis.line=element_line(linewidth=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",linewidth=0.5),
                            axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12)) + ylim(0,1.2) +
    geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 1, alpha = 0.6) + scale_fill_manual(values=c("skyblue1", "tomato")) + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
              aes(label = paste("Diff.  = ", round(Median_sub, digits = 2)), x = Cancer_Type, y = 1.15), size = 3, color = "red") +
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,4,5)]),
              aes(label = paste("p = ", format(adj_pvalue, scientific = TRUE, digits = 3)), x = Cancer_Type, y = 1.2), size = 3, color = "red") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("low = ", num_OXPHOS_low_cell_lines), x = Cancer_Type, y = 1.08), size = 3, color = "black") + 
    geom_text(data =unique(mutation_cancer_type_pancancer_KO_results[mutation_cancer_type_pancancer_KO_results$OXPHOS_state == "OXPHOS_low",c(1,2,6,7)]),
              aes(label = paste("high = ", num_OXPHOS_high_cell_lines), x = Cancer_Type, y = 1.12), size = 3, color = "black")
}
)

p <- ggarrange(plotlist = plot_list, ncol = 1)

ggsave(file.path(outDir, 'OXPHOS_low_vulnerabilities_boxplots.pdf'), plot = p,device = 'pdf',width =15, height = 25,dpi=300)





