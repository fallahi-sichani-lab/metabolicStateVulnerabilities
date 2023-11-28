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
library(pheatmap)
library(tidyr)
library(devtools)
library(ggplot2)
library(gridExtra)
library(glmnet)
library(fgsea)
library(ggrepel)
library(randomcoloR)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis")

#loading data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/organized_CRISPR_gene_effect.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/significant_OXPHOS_pan-cancer_targets.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/driver_mutations_OXPHOS_variable_cell_lines.RData")

#mutations for permutation analysis
mutation_candidates <- names(mutation_OXPHOS_variable_cell_lines_final.df)[32:42]

#calculating difference in medians between OXPHOS high and OXPHOS low cell lines for each mutation candidate
gene_knockout_mutation <- data.frame()
for (m in mutation_candidates) {
  
  #filtering out cell lines with a given mutation
  cell_lines_mutation <- mutation_OXPHOS_variable_cell_lines_final.df[,c("OXPHOS_state",as.character(unique(sig_candidates$Gene)),m)]
  cell_lines_mutation <- cell_lines_mutation[!cell_lines_mutation$OXPHOS_state == 0 & cell_lines_mutation[,m] == 1,]
  cell_lines_mutation <- na.omit(cell_lines_mutation)
  
  
  median_dependency <- cell_lines_mutation %>% dplyr::group_by(OXPHOS_state) %>% dplyr::summarise(dplyr::across(colnames(cell_lines_mutation)[2]:colnames(cell_lines_mutation)[28], median))
  
  median_diff_OXPHOS_high_low <- median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high",c(2:ncol(median_dependency))]-median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low",c(2:ncol(median_dependency))]
  
  gene_knockouts_across_mutations <- data.frame(Mutation = m, Gene_KO = colnames(median_diff_OXPHOS_high_low), Median_difference = t(median_diff_OXPHOS_high_low), OXPHOS_high_num_cell_lines = sum(cell_lines_mutation$OXPHOS_state == "OXPHOS_high"), 
                                                OXPHOS_low_num_cell_lines = sum(cell_lines_mutation$OXPHOS_state == "OXPHOS_low"), possible_unique_combinations = choose(nrow(cell_lines_mutation), sum(cell_lines_mutation$OXPHOS_state == "OXPHOS_high")))
  
  gene_knockout_mutation <- rbind(gene_knockout_mutation, gene_knockouts_across_mutations)
    

}


#permutation test to determine if knockout effect across each individual mutation type
#Note: code is computationally expensive, make sure to save results at the end of the code for analysis/don't need rerun permutation analysis 
permutation_iterations_mutation_type <- 1:4000
sig_OXPHOS_pancancer_targets <- as.character(unique(sig_candidates$Gene))
  
permutation_each_mutation_results.df <-data.frame()
for (m in mutation_candidates) {
  
  #filtering out cell lines with a given mutation
  cell_lines_mutation <- mutation_OXPHOS_variable_cell_lines_final.df[,c("OXPHOS_state",as.character(unique(sig_candidates$Gene)),m)]
  cell_lines_mutation <- cell_lines_mutation[!cell_lines_mutation$OXPHOS_state == 0 & cell_lines_mutation[,m] == 1,]
  cell_lines_mutation <- na.omit(cell_lines_mutation)
  
  #extracting mutation type and correct OXPHOS state labels
  OXPHOS_state_cell_lines_mutation_type <- cell_lines_mutation$OXPHOS_state
  
  #randomly shuffling OXPHOS labels for the correct number of iterations
  set.seed(123)
  shuffle_OXPHOS_state_list_mutation_type <- lapply(permutation_iterations_mutation_type,function(x) sample(length(OXPHOS_state_cell_lines_mutation_type)))
  
  #ensuring each iteration of the random shuffle are unique/not overlapping
  iterations.df <- data.frame(t(data.frame(shuffle_OXPHOS_state_list_mutation_type[permutation_iterations_mutation_type])))
  unique_iterations.df <- iterations.df %>% distinct()
    
    # Initialize an empty list to store the data frames
    permutation_each_gene_list <- list()
    
    # Iterate over the indices of unique iterations
    permutation_each_gene_list <- parallel::mclapply(1:nrow(unique_iterations.df), function(i) {
      iteration_numbers <- as.numeric(unique_iterations.df[i,])
      permutated_OXPHOS_states <- cell_lines_mutation[iteration_numbers, ]$OXPHOS_state
      permutated_df <- cell_lines_mutation
      permutated_df$OXPHOS_state <- permutated_OXPHOS_states
      
      df <- permutated_df[, c(sig_OXPHOS_pancancer_targets, "OXPHOS_state")]
      median_dependency <- df %>% dplyr::group_by(OXPHOS_state) %>% dplyr::summarise(dplyr::across(colnames(cell_lines_mutation)[2]:colnames(cell_lines_mutation)[28], median))
      median_diff_OXPHOS_high_low <- as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high", sig_OXPHOS_pancancer_targets]) - as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low", sig_OXPHOS_pancancer_targets])
      
      permutation_each_gene <- data.frame(Mutation = m, Gene_KO = colnames(median_dependency)[-c(1)], Iteration = i, Median_difference = median_diff_OXPHOS_high_low)
      
      # Return the data frame for each iteration
      return(permutation_each_gene)
    }, mc.cores = 16) #defining number of cores to use for parallel processing
    
    # Combine the list of data frames into a single data frame
    permutation_each_gene.df <- do.call(rbind, permutation_each_gene_list)
    
    
    
    #determining if KO effect is significant in a given mutation
    # Initialize an empty list to store the data frames
    permutation_mutation_type_list <- list()
    
    
    # Calculate the number of OXPHOS high and low cell lines
    num_OXPHOS_high <- as.numeric(table(OXPHOS_state_cell_lines_mutation_type)[1])
    num_OXPHOS_low <- as.numeric(table(OXPHOS_state_cell_lines_mutation_type)[2])
    
    # Iterate over significant pan-cancer OXPHOS gene knockout candidates
    permutation_mutation_type_list <- parallel::mclapply(sig_OXPHOS_pancancer_targets, function(v) {
      # Get the true difference in medians for the current mutation
      median_difference <- gene_knockout_mutation[gene_knockout_mutation$Mutation == m & gene_knockout_mutation$Gene_KO == v, ]$Median_difference
      
      # Calculate the adjusted p-value for each gene knockout
      adj_pval <- ifelse(median_difference > 0,
                         (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference > median_difference) + 1) / (max(permutation_iterations_mutation_type) + 1),
                         (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference < median_difference) + 1) / (max(permutation_iterations_mutation_type) + 1))
      
      # Create the permutation_mutation_type data frame
      permutation_mutation_type <- data.frame(
        Mutation = m,
        Gene_KO = v,
        adj_pvalue = adj_pval,
        Median_sub = median_difference,
        num_OXPHOS_high_cell_lines = num_OXPHOS_high,
        num_OXPHOS_low_cell_lines = num_OXPHOS_low
      )
      
      # Return the data frame for each gene
      return(permutation_mutation_type)
    }, mc.cores = 16)
    
    # Combine the list of data frames into a single data frame
    permutation_mutation_type.df <- do.call(rbind, permutation_mutation_type_list)
    
    
    #combining all dataframes across all mutations
    permutation_each_mutation_results.df <- rbind(permutation_each_mutation_results.df, permutation_mutation_type.df)
  
}


#saving permutation results for further analysis (DO NOT NEED TO RERUN CODE ABOVE)
save(permutation_each_mutation_results.df, file = "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/permutation_mutation_analysis_results.RData")

###
###

#loading permutation results 
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/permutation_mutation_analysis_results.RData")


#scatter plot of significant knockout effects across mutations

#creating color palette for  mutations with significant knockout effects
new_palette <- c("darkorange", "goldenrod1", "lightblue", "darkred", "yellow", "red", "turquoise", "navy", "steelblue2","dodgerblue", "turquoise4")

mat_colors <- list(Mutation = new_palette)
names(mat_colors$Mutation) <- sort(unique(permutation_each_mutation_results.df[permutation_each_mutation_results.df$adj_pvalue <= 0.05 & abs(permutation_each_mutation_results.df$Median_sub) >= 0.15, ]$Mutation))

mutation_pallette <- list()
for (i in sort(unique(permutation_each_mutation_results.df[permutation_each_mutation_results.df$adj_pvalue <= 0.05 & abs(permutation_each_mutation_results.df$Median_sub) >= 0.15, ]$Mutation))) {
  mutation_pallette[[i]] <- mat_colors[["Mutation"]][[i]]
}

#adding Pancancer to cancer palette
mutation_pallette[["Pancancer"]] <- "black"

#merging mutation and pancancer results
sig_OXPHOS_pancancer_targets <- sig_candidates[sig_candidates$Cancer_Type == "Pancancer",c(1:4)]
colnames(sig_OXPHOS_pancancer_targets) <- c("Mutation", "Gene_KO", "adj_pvalue", "Median_sub")

mutation_pancancer_KO.df <- rbind(sig_OXPHOS_pancancer_targets, permutation_each_mutation_results.df[,c(1:4)])

#ordering candidates from descending order
pan_cancer_gene_order <- arrange(sig_OXPHOS_pancancer_targets, Median_sub)
mutation_pancancer_KO.df$Gene_KO <- factor(mutation_pancancer_KO.df$Gene_KO, levels = pan_cancer_gene_order$Gene)



#subplot scatter plots of mutation knockout effects across different mutations
plot_list <- list()
for (g in unique(mutation_pancancer_KO.df$Mutation)[-1]) {

  df <- mutation_pancancer_KO.df[mutation_pancancer_KO.df$adj_pvalue <= 0.05 & abs(mutation_pancancer_KO.df$Median_sub) >= 0.15 & mutation_pancancer_KO.df$Mutation %in% c("Pancancer", g), ]
  
  #getting cell numbers for label/title of each plot
  num_cell_lines_OXPHOS_high <- unique(permutation_each_mutation_results.df[permutation_each_mutation_results.df$Mutation == g, c("num_OXPHOS_high_cell_lines")])
  num_cell_lines_OXPHOS_low <- unique(permutation_each_mutation_results.df[permutation_each_mutation_results.df$Mutation == g, c("num_OXPHOS_low_cell_lines")])
  
  title1 <- paste("Significant vulnerabilites in", g, "mutant cell lines")
  title2 <- paste("OXPHOS high =", num_cell_lines_OXPHOS_high)
  title3 <-  paste("OXPHOS low =", num_cell_lines_OXPHOS_low)
  title_final <- paste(title1, title2,  title3, sep = "\n")
  

  plot_list[[g]] <- ggplot(df, aes(x=Median_sub, y=Gene_KO, fill=Mutation)) + geom_point(color ="black",shape = 21, size = 4) + theme_bw() + labs(x = "Gene_dependency scores' difference (OXPHOS high - OXPHOS low)", y = "", title = title_final) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position="none", 
          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
          axis.line=element_line(linewidth=0.5,colour="black"),
          axis.ticks = element_line(colour = "black",linewidth=0.5),
          axis.text.x=element_text(colour="black", size = 8),
          axis.text.y=element_text(colour="black", size = 8),
          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
          strip.text=element_text(size=8), text = element_text(size =8))  + geom_vline(xintercept = 0, colour="red", linetype = "dashed", linewidth = 0.8) + scale_fill_manual(values =mutation_pallette)
}

p <- grid.arrange(grobs=plot_list,ncol=4)
ggsave(file.path(outDir, 'mutation_knockout_scatter_plot_individual_mutations_subplots.pdf'), plot = p,device = 'pdf',width =20, height = 22,dpi=300)




#scatter plot of all mutations/significant vulnerabilities on one plot
p <- ggplot(mutation_pancancer_KO.df[mutation_pancancer_KO.df$adj_pvalue <= 0.05 & abs(mutation_pancancer_KO.df$Median_sub) >= 0.15,], aes(x=Median_sub, y=Gene_KO, fill=Mutation)) + geom_point(color ="black",shape = 21, size = 4) + theme_bw() + labs(x = "Gene_dependency scores' difference (OXPHOS high - OXPHOS low)", y = "", title = "Significant OXPHOS vulnerabilites across mutations") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="right", 
        panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=10),
        axis.line=element_line(linewidth=0.5,colour="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.5),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
        strip.text=element_text(size=8), text = element_text(size =10))  + geom_vline(xintercept = 0, colour="red", linetype = "dashed", linewidth = 0.8) + scale_fill_manual(values =mutation_pallette)

ggsave(file.path(outDir, 'mutation_knockout_scatter_plot_all_mutations.pdf'), plot = p,device = 'pdf',width = 8, height = 10,dpi=300)







