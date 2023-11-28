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
library(ggvenn)
library(combinat)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS")

#loading Rdata files
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/Elastic_net_OXPHOS/cell_line_to_cell_line_metabolic_variability_scores.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/organized_CRISPR_gene_effect.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")

#loading PLSR OXPHOS variables 
PLSR_OXPHOS_variables.df <- read.table("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/OXPHOS_PLSR_selected_variables.txt", header = TRUE)

#selecting PLSR variables with VIP scores >= |1|
PLSR_OXPHOS_variables.df <- PLSR_OXPHOS_variables.df[abs(PLSR_OXPHOS_variables.df$VIP_signed) >= 1,]

#cleaning up CRISPR dataset
cell_lines_CRISPR <- data.frame(rownames(CRISPR_CCLE_data),CRISPR_CCLE_data)
rownames(cell_lines_CRISPR) <- NULL
colnames(cell_lines_CRISPR)[1] <- "cell_line"

#exploring missing data in CRISPR dataset and removing any cell lines/genes with no data
genes_NA <- as.data.frame(colSums(is.na(cell_lines_CRISPR[,c(2:17932)]))) 
genes_NA$empty <- 0
genes_NA <-genes_NA[genes_NA$`colSums(is.na(cell_lines_CRISPR[, c(2:17932)]))` == 0,]
genes_CRISPR_removed_NA<- cell_lines_CRISPR[,c("cell_line", rownames(genes_NA))]


#filtering CRISPR data to only selected pathways' cell lines
metabolic_pathways <- "Oxidative phosphorylation"

Gene_dependency_scores_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  metabolic_pathway_cell_lines <- cell_line_metabolic_variability.df[cell_line_metabolic_variability.df$Pathway == p,]
  metabolic_score_Gene_dependency <- merge(metabolic_pathway_cell_lines,genes_CRISPR_removed_NA, by = "cell_line") 
  
  Gene_dependency_scores_metabolic_pathways <- rbind(Gene_dependency_scores_metabolic_pathways, metabolic_score_Gene_dependency)
}


#correlation analysis of pan-cancer targets
Gene_dependency_scores_OXPHOS_PLSR_variables_correlation_analysis <- Gene_dependency_scores_metabolic_pathways[,c("cell_line", "mean_zscore", "Cancer_Type", PLSR_OXPHOS_variables.df$OXPHOS_genes)]

#correlation between OXPHOS state score and gene knockouts at pan-cancer level
OXPHOS_PLSR_correlation.test <- lapply(Gene_dependency_scores_OXPHOS_PLSR_variables_correlation_analysis[,-c(1,3)], function(x) cor.test(Gene_dependency_scores_OXPHOS_PLSR_variables_correlation_analysis$mean_zscore, x,  method = "pearson"))


#extracting correlation coefficient and significant from pan-cancer analysis
OXPHOS_PLSR_variables <- PLSR_OXPHOS_variables.df$OXPHOS_genes
correlation.test_pancancer.df <- data.frame()
for (v in OXPHOS_PLSR_variables) {
  correlation.test_pancancer <- data.frame(Cancer_Type = "Pancancer", Gene = v, Pearson_R = OXPHOS_PLSR_correlation.test[[v]][["estimate"]][["cor"]], p_value = OXPHOS_PLSR_correlation.test[[v]][["p.value"]], num_cell_lines = nrow(Gene_dependency_scores_OXPHOS_PLSR_variables_correlation_analysis))
  
  correlation.test_pancancer.df <- rbind(correlation.test_pancancer.df, correlation.test_pancancer)
}

#FDR adjustment for correlation analysis
correlation.test_pancancer.df$FDR <- p.adjust(correlation.test_pancancer.df$p_value, method = c("fdr"))


####
####


#statistical comparison of OXPHOS high versus OXPHOS low cell lines with respect to knockout effects using permutation test
#filtering datasets to remove samples with metabolic state scores close to 0 based on 33% and 67% distribution thresholds
Gene_dependency_scores_metabolic_pathways_filtered <- data.frame()
for (p in metabolic_pathways) {
  df <- Gene_dependency_scores_metabolic_pathways[Gene_dependency_scores_metabolic_pathways$Pathway == p, ]
  df_quants <- as.data.frame(t(quantile(df$mean_zscore, probs = c(0.33,0.67))))
  
  df_filtered <- df[df$mean_zscore < df_quants$`33%` | df$mean_zscore > df_quants$`67%`,]
  
  
  Gene_dependency_scores_metabolic_pathways_filtered <- rbind(Gene_dependency_scores_metabolic_pathways_filtered, df_filtered)
}



#filtering Gene_dependency dataset to PLSR selected gene knockouts
Gene_dependency_scores_OXPHOS_PLSR_variables <- Gene_dependency_scores_metabolic_pathways_filtered[,c("cell_line", "mean_zscore", "Cancer_Type",PLSR_OXPHOS_variables.df$OXPHOS_genes)]

Gene_dependency_scores_OXPHOS_PLSR_variables$OXPHOS_state<-ifelse(Gene_dependency_scores_OXPHOS_PLSR_variables$mean_zscore>0,"OXPHOS_high","OXPHOS_low")


#calculating median differences in knockout scores (OXPHOS high cell lines - OXPHOS low cell lines) for pan-cancer
OXPHOS_PLSR_variables <- PLSR_OXPHOS_variables.df$OXPHOS_genes
df <- Gene_dependency_scores_OXPHOS_PLSR_variables[,-c(1,2,3)]
median_dependency <- df %>% group_by(OXPHOS_state) %>% summarise(across(colnames(df)[1]:colnames(df)[ncol(df)-1], median))

median_difference_pancancer.df <- data.frame()
for (v in OXPHOS_PLSR_variables) {
  median_diff_OXPHOS_high_low <- as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high",v]) - as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low",v])
  gene_pancancer_median_difference <- data.frame(Gene = v, Cancer_Type = "Pancancer", Median_sub = median_diff_OXPHOS_high_low)
  gene_pancancer_median_difference$OXPHOS_high_count <- nrow(Gene_dependency_scores_OXPHOS_PLSR_variables)/2
  gene_pancancer_median_difference$OXPHOS_low_count <- nrow(Gene_dependency_scores_OXPHOS_PLSR_variables)/2
  gene_pancancer_median_difference$possible_unique_combinations <- choose(nrow(Gene_dependency_scores_OXPHOS_PLSR_variables), (nrow(Gene_dependency_scores_OXPHOS_PLSR_variables)/2))
  gene_pancancer_median_difference$lowest_pvalue <- 1/gene_pancancer_median_difference$possible_unique_combinations
  
  median_difference_pancancer.df <- rbind(median_difference_pancancer.df, gene_pancancer_median_difference)
}


#permutation test to determine if knockout effect is truly significant at pan-cancer level
#generating 8000 random sampling of pan-cancer OXPHOS high and low cell lines
permutation_iterations <- 1:8000
OXPHOS_state_cell_lines <- Gene_dependency_scores_OXPHOS_PLSR_variables$OXPHOS_state

set.seed(123)
shuffle_OXPHOS_state_list <- lapply(permutation_iterations,function(x) sample(length(OXPHOS_state_cell_lines))) 

iterations.df <- data.frame(t(data.frame(shuffle_OXPHOS_state_list[permutation_iterations ])))
unique_iterations.df <- iterations.df %>% distinct()
  
# Initialize an empty list to store the data frames
permutation_each_gene_list <- list()
  
# Iterate over the indices of unique iterations
permutation_each_gene_list <- parallel::mclapply(1:nrow(unique_iterations.df), function(i) {
        #running permutation for each gene
        iteration_numbers <- as.numeric(unique_iterations.df[i,])
        permutated_OXPHOS_states <- Gene_dependency_scores_OXPHOS_PLSR_variables[iteration_numbers,]$OXPHOS_state
        permutated_df <- Gene_dependency_scores_OXPHOS_PLSR_variables
        permutated_df$OXPHOS_state <- permutated_OXPHOS_states
        
        df <- permutated_df[,c(OXPHOS_PLSR_variables, "OXPHOS_state")]
        median_dependency <- df %>% group_by(OXPHOS_state) %>% summarise(across(colnames(df)[1]:colnames(df)[ncol(df)-1], median))
        median_diff_OXPHOS_high_low <- as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high",OXPHOS_PLSR_variables]) - as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low",OXPHOS_PLSR_variables])
        
        permutation_each_gene <- data.frame(Gene = OXPHOS_PLSR_variables, Iteration = i, Median_difference = median_diff_OXPHOS_high_low)
        
        # Return the data frame for each iteration
        return(permutation_each_gene)
}, mc.cores = 16) #defining number of cores to use for parallel processing
  
# Combine the list of data frames into a single data frame
permutation_each_gene.df <- do.call(rbind, permutation_each_gene_list)
  
  
  #calculating adjusted p-value for each gene knockout in pan-cancer
  permutation_pancancer.df <- data.frame()
  for (v in OXPHOS_PLSR_variables) {
    
    median_difference <- median_difference_pancancer.df[median_difference_pancancer.df$Gene == v,]$Median_sub
    
    adj_pval <- ifelse(median_difference > 0,
                       (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference > median_difference) + 1) / (max(permutation_iterations) + 1),
                       (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference < median_difference) + 1) / (max(permutation_iterations) + 1))
    
    
    permutation_pancancer <- data.frame(Cancer_Type = "Pancancer", Gene = v, adj_pvalue =adj_pval, Median_sub = median_difference_pancancer.df[median_difference_pancancer.df$Gene == v,]$Median_sub, 
                                        num_OXPHOS_high_cell_lines = nrow(df[df$OXPHOS_state == "OXPHOS_high",]), num_OXPHOS_low_cell_lines = nrow(df[df$OXPHOS_state == "OXPHOS_low",]), num_unique_iterations = nrow(unique_iterations.df),
                                        lowest_pvalue = 1/nrow(unique_iterations.df))
    
    #final dataframe
    permutation_pancancer.df <- rbind(permutation_pancancer.df, permutation_pancancer)
  }
  

  ####
  ####  
  
  
#comparing significant candidates between permutation and correlation analysis
significant_correlation_pancancer_targets <- correlation.test_pancancer.df[correlation.test_pancancer.df$FDR <= 0.05,]
significant_permutation_pancancer_targets <- permutation_pancancer.df[permutation_pancancer.df$adj_pvalue <= 0.05,]

#venn diagram comparing significant pancancer_targets from the two statisitcal methods
x <- list(
  Correlation_Pancancer_Vulnerabilites = significant_correlation_pancancer_targets$Gene, 
  Permutation_Pancancer_Vulnerabilites = significant_permutation_pancancer_targets$Gene 
)

p <- ggvenn(x, text_size = 8)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/venn_diagram_permutation_correlation_anlayses_comparison.pdf', plot = p,device = 'pdf',width =10, height = 8,dpi=300)


#using only significant pan-cancer targets in which permutated adjusted p-value and correlation p-value are <= 0.05 and median difference between OXPHOS high and OXPHOS low cell lines is at least 0.15
sig_OXPHOS_pancancer_targets <- permutation_pancancer.df[permutation_pancancer.df$adj_pvalue<= 0.05 & abs(permutation_pancancer.df$Median_sub) >= 0.15, ]

####
####

#permutation test for each individual cancer type
#calculating median difference for each cancer type that are variable in OXPHOS
median_difference_cancer_types.df <- data.frame()
cancer_types_OXPHOS <- unique(Gene_dependency_scores_OXPHOS_PLSR_variables$Cancer_Type)

for (c in cancer_types_OXPHOS[-c(10,11)]) { #removing Cervical Adenocarcinoma and Breast Ductal Carcinoma In Situ from analysis due to small number of cell lines
  df <- Gene_dependency_scores_OXPHOS_PLSR_variables[Gene_dependency_scores_OXPHOS_PLSR_variables$Cancer_Type == c,]
  
  median_dependency <- df[,-c(1,2,3)] %>% group_by(OXPHOS_state) %>% summarise(across(colnames(df)[4]:colnames(df)[ncol(df)-1], median))
  
  
    for (v in sig_OXPHOS_pancancer_targets$Gene) {
        median_diff_OXPHOS_high_low <- as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high",v]) - as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low",v])
        
        possible_unique_combinations_per_cancer_type <- choose(nrow(df), nrow(df[df$OXPHOS_state == "OXPHOS_high",]))
    
        median_difference_cancer_types <- data.frame(Gene = v, Cancer_Type = c, Median_sub = median_diff_OXPHOS_high_low, OXPHOS_high_count = nrow(df[df$OXPHOS_state == "OXPHOS_high",]), 
                                                     OXPHOS_low_count = nrow(df[df$OXPHOS_state == "OXPHOS_low",]), possible_unique_combinations = possible_unique_combinations_per_cancer_type, lowest_pvalue = 1/possible_unique_combinations_per_cancer_type)

        median_difference_cancer_types.df <- rbind(median_difference_cancer_types.df, median_difference_cancer_types)

           }
  }


#permutation test to determine if knockout effect is truly significant
#Note: code is computationally expensive, make sure to save results at the end of the code for analysis/don't need rerun permutation analysis 
permutation_iterations_cancer_type <- 1:8000

permutation_cancer_type_results.df <- data.frame()
for (c in unique(median_difference_cancer_types.df$Cancer_Type)) {
  
  #extracting cancer type and correct OXPHOS state labels
  cancer_type.df <- Gene_dependency_scores_OXPHOS_PLSR_variables[Gene_dependency_scores_OXPHOS_PLSR_variables$Cancer_Type == c,]
  OXPHOS_state_cell_lines_cancer_type <- cancer_type.df$OXPHOS_state
  
  #randomly shuffling OXPHOS labels for the correct number of iterations
  set.seed(123)
  shuffle_OXPHOS_state_list_cancer_type <- lapply(permutation_iterations_cancer_type,function(x) sample(length(OXPHOS_state_cell_lines_cancer_type)))
  iterations.df <- data.frame(t(data.frame(shuffle_OXPHOS_state_list_cancer_type[permutation_iterations ])))
  unique_iterations.df <- iterations.df %>% distinct()
    
    
    # Initialize an empty list to store the data frames
    permutation_each_gene_list <- list()
    
    # Iterate over the indices of unique iterations
    permutation_each_gene_list <- parallel::mclapply(1:nrow(unique_iterations.df), function(i) {
      iteration_numbers <- as.numeric(unique_iterations.df[i,])
      permutated_OXPHOS_states <- cancer_type.df[iteration_numbers,]$OXPHOS_state
      permutated_df <- cancer_type.df
      permutated_df$OXPHOS_state <- permutated_OXPHOS_states
      
      df <- permutated_df[,c(sig_OXPHOS_pancancer_targets$Gene, "OXPHOS_state")]
      median_dependency <- df %>% group_by(OXPHOS_state) %>% summarise(across(colnames(df)[1]:colnames(df)[ncol(df)-1], median))
      median_diff_OXPHOS_high_low <- as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_high",sig_OXPHOS_pancancer_targets$Gene]) - as.numeric(median_dependency[median_dependency$OXPHOS_state == "OXPHOS_low",sig_OXPHOS_pancancer_targets$Gene])
      
      permutation_each_gene <- data.frame(Cancer_Type = c, Gene = sig_OXPHOS_pancancer_targets$Gene, Iteration = i, Median_difference = median_diff_OXPHOS_high_low)
      
      # Return the data frame for each iteration
      return(permutation_each_gene)
    }, mc.cores = 16) #defining number of cores to use for parallel processing
    
    # Combine the list of data frames into a single data frame
    permutation_each_gene.df <- do.call(rbind, permutation_each_gene_list)
    
    
  
    # Initialize an empty list to store the data frames
    permutation_cancer_type_list <- list()
    
    
    # Calculate the number of OXPHOS high and low cell lines
    num_OXPHOS_high <- as.numeric(table(OXPHOS_state_cell_lines_cancer_type)[1])
    num_OXPHOS_low <- as.numeric(table(OXPHOS_state_cell_lines_cancer_type)[2])
    
    # Iterate over significant Opan-cancer OXPHOS gene knockout candidates
    permutation_cancer_type_list <- parallel::mclapply(sig_OXPHOS_pancancer_targets$Gene, function(v) {
      
      # Get the median difference for the current gene
      median_difference <- median_difference_cancer_types.df[median_difference_cancer_types.df$Cancer_Type == c & median_difference_cancer_types.df$Gene == v,]$Median_sub
      
      # Calculate the adjusted p-value for each gene knockout
      adj_pval <- ifelse(median_difference > 0,
                         (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference > median_difference) + 1) / (max(permutation_iterations_cancer_type) + 1),
                         (sum(permutation_each_gene.df[permutation_each_gene.df$Gene == v, ]$Median_difference < median_difference) + 1) / (max(permutation_iterations_cancer_type) + 1))
      
      # Create the permutation_mutation_type data frame
      permutation_cancer_type <- data.frame(
        Cancer_Type = c,
        Gene = v,
        adj_pvalue = adj_pval,
        Median_sub = median_difference,
        num_OXPHOS_high_cell_lines = num_OXPHOS_high,
        num_OXPHOS_low_cell_lines = num_OXPHOS_low,
        num_unique_iterations = nrow(unique_iterations.df)
      )
      
      # Return the data frame for each gene
      return(permutation_cancer_type)
    }, mc.cores = 16)
    
    # Combine the list of data frames into a single data frame
    permutation_cancer_type.df <- do.call(rbind, permutation_cancer_type_list)
    
    #combining all dataframes across all mutations
    permutation_cancer_type_results.df <- rbind(permutation_cancer_type_results.df, permutation_cancer_type.df)
    
    
  }


#merging permutation tests together from pan-cancer and cancer type analyses
cancer_type_pancancer_permutation_results.df <- rbind(permutation_cancer_type_results.df,permutation_pancancer.df[,-c(8)])

#saving permuation results from pan-cancer and cancer type analysis (DO NOT NEED TO RERUN PERMUTATION CODE ABOVE)
save(cancer_type_pancancer_permutation_results.df, file = "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/pancancer_cancer_type_permutation_results.RData")

####
####


#load permutation results without running code above
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/pancancer_cancer_type_permutation_results.RData")


#scatter plot of significant OXPHOS PLSR candidates across Pancancer and OXPHOS cancer types, selecting significant candidates that have adjusted p-value from permutation test <= 0.05 and median difference of at least |0.15|
sig_candidates <- cancer_type_pancancer_permutation_results.df[cancer_type_pancancer_permutation_results.df$adj_pvalue <= 0.05 & abs(cancer_type_pancancer_permutation_results.df$Median_sub) >= 0.15, ]

#ordering candidates from descending order
pan_cancer_gene_order <- arrange(median_difference_pancancer.df, Median_sub)

sig_candidates$Gene <- factor(sig_candidates$Gene, levels = pan_cancer_gene_order$Gene)

#creating color palette for cancers 
new_palette <- c("darkorange", "red", "blue", "hotpink", "darkorchid4", "green", "goldenrod1","darkgray", "brown4", "lightblue","lightpink","cyan", "yellow", "mediumpurple1", "darkgreen", "lightsalmon3")

mat_colors <- list(cancer_type = new_palette)
names(mat_colors$cancer_type) <- sort(unique(unique(median_difference_cancer_types.df$Cancer_Type)))

cancer_pallette <- list()
for (i in sort(unique(unique(median_difference_cancer_types.df$Cancer_Type)))) {
  cancer_pallette[[i]] <- mat_colors[["cancer_type"]][[i]]
}

#adding Pancancer to cancer palette
cancer_pallette[["Pancancer"]] <- "black"



#plotting significant PLSR candidates across all cancer types
p <- ggplot(data=sig_candidates, aes(x=Median_sub, y=Gene, fill=Cancer_Type)) + geom_point(color ="black",shape = 21, size = 4) + theme_bw() + labs(x = "Gene_dependency scores' differences (OXPHOS high cell lines - OXPHOS low cell lines)", y = "", title = "OXPHOS Elastic-Net/PLSR significant vulnerabilites across cancers") +
 theme(panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          legend.position="right", 
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=10),
                          axis.line=element_line(linewidth=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x=element_text(colour="black", size = 10),
                          axis.text.y=element_text(colour="black", size = 10),
                          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                          strip.text=element_text(size=8), text = element_text(size =10))  + geom_vline(xintercept = 0, colour="red", linetype = "dashed", linewidth = 0.8) + scale_fill_manual(values = cancer_pallette)

ggsave(file.path(outDir, 'significant_OXPHOS_pancancer_vulernabilities_across_cancers_scatter_plot.pdf'), plot = p,device = 'pdf',width =10, height = 10,dpi=300)


#plotting significant PLSR candidates across only pancancer
p <- ggplot(data=sig_candidates[sig_candidates$Cancer_Type == "Pancancer",], aes(x=Median_sub, y=Gene, fill=Cancer_Type)) + geom_point(color ="black",shape = 21, size = 4) + theme_bw() + labs(x = "Gene_dependency scores' differences (OXPHOS high cell lines - OXPHOS low cell lines)", y = "", title = "OXPHOS Elastic-Net/PLSR significant vulnerabilites in Pancancer") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="none", 
        panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=10),
        axis.line=element_line(linewidth=0.5,colour="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.5),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
        strip.text=element_text(size=8), text = element_text(size =10))  + geom_vline(xintercept = 0, colour="red", linetype = "dashed", linewidth = 0.8) + scale_fill_manual(values = cancer_pallette)

ggsave(file.path(outDir, 'significant_OXPHOS_pancancer_vulernabilities_across_pancancer_scatter_plot.pdf'), plot = p,device = 'pdf',width =8, height = 10,dpi=300)



#subplots for each cancer type/pancancer
plot_list <- list()

for (c in sort(unique(sig_candidates$Cancer_Type))) {
  
  df <- sig_candidates[sig_candidates$Cancer_Type %in% c("Pancancer", c),]
  
  #getting cell numbers for lable/title of each plot
  num_cell_lines_OXPHOS_high <- unique(df[df$Cancer_Type == c, c("num_OXPHOS_high_cell_lines")])
  num_cell_lines_OXPHOS_low <- unique(df[df$Cancer_Type == c, c("num_OXPHOS_low_cell_lines")])
  
  title1 <- paste("Significant vulnerabilites in", c)
  title2 <- paste("OXPHOS high =", num_cell_lines_OXPHOS_high)
  title3 <-  paste("OXPHOS low =", num_cell_lines_OXPHOS_low)
  title_final <- paste(title1, title2,  title3, sep = "\n")
  
  
  
  plot_list[[c]] <- ggplot(data=df, aes(x=Median_sub, y=Gene, fill=Cancer_Type)) + geom_point(color ="black",shape = 21, size = 4) + theme_bw() + labs(x = "Gene_dependency scores difference (OXPHOS high - OXPHOS low)", y = "", title = title_final) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position="none", 
          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
          axis.line=element_line(linewidth=0.5,colour="black"),
          axis.ticks = element_line(colour = "black",linewidth=0.5),
          axis.text.x=element_text(colour="black", size = 8),
          axis.text.y=element_text(colour="black", size = 8),
          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
          strip.text=element_text(size=8), text = element_text(size =8))  + geom_vline(xintercept = 0, colour="red", linetype = "dashed", linewidth = 0.8) + scale_fill_manual(values = cancer_pallette)
}

p <- grid.arrange(grobs=plot_list,ncol=4)

ggsave(file.path(outDir, 'significant_OXPHOS_pancancer_vulernabilities_scatter_plot_individual_cancers.pdf'), plot = p,device = 'pdf',width =20, height = 22,dpi=300)


#saving significant pancancer targets across pancancer/cancer types for mutation analysis
save(sig_candidates, file = "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/significant_OXPHOS_pan-cancer_targets.RData")


