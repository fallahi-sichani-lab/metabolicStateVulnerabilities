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

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/Elastic_net_OXPHOS")

#output directory for elastic net data to be used for PLSR modeling
outDirPLSR <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS")


#loading data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/Elastic_net_OXPHOS/cell_line_to_cell_line_metabolic_variability_scores.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/organized_CRISPR_gene_effect.RData")

#cleaning up CRISPR dataset
cell_lines_CRISPR <- data.frame(rownames(CRISPR_CCLE_data),CRISPR_CCLE_data)
rownames(cell_lines_CRISPR) <- NULL
colnames(cell_lines_CRISPR)[1] <- "cell_line"

#exploring missing data in CRISPR dataset, removing cell lines that do not have KO data
genes_NA <- as.data.frame(colSums(is.na(cell_lines_CRISPR[,c(2:17932)]))) 
genes_NA$empty <- 0
genes_NA <-genes_NA[genes_NA$`colSums(is.na(cell_lines_CRISPR[, c(2:17932)]))` == 0,]
genes_CRISPR_removed_NA<- cell_lines_CRISPR[,c("cell_line", rownames(genes_NA))]


#filtering CRISPR data to only selected pathways' cell lines (i.e. oxidative phosphorylation/OXPHOS)
metabolic_pathways <- "Oxidative phosphorylation"

Gene_dependency_scores_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  metabolic_pathway_cell_lines <- cell_line_metabolic_variability.df[cell_line_metabolic_variability.df$Pathway == p,]
  metabolic_score_Gene_dependency <- merge(metabolic_pathway_cell_lines,genes_CRISPR_removed_NA, by = "cell_line") 
  
  Gene_dependency_scores_metabolic_pathways <- rbind(Gene_dependency_scores_metabolic_pathways, metabolic_score_Gene_dependency)
}


#plotting histogram of OXPHOS state scores with 33/67 percentiles 
plot_list <- list()
for (p in metabolic_pathways) {
  df <- Gene_dependency_scores_metabolic_pathways[Gene_dependency_scores_metabolic_pathways$Pathway == p, ]
  df_quants <- as.data.frame(t(quantile(df$mean_zscore, probs = c(0.33,0.67))))
  
  plot_list[[p]]  <- ggplot(data = df,aes(x = mean_zscore))+
    geom_histogram(color="black", fill="lightpink", binwidth=0.05) + 
    xlab("Metabolic state score") +
    ylab("# cell lines") + 
    labs(title =  paste("Metabolic state score histogram ",p, sep = "\n"))+
    theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(), 
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
                            axis.line=element_line(size=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",linewidth=0.5),
                            axis.text.x=element_text(colour="black", size = 12),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",linewidth = 0.2,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12))  + ylim(0,50) + geom_vline(xintercept = c(df_quants$`33%`,df_quants$`67%`), colour="black", linetype = "dashed", linewidth = 0.9) +
                      annotate("text", x=df_quants$`33%` - 0.5, y=40, label= paste("33% cutoff = ",round(df_quants$`33%`, digits = 2)), col="black", size=5) + annotate("text", x=df_quants$`67%` + 0.5, y=40, label= paste("67% cutoff = ",round(df_quants$`67%`, digits = 2)), col="black", size=5) 
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'OXPHOS_state_scores_histogram.pdf'), plot = g,device = 'pdf',width =12, height = 9,dpi=300)



#filtering datasets to remove samples with metabolic state scores close to 0 based on 33% and 67% distribution thresholds
Gene_dependency_scores_metabolic_pathways_filtered <- data.frame()
for (p in metabolic_pathways) {
  df <- Gene_dependency_scores_metabolic_pathways[Gene_dependency_scores_metabolic_pathways$Pathway == p, ]
  df_quants <- as.data.frame(t(quantile(df$mean_zscore, probs = c(0.33,0.67))))
  
  df_filtered <- df[df$mean_zscore < df_quants$`33%` | df$mean_zscore > df_quants$`67%`,]
  
  
  Gene_dependency_scores_metabolic_pathways_filtered <- rbind(Gene_dependency_scores_metabolic_pathways_filtered, df_filtered)
}



#calculating |IQR| for each gene across OXPHOS high and low cell lines, want broad distribution of gene dependency scores to be broad
IQR_all_genes_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  
  df <- Gene_dependency_scores_metabolic_pathways_filtered[Gene_dependency_scores_metabolic_pathways_filtered$Pathway == p, ]
  
  IQR_each_gene <- t(as.data.frame(apply(df[,-c(1:4)], 2, quantile)))
  
  
  IQR_each_gene  <- as.data.frame(IQR_each_gene[,c(2:4)])
  colnames(IQR_each_gene) <- c("25_percentile", "50_percentile", "75_percentile")
  IQR_each_gene$IQR <- abs(IQR_each_gene$`75_percentile` - IQR_each_gene$`25_percentile`)
  IQR_each_gene$Pathway <- p
  IQR_each_gene$gene <- rownames(IQR_each_gene)
  rownames(IQR_each_gene) <- NULL
  
  IQR_all_genes_metabolic_pathways <- rbind(IQR_all_genes_metabolic_pathways, IQR_each_gene)
}


#establishing filters for IQR 
OXPHOS_IQR = 0.09


#plotting IQRs of all genes across metabolic pathways
plot_list <- list()
for (p in metabolic_pathways) {
  df <- IQR_all_genes_metabolic_pathways[IQR_all_genes_metabolic_pathways$Pathway == p, ]
  
  number_genes_remain_IQR <- round((nrow(df[df$IQR > OXPHOS_IQR, ])/nrow(df))*100, digits = 1)
  
  plot_list[[p]]  <- ggplot(data = df,aes(x = IQR))+
    geom_histogram(color="black", fill="lightblue", binwidth=0.005) + 
    xlab("IQR") +
    ylab("# genes") + 
    labs(title =  paste("Gene IQR histogram ",p, sep = "\n"))+
    theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(), 
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
                            axis.line=element_line(size=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",size=0.5),
                            axis.text.x=element_text(colour="black", size = 12),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",linewidth=0.2,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12)) + scale_x_continuous(breaks = seq(0, max(df$IQR), 0.1)) + scale_y_continuous(breaks = seq(0, 1800, 200)) + geom_vline(xintercept = OXPHOS_IQR, colour="red", linetype = "dashed", linewidth = 0.8) +
                      annotate("text", x=OXPHOS_IQR+0.1, y=750, label= paste(number_genes_remain_IQR, "% genes remain"), col="red", size=5) + annotate("text", x=OXPHOS_IQR+0.1, y=650, label= paste("IQR threshold =", OXPHOS_IQR), col="red", size=5)
    
    
    
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'IQR_histogram.pdf'), plot = g,device = 'pdf',width =12, height = 6,dpi=300)


#filtering data based on IQR threshold
filtered_genes <- IQR_all_genes_metabolic_pathways[IQR_all_genes_metabolic_pathways$IQR >  OXPHOS_IQR, ] #filtering based on IQR 
filtered_genes_across_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  genes  <- filtered_genes[filtered_genes$Pathway== p,]$gene
  print(length(genes))
  
  df2 <- Gene_dependency_scores_metabolic_pathways_filtered[Gene_dependency_scores_metabolic_pathways_filtered$Pathway == p, ]
  df2 <- df2[,colnames(df2) %in% c("mean_zscore", "Pathway", "Cancer_Type",genes)]
  
  filtered_genes_across_metabolic_pathways <- rbind.fill(filtered_genes_across_metabolic_pathways, df2) #rbind dataframes and input of NA for any missing data
}



#separating 90% of the data for training the model for elastic net/PLSR, 10% for testing for independent validation post elastic net/PLSR
train_datasets_metabolic_pathways <- data.frame()
test_datasets_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  
  metabolic_pathway.df <- filtered_genes_across_metabolic_pathways[filtered_genes_across_metabolic_pathways$Pathway == p, ]
  positive_Gene_dependency_score <- metabolic_pathway.df[metabolic_pathway.df$mean_zscore > 0,] #splitting samples into high/low OXPHOS scores to select samples for test dataset
  negative_Gene_dependency_score <- metabolic_pathway.df[metabolic_pathway.df$mean_zscore < 0,]
  
  number_observations_for_train <- 2 * round((nrow(metabolic_pathway.df)*0.1)/2) #10% for test dataset
  
  #random selection of cell lines for train/test datasets that have positive OXPHOS state scores
  set.seed(200)
  Gene_dependency_positive_samples <- sample(1:nrow(positive_Gene_dependency_score), 0.1*nrow(positive_Gene_dependency_score))
  test_positive_samples <- positive_Gene_dependency_score[Gene_dependency_positive_samples,]
  train_positive_samples <- positive_Gene_dependency_score[-Gene_dependency_positive_samples,]
  
  
  #random selection of cell lines for train/test datasets that have negative OXPHOS state scores
  set.seed(100)
  Gene_dependency_negative_samples <- sample(1:nrow(negative_Gene_dependency_score), 0.1*nrow(negative_Gene_dependency_score))
  test_negative_samples <- negative_Gene_dependency_score[Gene_dependency_negative_samples,]
  train_negative_samples <- negative_Gene_dependency_score[-Gene_dependency_negative_samples,]
  
  
  #combining positive and negative cell lines for train and test datasets
  train_metabolic_pathway.df <- rbind(train_positive_samples, train_negative_samples)
  test_metabolic_pathway.df <- rbind(test_positive_samples,test_negative_samples)
  
  train_datasets_metabolic_pathways <- rbind(train_datasets_metabolic_pathways, train_metabolic_pathway.df)
  test_datasets_metabolic_pathways <- rbind(test_datasets_metabolic_pathways, test_metabolic_pathway.df)
  
}


#histograms of test and train datasets
plot_list <- list()
for (p in metabolic_pathways) {
  df1 <- train_datasets_metabolic_pathways[train_datasets_metabolic_pathways$Pathway == p, ]
  df2 <- test_datasets_metabolic_pathways[test_datasets_metabolic_pathways$Pathway == p, ]
  
  plot_list[[p]]  <- ggplot() + geom_histogram(data = df1,aes(x = mean_zscore), color="black", fill="black", binwidth=0.05, alpha = 0.7) +
    geom_histogram(data = df2,aes(x = mean_zscore), color="black", fill="lightpink", binwidth=0.05) +
    xlab("Metabolic state score") +
    ylab("# cell lines") + 
    labs(title =  paste("Metabolic state score histogram ",p, sep = "\n"))+
    theme_classic() + theme(legend.position="right", panel.grid.major = element_blank(), 
                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
                            axis.line=element_line(linewidth=0.5,colour="black"),
                            axis.ticks = element_line(colour = "black",linewidth =0.5),
                            axis.text.x=element_text(colour="black", size = 12),
                            axis.text.y=element_text(colour="black", size = 12),
                            strip.background = element_rect(fill="white",size=0.2,colour = NULL),
                            strip.text=element_text(size=8), text = element_text(size = 12)) 
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'train_test_datasets_histograms.pdf'), plot = g,device = 'pdf',width =12, height = 6,dpi=300)


#optimizing alpha for Elastic net model
ptm <- proc.time() #timing alpha optimization

list.of.fits.alpha <- list()
alpha.results <- data.frame()
for (p in metabolic_pathways) {
  
  #filtering model to one metabolic pathway
  elastic_net.df <- train_datasets_metabolic_pathways[train_datasets_metabolic_pathways$Pathway == p, ]
  elastic_net.df <- elastic_net.df[ , colSums(is.na(elastic_net.df))==0] #removing genes with NA values
  
 for(m in 1:150) {
  for (a in 1:9) {
      
      set.seed(m*10) #setting same seed for alpha 0.1-0.9 so datasplits are the same 
    
      #running elastic net and varying alpha, 10-fold cross validation
      fit.name <- paste0(p, "seed", m, "alpha", a/10)
      seed.name <- paste0("seed",m)
      list.of.fits.alpha[[fit.name]]  <- cv.glmnet(as.matrix(elastic_net.df[,-c(1:3)]), elastic_net.df$mean_zscore, type.measure = "mse", alpha = a/10, nfolds = 10)
      
      #determining the smallest MSE model and its associated alpha and lambda.min
      mse.min <- min(list.of.fits.alpha[[fit.name]]$cvm)
      alpha.results <- rbind(alpha.results,data.frame(alpha = a/10, seed = seed.name,mse = mse.min, fit.name = fit.name, lambda.min = list.of.fits.alpha[[fit.name]]$lambda.min, pathway = p))
      
    }
  }
}

proc.time() - ptm  #reporting alpha optimization time



#boxplots of alphas for each metabolic pathway model (reported MSE based on lambda.min)
plot_list <- list()
for (p in metabolic_pathways) {
  df <- alpha.results[alpha.results$pathway == p, ]
  
  
  plot_list[[p]] <- ggplot(df,aes(x=as.character(alpha), y = mse)) +
    geom_boxplot(fill = "purple", alpha = 0.4, size=0.4,show.legend = F, outlier.shape = NA) +
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 14,angle=45,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 14),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          text = element_text(size=14),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))+ ylab("MSE") + xlab("alpha") + labs(title = p, size = 14) + geom_jitter(color="black", size=0.5, alpha=0.7, width = 0.2) 
  
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'boxplots_alphas_optimization_OXPHOS_elastic_net.pdf'), plot = g,device = 'pdf',width =10, height = 6,dpi=300)


#iterating elastic net over 150 times with random 10-fold cross validation splits using optimized alpha parameter and lambda.min
set.seed(NULL)
list.of.fits <- list()
results <- data.frame()
for (p in metabolic_pathways) {
  
  #filtering model to one metabolic pathway
  elastic_net.df <- train_datasets_metabolic_pathways[train_datasets_metabolic_pathways$Pathway == p, ]
  elastic_net.df <- elastic_net.df[ , colSums(is.na(elastic_net.df))==0] #removing genes with NA values
  
  #setting alpha = 0.1 (optimized alpha)
  for (a in 0.1) {
    for(m in 1:150) {
    
    #running elastic net and varying alpha
    fit.name <- paste0(p, "model", m)
    model.name <- paste0("model",m)
    list.of.fits[[fit.name]]  <- cv.glmnet(as.matrix(elastic_net.df[,-c(1:3)]), elastic_net.df$mean_zscore, type.measure = "mse", alpha = a, nfolds = 10)
    
    #determining the smallest MSE model and its associated alpha.min and lambda.min
    mse.min <- min(list.of.fits[[fit.name]]$cvm)
    results <- rbind(results,data.frame(alpha = a, model = model.name,mse = mse.min, fit.name = fit.name, lambda.min = list.of.fits[[fit.name]]$lambda.min, pathway = p))
    
    }
  }
}

#extracting coefficients for each model
elastic_net_model_coefficients <- data.frame()
model.num <- 0
for (r in results$fit.name) {
  
  alpha.min_lambda.min_coeffs <- coef(list.of.fits[[r]], s = list.of.fits[[r]]$lambda.min)
  alpha.min_lambda.min_coeffs.df <- data.frame(name = alpha.min_lambda.min_coeffs@Dimnames[[1]][alpha.min_lambda.min_coeffs@i + 1], coefficient =alpha.min_lambda.min_coeffs@x)
  alpha.min_lambda.min_coeffs.df <-  alpha.min_lambda.min_coeffs.df[-c(1),] #removing intercept
  
  model.num <- model.num + 1
  alpha.min_lambda.min_coeffs.df$model <- model.num
  alpha.min_lambda.min_coeffs.df$pathway <- results[results$fit.name == r, ]$pathway
  
 elastic_net_model_coefficients <- rbind(elastic_net_model_coefficients, alpha.min_lambda.min_coeffs.df)
  
}


#features selected across 150 random elastic models with alpha = 0.1 and lamba = lambda.min
plot_list <- list()
filtered_coefficeints_ENR_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  for (m in c(150)) {
  
  #determining number of times a variable appears in the 150 models
  counts <- elastic_net_model_coefficients[elastic_net_model_coefficients$pathway == p & elastic_net_model_coefficients$model %in% 1:m, ]
  counts <- as.data.frame(table(counts$name))
  counts$fraction <- counts$Freq/m
  
  # lock in factor level order
  counts <- counts[order(counts$fraction, decreasing = TRUE),] #ordering coefficient values
  counts$Var1 <- factor(counts$Var1 , levels = counts$Var1 )
  
  list_label = paste(p,m)
  plot_list[[list_label]] <-ggplot(counts, aes(x=1:nrow(counts), y=fraction, group = 1)) + labs(title =  paste(p,paste(m, "iterations"), sep = "\n"))+ 
    geom_point(size = 1) + geom_line()+theme_classic()  + theme(legend.position="none",
                                                                                        axis.text.x=element_text(colour="black", size = 10),
                                                                                        #axis.ticks.x=element_blank(),
                                                                                        axis.text.y=element_text(colour="black", size = 10),
                                                                                        axis.line=element_line(size=0.2,color="black"),
                                                                                        axis.ticks = element_line(colour = "black",size=0.2),
                                                                                        panel.border = element_blank(), panel.background = element_blank(),
                                                                                        axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("Fraction of occurence") + xlab("Genes") + 
                                                                geom_hline(yintercept = c(0.5), colour="red", linetype = "dashed", size = 0.6) + scale_x_continuous(breaks = seq(0,nrow(counts), by = 20))
  
  #selecting variables that appear in at least 50% of the models
  if (m == 150) {
    filtered_coefficients <- counts[counts$fraction >= 0.5, ]
    filtered_coefficients$pathway = p
    filtered_coefficeints_ENR_metabolic_pathways <- rbind(filtered_coefficeints_ENR_metabolic_pathways, filtered_coefficients)
  }
  
  }
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'elastic_net_150_model_iterations_fraction_0.5.pdf'), plot = g,device = 'pdf',width =12, height = 8,dpi=300)


#calculating mean and std. deviation for each gene's coefficient 
coefficient_mean_std.df <- data.frame()
for (p in metabolic_pathways) {
  pathway_ENR_genes <- filtered_coefficeints_ENR_metabolic_pathways[filtered_coefficeints_ENR_metabolic_pathways$pathway == p, ]$Var1
  df <- elastic_net_model_coefficients[elastic_net_model_coefficients$pathway == p, ]
  
  for (g in pathway_ENR_genes) {
  
  df2 <- df[df$name == g, ]
  df.mean <- mean(df2$coefficient)
  df.sd <- sd(df2$coefficient)
  
  coefficient_mean_std.df <- rbind(coefficient_mean_std.df, data.frame(gene = g, mean_coef = df.mean, std_coef = df.sd, pathway = p))
  }
}

# lock in factor level order
coefficient_mean_std.df <- coefficient_mean_std.df[order(coefficient_mean_std.df$mean_coef),] #ordering coefficient values
coefficient_mean_std.df$gene <- factor(coefficient_mean_std.df$gene , levels = coefficient_mean_std.df$gene)

#export elastic net mean coefficient and std
fwrite(coefficient_mean_std.df,file.path(outDir,"elastic_net_selected_variables_OXPHOS_model.csv"))


#plotting coefficients/variables identified from elastic net regression 
plot_list <- list()
for (p in metabolic_pathways) {
  
  df <- coefficient_mean_std.df
  
  plot_list[[p]] <-ggplot(df, aes(x=mean_coef, y=gene)) +
    geom_bar(stat="identity", fill="forestgreen", alpha = 0.7, width = 0.8)+theme_classic()  + theme(legend.position="none",
                                                       axis.text.x=element_text(colour="black", size = 10),
                                                       axis.text.y=element_text(colour="black", size = 4, angle = 30),
                                                       axis.line=element_line(size=0.2,color="black"),
                                                       axis.ticks = element_line(colour = "black",size=0.2),
                                                       panel.border = element_blank(), panel.background = element_blank(),
                                                       axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + ylab("Genes") + xlab("Coefficient") +
                              geom_errorbar( aes(xmin= mean_coef-std_coef, xmax = mean_coef+std_coef,y=gene), width=0.4, colour="black", size=0.3) + labs(title = p)
  
}
g <- grid.arrange(grobs=plot_list,ncol=1)
ggsave(file.path(outDir, 'elastic_net_coefficients_OXPHOS_model.pdf'), plot = g,device = 'pdf',width = 5, height = 8,dpi=300)


#exporting train and testing datasets with elastic net selected variables
training_datasets_final <- data.frame()
test_datasets_final<- data.frame()
for (p in metabolic_pathways) {
  variables <- filtered_coefficeints_ENR_metabolic_pathways[filtered_coefficeints_ENR_metabolic_pathways$pathway == p, ]$Var1
  train.df <- train_datasets_metabolic_pathways[train_datasets_metabolic_pathways$Pathway == p, c("mean_zscore", "Pathway", "Cancer_Type", as.character(variables))]
  test.df <- test_datasets_metabolic_pathways[test_datasets_metabolic_pathways$Pathway == p, c("mean_zscore", "Pathway", "Cancer_Type", as.character(variables))]
  
  training_datasets_final <- rbind(training_datasets_final, train.df)
  test_datasets_final <- rbind(test_datasets_final, test.df)
}


#export training and test datasets for PLSR model
fwrite(training_datasets_final,file.path(outDirPLSR,"metabolic_pathways_training_dataset.csv"))
fwrite(test_datasets_final,file.path(outDirPLSR,"metabolic_pathways_test_dataset.csv"))


