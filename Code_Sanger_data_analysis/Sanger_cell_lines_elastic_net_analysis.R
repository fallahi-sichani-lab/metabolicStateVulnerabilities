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
library(fgsea)
library(umap)
library(dplyr)
library(glmnet)

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger")

#importing OXPHOS genes
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"
pathways <-  gmtPathways(pathway_file)
OXPHOS_genes <- pathways[["Oxidative phosphorylation"]]
remove(pathway_file,pathways)

#importing sanger gene expression data 
gene_exp <- read_tsv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/CellLinesProject_CompleteGeneExpression_v99_GRCh38.tsv")
OXPHOS_gene_exp <- gene_exp[gene_exp$GENE_SYMBOL %in% OXPHOS_genes,]

OXPHOS_gene_exp <- OXPHOS_gene_exp %>%
  filter(GENE_SYMBOL %in% unique(OXPHOS_gene_exp$GENE_SYMBOL)) %>%
  group_by(SAMPLE_NAME) %>%
  mutate(OXPHOS_score = mean(Z_SCORE))

remove(gene_exp)

#importing Sanger KO data
KO_effect_sanger <-as.data.frame(read_tsv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/binaryDepScores.tsv"))
rownames(KO_effect_sanger) <- KO_effect_sanger$Gene; KO_effect_sanger <- KO_effect_sanger[,-c(1)]

#importing Sanger metadata
metadata <- read_tsv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/CellLinesProject_Sample_v99_GRCh38.tsv")
metadata$StrippedCellLineName <- toupper(gsub("-","", metadata$SAMPLE_NAME))

#importing depmap metadata
DepMap_metadata <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/DepMap_23Q2_data/Model.csv")

sanger_depmap_metadata <- merge(metadata[,c(1:3,35)], DepMap_metadata[,c(4,29)], by.x = "StrippedCellLineName", by.y = "StrippedCellLineName")
sanger_depmap_metadata$OncotreePrimaryDisease <- gsub(" ", "", sanger_depmap_metadata$OncotreePrimaryDisease)

remove(DepMap_metadata)

#importing pan-cancer OXPHOS variable cancer types
OXPHOS_cancer_types <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/Cell line to cell line variability/Sanger_OXPHOS_variable_cancertypes.csv")
  
OXPHOS_cancer_types <- unique(OXPHOS_cancer_types$CELLTYPE)
OXPHOS_cancer_types <- gsub(" ", "", OXPHOS_cancer_types)

#filtering Sanger cell lines to only OXPHOS variable cancer types
sanger_depmap_metadata <- sanger_depmap_metadata[sanger_depmap_metadata$OncotreePrimaryDisease %in% OXPHOS_cancer_types,]

sanger_depmap_metadata <- unique(merge(sanger_depmap_metadata, OXPHOS_gene_exp[,c(2,8)], by.x = "SAMPLE_NAME", by.y = "SAMPLE_NAME"))

OXPHOS_low_high_KO_effect <- KO_effect_sanger[,colnames(KO_effect_sanger) %in% unique(sanger_depmap_metadata$SAMPLE_NAME)]

sanger_depmap_metadata <- sanger_depmap_metadata[sanger_depmap_metadata$SAMPLE_NAME %in% colnames(OXPHOS_low_high_KO_effect),]

#splitting data into OXPHOS high/low based on score distributions
df_quants <- as.data.frame(t(quantile(sanger_depmap_metadata$OXPHOS_score, probs = c(0.33,0.67))))

p <- ggplot(data = sanger_depmap_metadata,aes(x = OXPHOS_score))+
  geom_histogram(color="black", fill="lightpink", binwidth=0.05) + 
  xlab("Metabolic state score") +
  ylab("# cell lines") + 
  labs(title =  "Metabolic state score histogram OXPHOS (Sanger)")+
  theme_classic() + theme(legend.position="right",panel.grid.major = element_blank(), 
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=8),
                          axis.line=element_line(size=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x=element_text(colour="black", size = 12),
                          axis.text.y=element_text(colour="black", size = 12),
                          strip.background = element_rect(fill="white",linewidth = 0.2,colour = NULL),
                          strip.text=element_text(size=8), text = element_text(size = 12))  + ylim(0,20) + 
  geom_vline(xintercept = c(df_quants$`33%`,df_quants$`67%`), colour="black", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x=df_quants$`33%` - 0.5, y=40, label= paste("33% cutoff = ",round(df_quants$`33%`, digits = 2)), col="black", size=5) + 
  annotate("text", x=df_quants$`67%` + 0.5, y=40, label= paste("67% cutoff = ",round(df_quants$`67%`, digits = 2)), col="black", size=5) 

ggsave(file.path(outDir, 'OXPHOS_state_scores_histogram_Sanger.pdf'), plot = p,device = 'pdf',width =12, height = 9,dpi=300)

OXPHOS_states_filtered <- sanger_depmap_metadata[sanger_depmap_metadata$OXPHOS_score < df_quants$`33%` | sanger_depmap_metadata$OXPHOS_score > df_quants$`67%`,]
  

#filtering out genes in which average KO score is 0 or 1
OXPHOS_mean_KO <- data.frame(avg_KO = rowMeans(OXPHOS_low_high_KO_effect), gene = rownames(OXPHOS_low_high_KO_effect))
gene_KOs_0 <- OXPHOS_mean_KO[OXPHOS_mean_KO$avg_KO == 0,]$gene
gene_KOs_1 <- OXPHOS_mean_KO[OXPHOS_mean_KO$avg_KO == 1,]$gene

OXPHOS_low_high_KO_effect_filtered <- OXPHOS_low_high_KO_effect[!(rownames(OXPHOS_low_high_KO_effect) %in% c(gene_KOs_0,gene_KOs_1)),colnames(OXPHOS_low_high_KO_effect) %in% OXPHOS_states_filtered$SAMPLE_NAME]


#final dataframe for elastic net (merging KO and OXPHOS scores data)
OXPHOS_cell_lines_gene_KO.df <- merge(OXPHOS_states_filtered[,c(1,6)], 
                                      t(OXPHOS_low_high_KO_effect_filtered), by.x = "SAMPLE_NAME", by.y = 0)


#optimizing alpha for Elastic net model
ptm <- proc.time() #timing alpha optimization

list.of.fits.alpha <- list()
alpha.results <- data.frame()
  
  #filtering model to one metabolic pathway
  elastic_net.df <- OXPHOS_cell_lines_gene_KO.df
  elastic_net.df <- elastic_net.df[ , colSums(is.na(elastic_net.df))==0] #removing genes with NA values
  
  for(m in 1:150) {
    for (a in 1:9) {
      
      set.seed(m*10) #setting same seed for alpha 0.1-0.9 so datasplits are the same 
      
      #running elastic net and varying alpha, 10-fold cross validation
      fit.name <- paste0("seed", m, "alpha", a/10)
      seed.name <- paste0("seed",m)
      list.of.fits.alpha[[fit.name]]  <- cv.glmnet(as.matrix(elastic_net.df[,-c(1:2)]), elastic_net.df$OXPHOS_score, type.measure = "mse", alpha = a/10, nfolds = 10)
      
      #determining the smallest MSE model and its associated alpha and lambda.min
      mse.min <- min(list.of.fits.alpha[[fit.name]]$cvm)
      alpha.results <- rbind(alpha.results,data.frame(alpha = a/10, seed = seed.name,mse = mse.min, fit.name = fit.name, lambda.min = list.of.fits.alpha[[fit.name]]$lambda.min))
      
    }
  }
  
proc.time() - ptm  #reporting alpha optimization time

remove(fit.name,seed.name,list.of.fits.alpha,mse.min)


#boxplots of alphas (reported MSE based on lambda.min)
  p <- ggplot(alpha.results,aes(x=as.character(alpha), y = mse)) +
    geom_boxplot(fill = "purple", alpha = 0.4, size=0.4,show.legend = F, outlier.shape = NA) +
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 14,angle=45,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 14),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          text = element_text(size=14),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))+ ylab("MSE") + xlab("alpha") + labs(title = "alpha optimization (Sanger)", size = 14) + geom_jitter(color="black", size=0.5, alpha=0.7, width = 0.2) 
  

ggsave(file.path(outDir, 'boxplots_alphas_optimization_OXPHOS_elastic_net_Sanger.pdf'), plot = p,device = 'pdf',width =10, height = 6,dpi=300)



#iterating elastic net over 150 times with random 10-fold cross validation splits using optimized alpha parameter and lambda.min
set.seed(NULL)
list.of.fits <- list()
results <- data.frame()
  
  #filtering model to one metabolic pathway
  elastic_net.df <- OXPHOS_cell_lines_gene_KO.df
  elastic_net.df <- elastic_net.df[ , colSums(is.na(elastic_net.df))==0]  #removing genes with NA values
  
  #setting alpha = 0.1 (optimized alpha)
  for (a in 0.1) {
    for(m in 1:150) {
      
      #running elastic net and varying alpha
      fit.name <- paste0("model", m)
      model.name <- paste0("model",m)
      list.of.fits[[fit.name]]  <- cv.glmnet(as.matrix(elastic_net.df[,-c(1:2)]), elastic_net.df$OXPHOS_score, type.measure = "mse", alpha = a, nfolds = 10)
      
      #determining the smallest MSE model and its associated alpha.min and lambda.min
      mse.min <- min(list.of.fits[[fit.name]]$cvm)
      results <- rbind(results,data.frame(alpha = a, model = model.name,mse = mse.min, fit.name = fit.name, lambda.min = list.of.fits[[fit.name]]$lambda.min))
      

  }
}
  
remove(fit.name,model.name,mse.min)


#extracting coefficients for each model
elastic_net_model_coefficients <- data.frame()
model.num <- 0
for (r in results$fit.name) {
  
  alpha.min_lambda.min_coeffs <- coef(list.of.fits[[r]], s = list.of.fits[[r]]$lambda.min)
  alpha.min_lambda.min_coeffs.df <- data.frame(name = alpha.min_lambda.min_coeffs@Dimnames[[1]][alpha.min_lambda.min_coeffs@i + 1], coefficient =alpha.min_lambda.min_coeffs@x)
  alpha.min_lambda.min_coeffs.df <-  alpha.min_lambda.min_coeffs.df[-c(1),] #removing intercept
  
  model.num <- model.num + 1
  alpha.min_lambda.min_coeffs.df$model <- model.num
  
  elastic_net_model_coefficients <- rbind(elastic_net_model_coefficients, alpha.min_lambda.min_coeffs.df)
  
}


#features selected across 150 random elastic models with alpha = 0.1 and lamba = lambda.min
filtered_coefficeints_ENR_metabolic_pathways <- data.frame()
  for (m in c(150)) {
    
    #determining number of times a variable appears in the 150 models
    counts <- elastic_net_model_coefficients[elastic_net_model_coefficients$model %in% 1:m, ]
    counts <- as.data.frame(table(counts$name))
    counts$fraction <- counts$Freq/m
    
    # lock in factor level order
    counts <- counts[order(counts$fraction, decreasing = TRUE),] #ordering coefficient values
    counts$Var1 <- factor(counts$Var1 , levels = counts$Var1 )
    
    list_label = paste(m)
    p <-ggplot(counts, aes(x=1:nrow(counts), y=fraction, group = 1)) + labs(title =  paste(m, "iterations"))+ 
      geom_point(size = 1) + geom_line()+theme_classic()  + theme(legend.position="none",
                                                                  axis.text.x=element_text(colour="black", size = 10),
                                                                  #axis.ticks.x=element_blank(),
                                                                  axis.text.y=element_text(colour="black", size = 10),
                                                                  axis.line=element_line(size=0.2,color="black"),
                                                                  axis.ticks = element_line(colour = "black",size=0.2),
                                                                  panel.border = element_blank(), panel.background = element_blank(),
                                                                  axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + 
                                                                  ylab("Fraction of occurence") + xlab("Genes") + 
                                                                  geom_hline(yintercept = c(0.5), colour="red", linetype = "dashed", linewidth = 0.6) + scale_x_continuous(breaks = seq(0,nrow(counts), by = 20))
    
    #selecting variables that appear in at least 50% of the models
    if (m == 150) {
      filtered_coefficients <- counts[counts$fraction >= 0.5, ]
      filtered_coefficeints_ENR_metabolic_pathways <- rbind(filtered_coefficeints_ENR_metabolic_pathways, filtered_coefficients)
    }
    
  }


ggsave(file.path(outDir, 'elastic_net_150_model_iterations_fraction_0.5_Sanger.pdf'), plot = p,device = 'pdf',width =12, height = 8,dpi=300)


  

#calculating mean and std. deviation for each gene's coefficient 
  coefficient_mean_std.df <- data.frame()
  pathway_ENR_genes <- filtered_coefficeints_ENR_metabolic_pathways$Var1
  df <- elastic_net_model_coefficients
  
  for (g in pathway_ENR_genes) {
    
    df2 <- df[df$name == g, ]
    df.mean <- mean(df2$coefficient)
    df.sd <- sd(df2$coefficient)
    
    coefficient_mean_std.df <- rbind(coefficient_mean_std.df, data.frame(gene = g, mean_coef = df.mean, std_coef = df.sd))
  }
  
  remove(pathway_ENR_genes,df,df2,df.mean,df.sd,g)

# lock in factor level order
coefficient_mean_std.df <- coefficient_mean_std.df[order(coefficient_mean_std.df$mean_coef),] #ordering coefficient values
coefficient_mean_std.df$gene <- factor(coefficient_mean_std.df$gene , levels = coefficient_mean_std.df$gene)


#exporting coefficient results
write.csv(coefficient_mean_std.df, file.path(outDir, "Sanger_elastic_net_coefficients.csv"))

#plotting coefficients/variables identified from elastic net regression 
  p <-ggplot(coefficient_mean_std.df, aes(x=mean_coef, y=gene)) +
    geom_bar(stat="identity", fill="forestgreen", alpha = 0.7, width = 0.8)+theme_classic()  + theme(legend.position="none",
                                                                                                     axis.text.x=element_text(colour="black", size = 10),
                                                                                                     axis.text.y=element_text(colour="black", size = 4, angle = 30),
                                                                                                     axis.line=element_line(size=0.2,color="black"),
                                                                                                     axis.ticks = element_line(colour = "black",size=0.2),
                                                                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                                                                     axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + 
                                                                                                      ylab("Genes") + xlab("Coefficient") +
                                                                                                      geom_errorbar( aes(xmin= mean_coef-std_coef, xmax = mean_coef+std_coef,y=gene), width=0.4, colour="black", size=0.3) + 
                                                                                                labs(title = "Selected elastic net variables (Sanger)")
  

ggsave(file.path(outDir, 'elastic_net_coefficients_OXPHOS_model_Sanger.pdf'), plot = p,device = 'pdf',width = 7, height = 10,dpi=300)


  