rm(list=ls()) #clear all
cat("\014") #clc


## Libraries Required
library(ggplot2)
library(ggpubr)
library(purrr)
library(factoextra)
library(dplyr)
library(jcolors)
library(umap)
library(randomcoloR)
library(NbClust)
library(gridExtra)
library(fgsea)
library(tidyverse)

#load in KEGG metabolism data
KEGG_metabolic_pathways <- gmtPathways("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt")

KEGG_metabolic_genes = as.data.frame(unique(unlist(KEGG_metabolic_pathways)))

colnames(KEGG_metabolic_genes) = c("Metabolic_gene")

#loading z-scored TPM data,media, and mutation data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")



#filtering genes for only metabolic KEGG genes
CCLE_TPM_metabolic_genes = as.data.frame(CCLE_TPM_zscore[,colnames(CCLE_TPM_zscore) %in% KEGG_metabolic_genes$Metabolic_gene])

#selecting blood cancers 
blood_cancers <- CCLE_annotation[CCLE_annotation$OncotreeLineage %in% c("Myeloid", "Lymphoid"),]
blood_cells_metabolic_genes <- CCLE_TPM_metabolic_genes[rownames(CCLE_TPM_metabolic_genes) %in% unique(blood_cancers$StrippedCellLineName),]

#run PCA 
pca <- prcomp(blood_cells_metabolic_genes, scale = FALSE)
loadings <- pca$rotation
scores <- pca$x

#calculating % variance
eigs <- pca$sdev^2
data_variance_explained <- 100*(eigs / sum(eigs))
cum_var <- cumsum(data_variance_explained)
select_pcs <- which(cum_var>=80)[1] #selecting scores with 80% variance explained

PCs<- merge(scores[,c(1:select_pcs)], blood_cancers[,c(4,29,30)], by.x = 0, by.y = "StrippedCellLineName")


#plotting PCA across lineage types
plot_embedding_categorical <- function(category,.df){
  dpal_func = jcolors_contin("pal4",reverse = TRUE,bias = 0.5)
  dpal = rev(dpal_func(length(unique(.df[[category]]))+1))
  ggplot(data = .df,aes(x = `PC1`,y = `PC2`,fill = .data[[category]]))+
    geom_point(stroke = 0.2 ,shape = 21, size = 6)+
    guides(fill= guide_legend(override.aes = list(size=3)))+
    xlab(paste("PC1", "(", round(data_variance_explained[1], digits = 2), "%)")) +
    ylab(paste("PC2", "(", round(data_variance_explained[2], digits = 2), "%)")) +
    theme_classic() + scale_fill_manual(values = c("#D95F02","#7570B3"))
}


g = plot_embedding_categorical('OncotreeLineage',PCs)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/PCA_blood_lineage_type.pdf',plot = g,device = 'pdf',width = 9, height = 7,dpi=300)


#plotting PCA across cancer types
plot_embedding_categorical <- function(category,.df){
  dpal_func = jcolors_contin("pal4",reverse = TRUE,bias = 0.5)
  dpal = rev(dpal_func(length(unique(.df[[category]]))+1))
  ggplot(data = .df,aes(x = `PC1`,y = `PC2`,fill = .data[[category]]))+
    geom_point(stroke = 0.2 ,shape = 21, size = 6)+
    guides(fill= guide_legend(override.aes = list(size=3)))+
    xlab(paste("PC1", "(", round(data_variance_explained[1], digits = 2), "%)")) +
    ylab(paste("PC2", "(", round(data_variance_explained[2], digits = 2), "%)")) +
    theme_classic() + scale_fill_manual(values = palette_final[c(2,4,8,11,20,22)])
}


g = plot_embedding_categorical('OncotreePrimaryDisease',PCs)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/PCA_blood_cancer_type.pdf',plot = g,device = 'pdf',width = 9, height = 7,dpi=300)
