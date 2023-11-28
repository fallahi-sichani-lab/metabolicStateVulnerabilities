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
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE_mutation_data.RData")



#filtering genes for only metabolic KEGG genes
CCLE_TPM_metabolic_genes = as.data.frame(CCLE_TPM_zscore[,colnames(CCLE_TPM_zscore) %in% KEGG_metabolic_genes$Metabolic_gene])


#Set up R UMAP, optimizing parameters
for (nn in c(50)){
  for (md in c(0.5)){
      for (metric in c("pearson")){

seed = 100
ncomp = 2

custom.config = umap.defaults
custom.config$random_state = seed
custom.config$n_neighbors = as.integer(nn)
custom.config$n_components = as.integer(ncomp)
custom.config$metric = metric
custom.config$min_dist = md
print(custom.config)

# run umap
set.seed(seed)

umap_output <- umap(as.matrix(CCLE_TPM_metabolic_genes),config = custom.config)


#Collect Output
output = list(umap = umap_output$layout,
              data.zscore = CCLE_TPM_metabolic_genes,
              parameters = data.frame(nn=nn,md=md,ncomp=ncomp,metric=metric))


df = output[["data.zscore"]]
df$X1=output[["umap"]][,1]
df$X2=output[["umap"]][,2]
df$cancer_type = CCLE_cell_line_tumor_type$OncotreePrimaryDisease



#plotting UMAP across all cancer subtypes
plot_embedding_categorical <- function(category,.df){
  dpal_func = jcolors_contin("pal4",reverse = TRUE,bias = 0.5)
  dpal = rev(dpal_func(length(unique(.df[[category]]))+1))
  ggplot(data = .df,aes(x = `X1`,y = `X2`,fill = .data[[category]]))+
    geom_point(stroke = 0.2 ,shape = 21, size = 4)+
    guides(fill= guide_legend(override.aes = list(size=3)))+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("Cancer Type")+
    theme_classic() +scale_fill_manual(values = palette_final)
}


g = plot_embedding_categorical('cancer_type',df)
annotate_figure(g, fig.lab = "Metabolic clustering CCLE cell lines", fig.lab.face = "bold",fig.lab.size = 14,fig.lab.pos = "top.left")
ggsave(sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/test_UMAP-nn=%0.0f-md=%0.0e-%s.pdf',nn,md,metric),plot = g,device = 'pdf',width = 15, height = 7,dpi=300)

      }
    }
}


#UMAP subplots for each cancer type
cancer_list = unique(df$cancer_type)

plot_list <- list() 
for (i in 1:length(cancer_list)){
  kk = cancer_list[i]
  color = palette_final[which(cancer_list == kk)]
  
  
  df2 = filter(df,cancer_type%in%kk)
  
  plot_list[[i]] <-ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 0.1,color = "#c6c6c6")  +
    geom_point(data = df2,aes(x = X1, y = X2,fill = color),color ="black",shape = 21,stroke = 0.2, size = 4)+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(title = kk)+
    scale_fill_manual(values = color)+
    theme_classic() + theme(legend.position="none", plot.title = element_text(size=12))

}
g <- grid.arrange(grobs=plot_list,ncol=9)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/UMAP_cancer_type.pdf',plot = g,device = 'pdf',width = 30, height = 18,dpi=300)


#plotting example cancer types UMAPs
for (kk in c("Non-Small Cell Lung Cancer", "Melanoma", "Acute Myeloid Leukemia", "Ovarian Epithelial Tumor", "Lung Neuroendocrine Tumor", "Neuroblastoma", "Ocular Melanoma")) {
  
  color = palette_final[which(cancer_list == kk)]
  df2 = filter(df,cancer_type%in%kk)
  
  p <-ggplot(df, aes(x = X1, y = X2)) +
  geom_point(size = 2,color = "#c6c6c6")  +
  geom_point(data = df2,aes(x = X1, y = X2,fill = color),color ="black",shape = 21,stroke = 0.2, size = 6)+
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  labs(title = kk)+
  scale_fill_manual(values = color)+
  theme_classic() + theme(legend.position="none")
  
  ggsave(sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/UMAP_%s.pdf',kk),plot = p,device = 'pdf',width = 9, height = 9,dpi=300)
}


#plotting different cancer types under the umbrella of NSCLC cancer type
for (kk in c("Non-Small Cell Lung Cancer")) {

  df2 = filter(df,cancer_type %in% kk)
  df3 <- CCLE_annotation[CCLE_annotation$OncotreePrimaryDisease == kk, ]
  rownames(df3) <- df3$StrippedCellLineName
  merged_df <- merge(df3, df2, by = 0)

  p <-ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 2,color = "#c6c6c6")  +
    geom_point(data = merged_df,aes(x = X1, y = X2,fill = OncotreeSubtype),color ="black",shape = 21,stroke = 0.2, size = 6)+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(title = kk)+
    theme_classic() + theme(legend.position="right")

  ggsave(sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/UMAP_%s_cancer_classifcations.pdf',kk),plot = p,device = 'pdf',width = 11, height = 9, dpi=300)
}

#UMAP of NSCLC mutations
plot_list <- list()
NSCLC_mutations <- c("EGFR", "KRAS", "TP53")
mutation_palette <- c("blue", "red", "darkgreen")
for (i in 1:length(NSCLC_mutations)){
  
  m <- NSCLC_mutations[i]
  kk <- c("Non-Small Cell Lung Cancer")
  
  df2 = filter(df,cancer_type %in% kk)
  df3 <-CCLE_mutation_cell_line[CCLE_mutation_cell_line$OncotreePrimaryDisease ==kk & CCLE_mutation_cell_line$HugoSymbol == m, c(4,8,9,11)]
  
  #keeping mutations that are likely gain of function or loss of function
  df3 <- unique(df3[!(df3$LikelyGoF == "False" & df3$LikelyLoF == "False"),])
  
  rownames(df3) <- df3$StrippedCellLineName
  merged_df <- merge(df3, df2, by = 0)
  
  plot_list[[i]] <-ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 2,color = "#c6c6c6")  +
    geom_point(data = merged_df,aes(x = X1, y = X2,fill = HugoSymbol),color ="black",shape = 21,stroke = 0.2, size = 4)+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(title = paste(kk, m, "n =", nrow(merged_df), "/", nrow(df2), "cell lines"))+
    scale_fill_manual(values = mutation_palette[i])+
    theme_classic() + theme(legend.position="none", plot.title = element_text(size=14))
  
}
g <- grid.arrange(grobs=plot_list,ncol=3)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/NSCLC_mutations_UMAP.pdf',plot = g,device = 'pdf',width = 15, height = 6,dpi=300)

#####
#####

#UMAP analysis with different growth medias

#UMAP with media compositions overlaid on top
count_medium = as.data.frame(table(medium_metadata_unfiltered$FormulationID))
count_medium = count_medium[count_medium$Freq >= 10,] #filering medium with at least 10 cell lines cultured in 

#keeping only culture media in which at least 10 cell lines are cultured in 
df2 = medium_metadata_unfiltered[medium_metadata_unfiltered$FormulationID %in% count_medium$Var1,]
df2 <- unique(df2)
df2 <- df2[!df2$StrippedCellLineName == "JIMT1",] #removing redudant cell lines from data
rownames(df2) = df2$StrippedCellLineName
merged_df <- merge(df2, df, by = 0)

#color palette
new_palette <- c("darkorange", "red", "blue", "hotpink", "darkorchid4", "green", "goldenrod1","darkgray", "brown4", "lightblue","lightpink","cyan", "yellow", "mediumpurple1", 
                 "darkgreen", "lightsalmon3", "lightgreen", "coral2", "lightcyan3", "firebrick2", "blueviolet")

#subplots of UMAPs with different growth media conditions
plot_list <- list() 
for (i in 1:length(count_medium$Var1)) {
  
  color <- new_palette[i]
  m <- count_medium$Var1[i]
  df3 <- merged_df[merged_df$FormulationID == m,]
  
  if (nrow(df3) < 10) next
  
  plot_list[[i]] <-ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 0.5,color = "#c6c6c6")  +
    geom_point(data = df3,aes(x = X1, y = X2,fill = color),color ="black",shape = 21,stroke = 0.2, size = 3)+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(title = paste(m, "n =",nrow(df3), "cell lines"))+
    theme_classic() + theme(legend.position="none", plot.title = element_text(size=12)) + scale_fill_manual(values = color)
  
}

g <- grid.arrange(grobs=plot_list,ncol=7)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/UMAP_growth_mediums_individual_subplots.pdf',plot = g,device = 'pdf',width = 20, height = 10,dpi=300)


#media compiled onto one UMAP plot
p <- ggplot(df, aes(x = X1, y = X2)) +
  geom_point(size = 1,color = "#c6c6c6")  +
  geom_point(data = merged_df,aes(x = X1, y = X2,fill = FormulationID ),color ="black",shape = 21,stroke = 0.2, size = 3)+
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  labs(title = "Cell line growth media")+
  theme_classic() + theme(legend.position="right") + scale_fill_manual(values = new_palette)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/UMAP_optimization/UMAP_growth_mediums.pdf',plot = p,device = 'pdf',width = 20, height = 7,dpi=300)


#counting number of blood vs. solid tissue cancers cultured in each media
media_type_cancer_type.df <- merge(merged_df[,c(2,3,4)],CCLE_annotation[,c(4,30)], by = "StrippedCellLineName")
media_cancer_count.df <- data.frame()
for (m in unique(merged_df$FormulationID)) {
  
   cancer_types_media <- media_type_cancer_type.df[media_type_cancer_type.df$FormulationID == m,]
   blood_cancer_count <- sum(cancer_types_media$OncotreeLineage %in% c("Myeloid", "Lymphoid")) #counting number of blood cancer cell lines cultured in each cell growht m edia
   
   media_cancer_count.df <- rbind(media_cancer_count.df, data.frame(media = m, num_blood_cancer_cell_lines = blood_cancer_count, 
                                                                    total_cell_lines = nrow(cancer_types_media), 
                                                                    blood_versus_solid_cancer_cell_lines = blood_cancer_count/nrow(cancer_types_media)))
}


