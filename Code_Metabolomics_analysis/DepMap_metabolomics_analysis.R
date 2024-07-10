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
library(plyr)
library(tidyr)



#load metabolomics data and reference metadata for clean up
metabolomics.df <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/DepMap_23Q2_data/Metabolomics.csv", header=FALSE)
metabolomics.df[c(1),c(86)]<- "hydroxyproline"
metabolomics.df[c(1),] <- gsub("\\/.*","",metabolomics.df[c(1),])
colnames(metabolomics.df) <- metabolomics.df[c(1),] 
metabolomics.df <- metabolomics.df[-c(1),] 

#uploading metadata for naming of lipids/metabolites
ref_met.txt <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/metabolite_lipid_HMDB.txt")

ref_met.txt <- ref_met.txt %>%
  mutate(HMDB = case_when(str_detect(HMDB, "HMDB") ~ HMDB, 
.default = Standardized.name)) %>%
  mutate(HMDB = case_when(str_detect(HMDB, "-") ~ Input.name, 
                          HMDB == "CAR 18:1" ~ Input.name,
                          .default = HMDB)) %>%
  mutate(lipid_metabolite = case_when(str_detect(Input.name, ":") ~ "lipid", .default = "metabolite"))


colnames(metabolomics.df)[c(9:233)] <- ref_met.txt$HMDB[c(1:225)]
metabolomics.df <- metabolomics.df[,colnames(metabolomics.df) %in% c("cell_line_display_name", "lineage_2", ref_met.txt[ref_met.txt$lipid_metabolite == "metabolite",]$HMDB)]

#loading z-scored TPM data and CCLE data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")
remove(CCLE_TPM, medium_metadata_unfiltered)

#filtering data to same cell lines that were analyzed from RNA-seq
metabolomics.df <- metabolomics.df[metabolomics.df$cell_line_display_name %in% rownames(CCLE_TPM_zscore),]
rownames(metabolomics.df) <- metabolomics.df$cell_line_display_name
metabolomics.df <- as.data.frame(metabolomics.df[,-c(1:2)]); cell_lines <- rownames(metabolomics.df)
metabolomics.df <- apply(metabolomics.df,2,as.numeric)
rownames(metabolomics.df) <- cell_lines

#scale metabolimics data
metabolomics_zscored.df <- as.data.frame(scale(metabolomics.df))
metabolomics_zscored.df <- merge(CCLE_annotation[,c(4,29,30)], metabolomics_zscored.df, by.x = "StrippedCellLineName",by.y = 0)
metabolomics_zscored.df$OncotreeLineage <- ifelse(metabolomics_zscored.df$OncotreeLineage == "Lymphoid" | 
                                                    metabolomics_zscored.df$OncotreeLineage == "Myeloid", "Hematopoietic Cancers", "Solid Tissue Cancers")

###
###

#UMAP clustering of metabolomics expression
#Set up R UMAP, optimizing parameters
for (nn in c(20)){
  for (md in c(0.3)){
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
      
      umap_output <- umap(as.matrix(metabolomics_zscored.df[,-c(1:3)]),config = custom.config)
      
      
      #Collect Output
      output = list(umap = umap_output$layout,
                    data.zscore = metabolomics_zscored.df,
                    parameters = data.frame(nn=nn,md=md,ncomp=ncomp,metric=metric))
      
      
      df = output[["data.zscore"]]
      df$X1=output[["umap"]][,1]
      df$X2=output[["umap"]][,2]
      df$cancer_type = metabolomics_zscored.df$OncotreePrimaryDisease
      df$tissue_lineage = metabolomics_zscored.df$OncotreeLineage
      
      
      
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
      ggsave(sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/UMAP-nn=%0.0f-md=%0.0e-%s.pdf',nn,md,metric),plot = g,device = 'pdf',width = 13, height = 7,dpi=300)
      
    }
  }
}

#UMAP for blood vs. solid tissue cancers
g = plot_embedding_categorical('tissue_lineage',df)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/blood_solid_tissue_cancers_UMAP.pdf',plot = g,device = 'pdf',width = 9, height = 7,dpi=300)

#UMAP subplots for each cancer type
cancer_list = unique(df$cancer_type)

plot_list <- list() 
for (i in 1:length(cancer_list)){
  kk = cancer_list[i]
  color = palette_final[which(cancer_list == kk)]
  
  
  df2 = filter(df,cancer_type%in%kk)
  
  plot_list[[i]] <-ggplot(df, aes(x = X1, y = X2)) +
    geom_point(size = 0.1,color = "#c6c6c6")  +
    geom_point(data = df2,aes(x = X1, y = X2,fill = color),color ="black",shape = 21,stroke = 0.2, size = 2)+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    labs(title = kk)+
    scale_fill_manual(values = color)+
    theme_classic() + theme(legend.position="none", plot.title = element_text(size=12))
  
}
g <- grid.arrange(grobs=plot_list,ncol=9)
ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/UMAP_cancer_types.pdf',plot = g,device = 'pdf',width = 20, height = 13,dpi=300)

###
###

#enrichment analysis of variable and non-variable metabolites across cancer types
metabolomics.df <- merge(CCLE_annotation[,c(4,29,30)], as.data.frame(metabolomics.df), by.x = "StrippedCellLineName",by.y = 0)
num_cancer_types <- as.data.frame(table(metabolomics.df$OncotreePrimaryDisease))
num_cancer_types <- num_cancer_types[num_cancer_types$Freq >= 5,] #filtering cancer types that have at least 5 cell lines profile

metabolomics_filtered.df <- metabolomics.df[metabolomics.df$OncotreePrimaryDisease %in% num_cancer_types$Var1,]


median_metabolite_expression_per_cancer_type <- metabolomics_filtered.df[,-c(1,3)] %>%
  group_by(OncotreePrimaryDisease) %>%
  summarise_all(median) %>%
  pivot_longer(!OncotreePrimaryDisease, names_to = "HMDB", values_to = "median_value") %>%
  ungroup()

IQR.df <- data.frame()
for (v in unique(median_metabolite_expression_per_cancer_type$HMDB)) {
  df <- median_metabolite_expression_per_cancer_type[median_metabolite_expression_per_cancer_type$HMDB == v,]
  IQR_value <- as.numeric(quantile(df$median_value, 0.75)) - as.numeric(quantile(df$median_value, 0.25))
  
  IQR.df <- rbind(IQR.df, data.frame("HMDB" = v, "IQR" = IQR_value))
}

median_metabolite_expression_per_cancer_type <- merge(median_metabolite_expression_per_cancer_type, IQR.df, by.x = "HMDB", by.y = "HMDB")
median_metabolite_expression_per_cancer_type <- merge(median_metabolite_expression_per_cancer_type, ref_met.txt[,c(8,9)], by.x = "HMDB", by.y = "HMDB")

  
df <- unique(median_metabolite_expression_per_cancer_type[,c(1,4:5)])
p <- ggplot(df, aes(x = IQR)) +
  ggtitle("IQR of metabolite median expression across cancer types") +
  geom_histogram(fill = "lightgray",  color = "black", linewidth = 0.3) +
  theme_classic() + theme(legend.position="none",panel.grid.major = element_blank(),
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=26),
                          axis.line=element_line(linewidth=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                          axis.text.y=element_text(colour="black", size = 12),
                          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                          strip.text=element_text(size=12), text = element_text(size = 12)) +
  geom_vline(data = ddply(median_metabolite_expression_per_cancer_type, .(lipid_metabolite), summarize, x_line = as.numeric(quantile(IQR, 0.5))), 
             aes(xintercept = x_line), color = "red", linetype = "dashed")

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/IQR_histogram_metabolite.pdf',plot = p,device = 'pdf',width = 6, height = 5,dpi=300)



#finding top 50 percentile list of lipids or metabolites
filtered_metabolites.df <- data.frame()
for (v in c("metabolite")) {
  df <- median_metabolite_expression_per_cancer_type[median_metabolite_expression_per_cancer_type$lipid_metabolite == v,]
  threshold <- quantile(unique(df$IQR), c(0.5))
  variable_metabolite <- df[df$IQR >= as.numeric(threshold),]
  variable_metabolite$variable_or_non_variable <- "variable"
  non_variable_metabolite <- df[df$IQR < as.numeric(threshold),]
  non_variable_metabolite$variable_or_non_variable <- "not variable"
  
  df_filtered <- rbind(variable_metabolite, non_variable_metabolite)
  
  filtered_metabolites.df <- rbind(filtered_metabolites.df, df_filtered)
}

metabolites_varied_aross_cancer_types.df <- unique(filtered_metabolites.df[filtered_metabolites.df$lipid_metabolite == "metabolite",c(1,4,6)])
write.csv(metabolites_varied_aross_cancer_types.df, "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/varied_metabolites_across_cancer_types.csv")


#Metaboanalyst results 
pathway_activity_variability <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity/PCA_thresholds_summed_loadings_pathways.txt")
pathway_activity_variability <- pathway_activity_variability[,c(1,3)]
variable_metabolites_pathway_results <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/Metaboanalyst results/msea_ora_result_variable_metabolites_50_threshold.csv")

variable_metabolite_results.df <- merge(variable_metabolites_pathway_results, pathway_activity_variability, by.x = "X", by.y = "pathways")
colnames(variable_metabolite_results.df)[1] <- "pathways"
variable_metabolite_results.df$significant <- ifelse(variable_metabolite_results.df$Raw.p <= 0.05, "Yes", "No")
variable_metabolite_results.df$significant <- factor(variable_metabolite_results.df$significant, levels = c("Yes", "No"))


variable_metabolite_results.df$Enrichment.Ratio <- variable_metabolite_results.df$hits/variable_metabolite_results.df$expected
bar_plot.df <- variable_metabolite_results.df[variable_metabolite_results.df$Raw.p <= 0.05,]

p <- ggplot(bar_plot.df, aes(x= Enrichment.Ratio, y= reorder(pathways, Enrichment.Ratio))) +  
    geom_bar(stat="identity")+ theme_classic() +  theme(legend.position="none",
                                                       axis.text.x=element_text(colour="black", size = 10),
                                                       axis.text.y=element_text(colour="black", size = 10),
                                                       axis.line=element_line(size=0.2,color="black"),
                                                       axis.ticks = element_line(colour = "black",size=0.2),
                                                       panel.border = element_blank(), panel.background = element_blank(),
                                                       axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + 
    ylab("") + xlab('Enrichment Ratio') + ggtitle(label = "Significant variable pathways (Metaboanalyst)") 
  
ggsave(file.path('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/cancer_type_metabolite_variability_enrichment_barplot.pdf'), plot = p,device = 'pdf',width =8, height = 4,dpi=300)

#correlation of enrichment ratio and gene loadings (variable metabolites)
correlation_metabolite_gene_pathway <- cor.test(variable_metabolite_results.df$Enrichment.Ratio,variable_metabolite_results.df$X80.)

variable_metabolite_results.df$significant <- factor(variable_metabolite_results.df$significant, levels = c("No","Yes"))
variable_metabolite_results.df$pathways <- ifelse(variable_metabolite_results.df$Enrichment.Ratio < 3 & 
                                                    variable_metabolite_results.df$X80. < 0.6, "", variable_metabolite_results.df$pathways)
variable_metabolite_results.df$label <- ifelse(variable_metabolite_results.df$pathways == "", "no_label","label")

p <- ggplot(variable_metabolite_results.df, aes(x = Enrichment.Ratio, y = X80., label = pathways)) +
  ggtitle("Variable pathway correlation", subtitle = paste0(sprintf("r = %s", 
                                                            round(correlation_metabolite_gene_pathway[["estimate"]][["cor"]], digits = 2)),
                                                            sprintf(", p = %s", round(correlation_metabolite_gene_pathway[["p.value"]], digits = 3)))) +
  geom_point(aes(color = significant, size = label)) + xlab("Enrichment ratio (Metaboanalyst metabolite analysis)") + ylab("Summed loading |PC1-PC8| loadings") +
  theme_classic() + theme(legend.position="none",panel.grid.major = element_blank(),
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=26),
                          axis.line=element_line(linewidth=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x = element_text(colour="black", size = 12),
                          axis.text.y=element_text(colour="black", size = 12),
                          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                          strip.text=element_text(size=12), text = element_text(size = 12)) +
  scale_color_manual(values=c("black", "red")) + geom_text(size = 3, hjust=0.2,vjust=-1) + 
  scale_size_manual(values=c(5, 2)) +
  geom_smooth(method = glm, se = FALSE, color = "blue", linetype = "dashed")

ggsave(file.path('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/cancer_type_variability_metabolite_versus_gene_expression_pathway_analysis.pdf'), plot = p,device = 'pdf',width =7, height = 7,dpi=300)




#Metaboanalyst results (metabolite enrichemnt analysis)
non_variable_metabolites_pathway_results <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/Metaboanalyst results/msea_ora_result_non_variable_metabolites_50_threshold.csv")

non_variable_metabolite_results.df <- merge(non_variable_metabolites_pathway_results, pathway_activity_variability, by.x = "X", by.y = "pathways")
colnames(non_variable_metabolite_results.df)[1] <- "pathways"
non_variable_metabolite_results.df$significant <- ifelse(non_variable_metabolite_results.df$Raw.p <= 0.05, "Yes", "No")
non_variable_metabolite_results.df$significant <- factor(non_variable_metabolite_results.df$significant, levels = c("Yes", "No"))

non_variable_metabolite_results.df$Enrichment.Ratio <- non_variable_metabolite_results.df$hits/non_variable_metabolite_results.df$expected
bar_plot.df <- non_variable_metabolite_results.df[non_variable_metabolite_results.df$Raw.p <= 0.05,]

p <- ggplot(bar_plot.df, aes(x= Enrichment.Ratio, y= reorder(pathways, Enrichment.Ratio))) +  
  geom_bar(stat="identity")+ theme_classic() +  theme(legend.position="none",
                                                      axis.text.x=element_text(colour="black", size = 10),
                                                      axis.text.y=element_text(colour="black", size = 10),
                                                      axis.line=element_line(size=0.2,color="black"),
                                                      axis.ticks = element_line(colour = "black",size=0.2),
                                                      panel.border = element_blank(), panel.background = element_blank(),
                                                      axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 12, colour = "black")) + 
  ylab("") + xlab('Enrichment Ratio') + ggtitle(label = "Significant non-variable pathways (Metaboanalyst)")

ggsave(file.path('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Metabolomics analysis/cancer_type_metabolite_non_variable_enrichment_barplot.pdf'), plot = p,device = 'pdf',width =8, height = 4,dpi=300)


