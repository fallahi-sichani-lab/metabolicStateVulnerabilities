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
library(overlapping)
library(lattice)
library(gridExtra)

#importing single cell data and metadata and DepMap cell line info
metadata <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Tirosh_pancancer_cell_line_heterogeneity/Metadata.txt")
metadata <- metadata[-c(1),]
metadata$Cell_line <- gsub("\\_.*","",metadata$Cell_line)

CPM_data <- readRDS("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Tirosh_pancancer_cell_line_heterogeneity/CPM_data.rds")
rownames(CPM_data) <- CPM_data$GENE; CPM_data <- CPM_data[,-c(1)]
colnames(CPM_data) <- gsub("\\.", "-", colnames(CPM_data))

cell_line_data <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Tirosh_pancancer_cell_line_heterogeneity/OXPHOS_variable_cell_lines_mutation_KO_data.csv")

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Tirosh_pancancer_cell_line_heterogeneity")

#importing gene datasets
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"
pathways <-  gmtPathways(pathway_file)
KEGG_OXPHOS_genes <- pathways[["Oxidative phosphorylation"]]


#extracting cancer types and cell lines that are OXPHOS high and OXPHOS low
OXPHOS_high_low_cell_lines <- cell_line_data %>%
  filter(OXPHOS_state %in% c("OXPHOS_high", "OXPHOS_low"))

OXPHOS_cell_lines_single_cell <- merge(OXPHOS_high_low_cell_lines[,c(1:4),], metadata[,c(1:4)], by.x = "StrippedCellLineName", by.y = "Cell_line")

OXPHOS_cell_lines_single_cell <- OXPHOS_cell_lines_single_cell %>%
  dplyr::group_by(Cancer_Type) %>%
  filter(all(c("OXPHOS_high", "OXPHOS_low") %in% OXPHOS_state))
  

#creating histograms for OXPHOS variable cancer types
for (c in unique(OXPHOS_cell_lines_single_cell$Cancer_Type)) {
  
  df <- OXPHOS_cell_lines_single_cell[OXPHOS_cell_lines_single_cell$Cancer_Type == c,]
  df <- df[order(-df$mean_zscore),]
  
  cell_lines <- unique(df$StrippedCellLineName)
  
  plot_list <- list()
  for (t in cell_lines) {
    cancer_type_mean_OXPHOS.df <- data.frame()
    
    single_cells <- unique(OXPHOS_cell_lines_single_cell[OXPHOS_cell_lines_single_cell$StrippedCellLineName == t,]$NAME)
    single_cell_genes <- CPM_data[,colnames(CPM_data) %in% single_cells]
    single_cell_OXPHOS_genes <- single_cell_genes[rownames(single_cell_genes) %in% KEGG_OXPHOS_genes,]
    mean_OXPHOS <- colMeans(single_cell_OXPHOS_genes)
    
    #percent OXPHOS high/low
    percent_high = 100*sum(as.numeric(mean_OXPHOS) > 200)/length(mean_OXPHOS)
    percent_low= 100*sum(as.numeric(mean_OXPHOS) < 200)/length(mean_OXPHOS)
    
    cancer_type_mean_OXPHOS.df <- rbind(cancer_type_mean_OXPHOS.df, data.frame(Cell_line = t, OXPHOS_expression = mean_OXPHOS, 
                                                                             OXPHOS_state = ifelse(unique(OXPHOS_cell_lines_single_cell[OXPHOS_cell_lines_single_cell$StrippedCellLineName == t,]$mean_zscore) > 0, "OXPHOS_high", "OXPHOS_low"),
                                                                             OXPHOS_score = unique(OXPHOS_cell_lines_single_cell[OXPHOS_cell_lines_single_cell$StrippedCellLineName == t,]$mean_zscore)))
  

  plot_list[[t]] <- ggplot(cancer_type_mean_OXPHOS.df, aes(x = OXPHOS_expression)) + 
                    ggtitle(unique(cancer_type_mean_OXPHOS.df$Cell_line)) +
                    geom_histogram(fill = "lightgray", aes(y = after_stat((count)/sum(count))), 
                                    binwidth = 5, color = "black", linewidth = 0.3) +
                    scale_y_continuous(labels = scales::percent, limits = c(0, 0.15)) + 
                    theme_classic() + theme(legend.position="none",panel.grid.major = element_blank(),
                                            panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=26),
                                            axis.line=element_line(linewidth=0.5,colour="black"),
                                            axis.ticks = element_line(colour = "black",linewidth=0.5),
                                            axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                                            axis.text.y=element_text(colour="black", size = 12),
                                            strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                                            strip.text=element_text(size=12), text = element_text(size = 12)) + 
                    xlim(50,450) + xlab(NULL) + ylab(NULL) + geom_vline(xintercept = 200, colour="red", linetype = "dashed", linewidth = 0.8) +
                    annotate("text", x=300, y=0.12, label= paste("% cells OXPHOS high =", round(percent_high, digits = 2)), col="red", size=3) +
                    annotate("text", x=120, y=0.12, label= paste("% cells OXPHOS low =", round(percent_low, digits =2)), col="blue", size=3)
  
  }
 
  p <- grid.arrange(grobs=plot_list,nrow=length(cell_lines), bottom = "single cell mean OXPHOS expression (CPM)", left = "% single cells")
  
  ggsave(file.path(outDir, sprintf('%s_cell_lines_OXPHOS_single_cell_hist.pdf',c)), plot = p,device = 'pdf',width =5, height = 1.8*length(cell_lines),dpi=300)
  
  
}


