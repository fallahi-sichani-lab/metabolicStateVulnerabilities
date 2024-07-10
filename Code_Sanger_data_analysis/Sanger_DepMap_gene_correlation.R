rm(list=ls()) #clear all
cat("\014") #clc


#importing libraries needed
library(stringr)
library(pheatmap)
library(gtools)
library(matrixStats)
library(fgsea)
library(ggplot2)
library(data.table)
library(reshape2)
library(randomcoloR)
library(plyr)
library(grid)
library(gridExtra)
library(tidyverse)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/Gene correlation")

#loading metabolic pathway file
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"

#extracting metabolic pathways and associated genes
pathways <-  gmtPathways(pathway_file)
KEGG_metabolic_genes = as.data.frame(unique(unlist(pathways)))
colnames(KEGG_metabolic_genes) = c("Metabolic_gene")

#importing sanger gene expression data 
Sanger_gene_exp <- read_tsv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/CellLinesProject_CompleteGeneExpression_v99_GRCh38.tsv")

Sanger_metabolic_gene_exp <- Sanger_gene_exp[,c(2,4,6)] %>%
  filter(GENE_SYMBOL %in% unique(KEGG_metabolic_genes$Metabolic_gene)) 

Sanger_metabolic_gene_exp$SAMPLE_NAME <- toupper(gsub("-","", Sanger_metabolic_gene_exp$SAMPLE_NAME))

#importing DepMap gene expression
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")
remove(medium_metadata_unfiltered, CCLE_TPM, CCLE_cell_line_tumor_type)

CCLE_TPM_zscore <- CCLE_TPM_zscore[,colnames(CCLE_TPM_zscore) %in% KEGG_metabolic_genes$Metabolic_gene] %>%
  as.data.frame()

CCLE_TPM_zscore$SAMPLE_NAME <- rownames(CCLE_TPM_zscore)

CCLE_TPM_zscore_long <- CCLE_TPM_zscore %>%
  pivot_longer(
  cols = colnames(CCLE_TPM_zscore)[1]:colnames(CCLE_TPM_zscore)[1620], 
  names_to = "GENE_SYMBOL",
  values_to = "DepMap_Z_SCORE",
)

#merging gene expression dataframes from DepMap and Sanger
CCLE_Sanger_gene_exp <- merge(CCLE_TPM_zscore_long, Sanger_metabolic_gene_exp, by.x = c("SAMPLE_NAME", "GENE_SYMBOL"),
                              by.y = c("SAMPLE_NAME", "GENE_SYMBOL"))

CCLE_Sanger_gene_exp_correlation <- CCLE_Sanger_gene_exp %>%
  group_by(GENE_SYMBOL) %>%
  summarise(correlation = cor.test(DepMap_Z_SCORE, Z_SCORE)$estimate,
            p_value = cor.test(DepMap_Z_SCORE, Z_SCORE)$p.value)


#histogram looking at correlation of each metabolic gene across same cell lines in Depmap and Sanger datasets

p <- ggplot(CCLE_Sanger_gene_exp_correlation, aes(x = correlation)) + 
  ggtitle("Per gene correlation coefficients (DepMap vs. Sanger expression)") +
  geom_histogram(fill = "lightgray",  color = "black", linewidth = 0.3) +
  theme_classic() + theme(legend.position="none",panel.grid.major = element_blank(),
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=26),
                          axis.line=element_line(linewidth=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                          axis.text.y=element_text(colour="black", size = 12),
                          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                          strip.text=element_text(size=12), text = element_text(size = 12)) + 
  xlab("Correlation coefficient") + ylab("# genes") + geom_vline(xintercept = 0.5, colour="red", linetype = "dashed", linewidth = 0.8) + 
  annotate("text", x=0.2, y=130, label= paste("# genes (r > 0.5) =", sum(CCLE_Sanger_gene_exp_correlation$correlation >= 0.5)), col="red", size=5) +
  annotate("text", x=0.2, y=100, label= paste("# genes (r < 0.5) =", sum(CCLE_Sanger_gene_exp_correlation$correlation < 0.5)), col="blue", size=5)

ggsave(file.path(outDir, 'histogram_DepMap_Sanger_gene_correlation.pdf'), plot = p,device = 'pdf',width =7, height = 6,dpi=300)


#determining what genes are not well correlated and what pathways the genes are enriched in
CCLE_Sanger_gene_exp_corr <- CCLE_Sanger_gene_exp_correlation[CCLE_Sanger_gene_exp_correlation$correlation < 0.5,]

genes_not_correlated_pathway_enrichment.df <- data.frame()
for (p in names(pathways)) {
  
  pathway_genes <- pathways[[p]]
  percent_genes_in_pathways <- sum(pathway_genes %in% CCLE_Sanger_gene_exp_corr$GENE_SYMBOL)/length(pathway_genes)
  
  genes_not_correlated_pathway_enrichment.df <- rbind(genes_not_correlated_pathway_enrichment.df, 
                                                      data.frame(Pathway = p, percent_enrich = percent_genes_in_pathways,
                                                                 num_genes_in_pathway = length(pathway_genes)))
}

genes_not_correlated_pathway_enrichment.df <- genes_not_correlated_pathway_enrichment.df[genes_not_correlated_pathway_enrichment.df$num_genes_in_pathway >= 10,]

p <- ggplot(genes_not_correlated_pathway_enrichment.df[genes_not_correlated_pathway_enrichment.df$percent_enrich >= 0.29,], aes(x = 100*percent_enrich, y= reorder(Pathway, percent_enrich))) + 
  geom_bar(stat="identity") +theme_classic() + theme(plot.title = element_text(size = 14),
                                                                   axis.text.x=element_text(colour="black", size = 14),
                                                                   axis.text.y=element_text(colour="black", size = 12),
                                                                   axis.line=element_line(linewidth=0.2,color="black"),
                                                                   axis.ticks = element_line(colour = "black",linewidth=0.2),
                                                                   panel.border = element_blank(), panel.background = element_blank(),
                                                                   axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 14, colour = "black")) + 
                                                ylab("") + xlab("% of genes in pathway with variability") + 
                                                ggtitle("Pathways with significant gene expression variability in Sanger and DepMap")

ggsave(file.path(outDir, 'barplot_percent_genes_varied_in_pathways.pdf'), plot = p,device = 'pdf',width =11, height = 5,dpi=300)

#difference in metabolic genes profiled in Sanger vs. DepMap
diff_genes <- setdiff(unique(CCLE_TPM_zscore_long$GENE_SYMBOL), unique(Sanger_metabolic_gene_exp$GENE_SYMBOL))

genes_not_in_Sanger_enrichment.df <- data.frame()
for (p in names(pathways)) {
  
  pathway_genes <- pathways[[p]]
  percent_genes_in_pathways <- sum(pathway_genes %in% diff_genes)/length(pathway_genes)
  
  genes_not_in_Sanger_enrichment.df <- rbind(genes_not_in_Sanger_enrichment.df, 
                                                      data.frame(Pathway = p, percent_enrich = percent_genes_in_pathways,
                                                                 num_genes_in_pathway = length(pathway_genes)))
}

genes_not_in_Sanger_enrichment.df <- genes_not_in_Sanger_enrichment.df[genes_not_in_Sanger_enrichment.df$num_genes_in_pathway >= 10,]

p <- ggplot(genes_not_in_Sanger_enrichment.df[genes_not_in_Sanger_enrichment.df$percent_enrich >= 0.29,], aes(x = 100*percent_enrich, y= reorder(Pathway, percent_enrich))) + 
  geom_bar(stat="identity") +theme_classic() + theme(plot.title = element_text(size = 14),
                                                     axis.text.x=element_text(colour="black", size = 14),
                                                     axis.text.y=element_text(colour="black", size = 12),
                                                     axis.line=element_line(linewidth=0.2,color="black"),
                                                     axis.ticks = element_line(colour = "black",linewidth=0.2),
                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                     axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 14, colour = "black")) + 
  ylab("") + xlab("% of genes in pathway that are missing") + 
  ggtitle("Pathways where Sanger data has missing genes")

ggsave(file.path(outDir, 'barplot_missing_genes_in_pathways_from_sanger.pdf'), plot = p,device = 'pdf',width =8, height = 4,dpi=300)
