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
library(dplyr)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/Cell line to cell line variability")

#loading metabolic pathway file
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"

#extracing metabolic pathways and associated genes
pathways <-  gmtPathways(pathway_file)

KEGG_metabolic_genes = as.data.frame(unique(unlist(pathways)))

colnames(KEGG_metabolic_genes) = c("Metabolic_gene")

#importing sanger gene expression data 
gene_exp <- read_tsv("/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/CellLinesProject_CompleteGeneExpression_v99_GRCh38.tsv")

metabolic_gene_exp <- gene_exp[,c(2,4,6)] %>%
  filter(GENE_SYMBOL %in% unique(KEGG_metabolic_genes$Metabolic_gene)) 

metabolic_gene_exp$SAMPLE_NAME <- toupper(gsub("-","", metabolic_gene_exp$SAMPLE_NAME))

#importing depmap metadata, assigning Sanger cell lines to cancer types based on DepMap data
DepMap_metadata <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/DepMap_23Q2_data/Model.csv")

sanger_depmap_metadata <- merge(metabolic_gene_exp, DepMap_metadata[,c(4,29)], by.x = "SAMPLE_NAME", by.y = "StrippedCellLineName")
sanger_depmap_metadata <- unique(sanger_depmap_metadata[,c(1,4)])


#filtering cancer types that have at least 5 cell lines
num_cell_lines_per_cancer <- sanger_depmap_metadata %>%
  group_by(OncotreePrimaryDisease) %>%
  dplyr::count() %>%
  filter(n >= 5)

sanger_depmap_metadata <- sanger_depmap_metadata[sanger_depmap_metadata$OncotreePrimaryDisease %in% num_cell_lines_per_cancer$OncotreePrimaryDisease,]


#organizing metabolic gene expression
metabolic_gene_exp <- metabolic_gene_exp[metabolic_gene_exp$SAMPLE_NAME %in% sanger_depmap_metadata$SAMPLE_NAME,]

metabolic_gene_exp.wide <- metabolic_gene_exp %>%
  unique() %>%
  pivot_wider(names_from = GENE_SYMBOL, values_from = Z_SCORE) %>%
  as.data.frame()

metabolic_gene_exp.wide <- metabolic_gene_exp.wide[order(metabolic_gene_exp.wide$SAMPLE_NAME),]
rownames(metabolic_gene_exp.wide) <- metabolic_gene_exp.wide$SAMPLE_NAME; metabolic_gene_exp.wide <- metabolic_gene_exp.wide[,-c(1)]

sanger_depmap_metadata <- sanger_depmap_metadata[order(sanger_depmap_metadata$SAMPLE_NAME),]

remove(DepMap_metadata,num_cell_lines_per_cancer)


#getting unique cancer types
metabolic_gene_exp.wide$cellType <- sanger_depmap_metadata$OncotreePrimaryDisease
metabolic_gene_exp.wide <- metabolic_gene_exp.wide[!(metabolic_gene_exp.wide$cellType == "Non-Cancerous"),]

celltypes <- unique(metabolic_gene_exp.wide$cellType)


#cancer cell line to cell line variability, code adapted from Xiao et al (2019, Nature Commn, Metabolic landscape of the tumor microenvironment at single cell resolution)
#changed from single-cell data to cell line data
enrich_data_df <- data.frame(CELLTYPE=NULL,PATHWAY=NULL,NES=NULL,PVAL=NULL, ADJPVAL = NULL)
pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),types=character())

for (t in celltypes){
  t2 <- str_replace(t," ","")
  each_metabolic_tpm <- metabolic_gene_exp.wide[metabolic_gene_exp.wide$cellType==t,]
  each_metabolic_tpm <- each_metabolic_tpm[ , -which(names(each_metabolic_tpm) %in% c("cellType"))]
  
  x <- as.matrix(each_metabolic_tpm)
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(x)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             types=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)

  
  #GSEA pre-rank analysis of absolute values of PC loadings
  set.seed(100) #setting seed for GSEA
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  colnames(pre_rank_matrix) <- c("GeneRank")
  pre_rank_vector <- as.numeric(pre_rank_matrix)
  names(pre_rank_vector) <- rownames(pre_rank_matrix)
  
  metabolic_pathways <-  gmtPathways(pathway_file)
  gseaRes <- fgseaSimple(pathways = metabolic_pathways, 
                    stats = pre_rank_vector,
                    nperm = 5000, minSize  = 10,
                    maxSize  = 500, scoreType = "pos")
  
  

  enrich_data_df <- rbind(enrich_data_df,data.frame(CELLTYPE=t2,PATHWAY=gseaRes$pathway,NES=gseaRes$NES,PVAL=gseaRes$pval, ADJPVAL = gseaRes$padj))
}

#export the result
fwrite(enrich_data_df,"/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/Cell line to cell line variability/pre_rank_GSEA_KEGG_metabolic_pathways_Sanger.csv")

##plot PCA variance for each cancer type
p <- ggplot(pc_plotdata) + geom_point(aes(x,100*y,colour=factor(sel)),size=0.75) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~factor(types),scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=12),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",linewidth =0.2),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        strip.background = element_rect(fill="white",linewidth=0.2,colour = NULL),
        strip.text=element_text(size=8))

ggsave(file.path(outDir,"PC_variance_pancancer.pdf"),p,width = 10,height=30,units="in",device="pdf",useDingbats=FALSE)

#evaluate adjusted pvalue <= 0.05 pathways only
min_pval <- by(enrich_data_df$ADJPVAL, enrich_data_df$PATHWAY, FUN=min)
select_pathways <- names(min_pval)[(min_pval<=0.05)]
select_enrich_data_df <- enrich_data_df[enrich_data_df$PATHWAY%in% select_pathways,]

#conver adjusted p-value to -log10
pvals <- select_enrich_data_df$ADJPVAL
pvals[pvals<=0] = 1e-10
select_enrich_data_df$ADJPVAL <- -log10(pvals)


#heatmap of -log10 adjusted p-values
color <- colorRampPalette(c("white","red"))(100)
dat <-  select_enrich_data_df[ , -which(names(select_enrich_data_df) %in% c("NES"))]
dat <- dat[,-c(3)]
dat <- reshape(dat, idvar = "PATHWAY", timevar = "CELLTYPE", direction = "wide")
rownames(dat)<- dat$PATHWAY
dat <- dat[,-c(1)]

cancer_character_names <- data.frame(celltype = character())
for (i in 1:ncol(dat)) {
  cancer_character_names[i,1] <- substr(colnames(dat)[i], 9, nchar(colnames(dat)[i]))
}

colnames(dat) <- cancer_character_names$celltype
test_labels <- matrix(as.numeric(unlist(dat)),nrow=nrow(dat))
test_labels[test_labels < -log10(0.05)] <- ""
test_labels[test_labels >= -log10(0.05)] <- "*"

p <- pheatmap(dat,cluster_cols = T,cluster_rows = T,color=color, display_numbers = test_labels, fontsize_number = 15, number_color = "black")
ggsave(file.path(outDir,"pathway_variability_across_cancers.pdf"),p,width = 9, height=6,units="in",device="pdf",useDingbats=FALSE)

#exporting variable cancer types
OXPHOS_variable_cancers <- select_enrich_data_df[select_enrich_data_df$ADJPVAL >= -log10(0.05) & select_enrich_data_df$PATHWAY == "Oxidative phosphorylation",]
fwrite(OXPHOS_variable_cancers,"/Volumes/FallahiLab/Maize-Data/Data/Cara/Sanger/Cell line to cell line variability/Sanger_OXPHOS_variable_cancertypes.csv")


#determining percentage of cancers that exhibit significant cell line-cell line heterogeneity across metabolic pathways
pathway_counts <- rowCounts(dat >= -log10(0.05))
pathway_counts.df <- as.data.frame(pathway_counts)
pathway_counts.df$pathway_name <- rownames(dat)
pathway_counts.df <- pathway_counts.df[order(pathway_counts.df$pathway_counts, decreasing = FALSE),]
pathway_counts.df$pathway_counts <- (pathway_counts.df$pathway_counts/length(celltypes))*100
pathway_counts.df$pathway_name <- factor(pathway_counts.df$pathway_name, levels = as.character(pathway_counts.df$pathway_name))



#bar plot of cancers that exhibit significant cell line-cell line heterogeneity across metabolic pathways
p<-ggplot(pathway_counts.df, aes(x = pathway_counts, y=pathway_name))  + 
geom_bar(stat="identity") +theme_classic()  + xlim(0,50) + theme(plot.title = element_text(size = 10),
axis.text.x=element_text(colour="black", size = 8),
axis.text.y=element_text(colour="black", size = 8),
axis.line=element_line(linewidth=0.2,color="black"),
axis.ticks = element_line(colour = "black",linewidth=0.2),
panel.border = element_blank(), panel.background = element_blank(),
axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 10, colour = "black")) + ylab("") + 
xlab("% of cancers with cell-cell pathway variability") + ggtitle("Pathway variability across Pancancer (Sanger)") 

ggsave(file.path(outDir,"Percent_cancers_pathway_variability.pdf"),p,width = 8,height=6,units="in",device="pdf",useDingbats=FALSE)

