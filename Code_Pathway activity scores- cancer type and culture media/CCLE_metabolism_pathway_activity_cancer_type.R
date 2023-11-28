#following pathway activity score code was developed by Xiao et al (2019, Nature Commn, Metabolic landscape of the tumor microenvironment at single cell resolution)
#changed from single-cell data to cell lines across different cancer types

rm(list=ls()) #clear all
cat("\014") #clc

##libraries
library(stringr)
library(reshape2)
library(scales)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(matrixStats)
library(gridExtra)
library(randomcoloR)
library(fgsea)

#functions required for pathway activity analysis
#calculate how many pathways of one gene involved function (developed by Xiao et al)
num_of_pathways <- function (gmtfile,overlapgenes){
  pathways <- gmtPathways(gmtfile)
  pathway_names <- names(pathways)
  filter_pathways <- list()
  for (p in pathway_names){
    genes <- pathways[[p]]
    common_genes <- intersect(genes,overlapgenes)
    if(length(common_genes>=5)){
      filter_pathways[[p]] <- common_genes
    }
  }
  
  all_genes <- unique(as.vector(unlist(filter_pathways)))
  gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
  for(p in pathway_names){
    for(g in filter_pathways[[p]]){
      gene_times[g,"num"] = gene_times[g,"num"]+1
    }
  }
  gene_times
} 


#output directory 
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity")

#metabolism KEGG gmt file upload
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"

#loading z-scored TPM data and CCLE data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")



#getting metabolic pathway names
pathways <-  gmtPathways(pathway_file)
pathway_names <- names(pathways)

#getting all cancer cell line type
all_cell_types <- as.vector(CCLE_cell_line_tumor_type$OncotreePrimaryDisease)
cell_types <- unique(all_cell_types)

KEGG_metabolic_genes = as.data.frame(unique(unlist(pathways)))
colnames(KEGG_metabolic_genes) = c("Metabolic_gene")

#extracting metabolic genes from KEGG data base for pathway activity analysis
CCLE_TPM_metabolic_genes <- as.data.frame(t(CCLE_TPM[,colnames(CCLE_TPM) %in% KEGG_metabolic_genes$Metabolic_gene]))

#some genes occur in multiple pathways, identifying number of pathways a given gene is involved
gene_pathway_number <- num_of_pathways(pathway_file,rownames(CCLE_TPM_metabolic_genes))

#using log2(TPM + 1) values provided by CCLE
norm_tpm <- CCLE_TPM_metabolic_genes

##Calculate the pathway activities 
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))

#calculate the pvalues using shuffle/permutation
set.seed(123)
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 10) next

  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))

  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]


  if(length(keep)<5) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))

  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 5) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- (sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) + 1) / (5000 +1)
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- (sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) + 1) / (5000 + 1)
    }
    if(pval>0.05) mean_expression_shuffle[p, c] <- NA  ### if p-value is greater than 0.05, NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}

#removing pathways in which all values are NA
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]


#heatmap of permutated pathway activity scores
dat <- mean_expression_shuffle

dat[is.na(dat)] <- 1

mybreaks <- c(
  seq(0.6, .97, length.out=33),
  seq(0.98, 1.02, length.out=1),
  seq(1.03, 1.4,length.out=33))

color <- colorRampPalette(c("blue","white","red"))(67)


p <- pheatmap(dat,cluster_cols = T,cluster_rows = T,color=color,breaks=mybreaks,vannotation_row = NA, annotation_colors = mat_colors)
ggsave(file.path(outDir, "KEGGpathway_activity_heatmap.pdf"),p,width=16,height=15,units="in",device="pdf",useDingbats=FALSE)

#exporting results from pathway activity calculations
write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")



#boxplot show the distribution of pathway activity across cancer types (non-permutated pathway activity scores)
RNA_dat <- as.data.frame(mean_expression_noshuffle)
RNA_dat$X <- NULL

RNA_df <- melt(RNA_dat)
RNA_df <- RNA_df[!is.na(RNA_df$value),]
g <- ggplot(RNA_df,aes(x=reorder(variable,-value),y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,3),breaks=0:3,labels=0:3)+
  geom_boxplot(size=0.4,show.legend = F, outlier.size = 1)  + labs(y=NULL,x=NULL) + 
  scale_fill_manual(values = palette_final) +
  theme_classic() + ylim(0.5,2) + 
  theme(legend.position="none",
  axis.text.x=element_text(colour="black", size = 12,angle=45,hjust=1,vjust=1),
  axis.text.y=element_text(colour="black", size = 12),
  axis.line=element_line(size=0.2,color="black"),
  axis.ticks = element_line(colour = "black",size=0.2),
  panel.border = element_blank(), panel.background = element_blank(),
  axis.ticks.length= unit(.5, "mm")) + geom_hline(yintercept=1, linetype="dashed", color = "black") + ylab("Pathway Activity Score")

ggsave(file.path(outDir,"pathway_activity_boxplot.pdf"),g,width = 15,height=8,units="in",device="pdf",useDingbats=FALSE)



#using PCA to determine the most variable metabolic pathways across all cancer types
PCA_pathway <- prcomp(t(dat))

#calculating data variance via PCA
eigs <- PCA_pathway$sdev^2
data_variance_explained <- 100*(eigs / sum(eigs))
cum_var <- cumsum(data_variance_explained)

#plotting PCA data explained variance
data_variance_explained.df <- as.data.frame(data_variance_explained[1:8])
colnames(data_variance_explained.df) <- c("Data_Variance_Explained")
data_variance_explained.df$Variance <- c("Each_PC_Variance")

cum_var.df <- as.data.frame(cum_var[1:8])
colnames(cum_var.df) <- c("Data_Variance_Explained")
cum_var.df$Variance <- c("Cummulative_Variance")
variance_bind.df <- rbind(data_variance_explained.df,cum_var.df)
variance_bind.df$PC <- c(1:8)


p <- ggplot(variance_bind.df,aes(x = PC, y = Data_Variance_Explained, fill = Variance)) + geom_bar(stat = "identity",position = position_dodge(), width = 0.7) +  xlab("Prinicpal Components") +
  ylab("Data Variance (%)") +
  ggtitle("Variance across PCs")+
  theme_classic() + theme(axis.text.x=element_text(colour="black", size = 12), axis.text.y=element_text(colour="black", size = 12), axis.title = element_text(size = 13, colour = "black"), title = element_text(size = 15, colour = "black")) + ylim(0,100)

ggsave(file.path(outDir,"PCA_variance_explained.pdf"),p,width = 8,height=6,units="in",device="pdf",useDingbats=FALSE)


#calculating loadings across all metbaolic pathway across PCs that add up ~80% variance explained
select_pcs <- which(cum_var>=80)[1] #selecting loadings with 80% variance explained
loadings <- abs(PCA_pathway$rotation)
pathway_loading_score <- as.data.frame(sort(apply(abs(loadings[,1:select_pcs]), 1, sum)))

pathway_loading_score_condensed <- as.data.frame(pathway_loading_score)
pathway_loadings_condensed_rownames <- rownames(pathway_loading_score)

rownames(pathway_loading_score_condensed) <- pathway_loadings_condensed_rownames
pathway_loading_score_condensed$pathway <- pathway_loadings_condensed_rownames
colnames(pathway_loading_score_condensed)[1] <- c("Loading")

# lock in factor level order
pathway_loading_score_condensed$pathway <- factor(pathway_loading_score_condensed$pathway, levels = pathway_loading_score_condensed$pathway)

#plotting summed PCA loadings across all pathways
p<-ggplot(pathway_loading_score_condensed, aes(x=Loading, y=pathway)) +
  geom_bar(stat="identity")+theme_classic() +  theme(legend.position="none",
                                                     axis.text.x=element_text(colour="black", size = 12),
                                                     axis.text.y=element_text(colour="black", size = 12),
                                                     axis.line=element_line(size=0.2,color="black"),
                                                     axis.ticks = element_line(colour = "black",size=0.2),
                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                     axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 15, colour = "black")) + ylab("") + xlab('Summed loadings |PC1-PC8|') + xlim(0,3)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity/summed_PC1_PC8_variance.pdf',plot = p,width = 12,height=12,units="in",device="pdf",useDingbats=FALSE)


