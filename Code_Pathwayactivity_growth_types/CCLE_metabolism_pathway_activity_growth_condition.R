#following code was developed by Xiao et al (2019, Nature Commn, Metabolic landscape of the tumor microenvironment at single cell resolution)
#changed from single-cell data to growth conditions across cell lines


rm(list=ls()) #clear all
cat("\014") #clc

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



outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity")
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"

#loading z-scored TPM data and CCLE data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")

#getting metabolic pathway names
pathways <-  gmtPathways(pathway_file)
pathway_names <- names(pathways)

KEGG_metabolic_genes = as.data.frame(unique(unlist(pathways)))
colnames(KEGG_metabolic_genes) = c("Metabolic_gene")

#extracting metabolic genes from KEGG data base for pathway activity analysis
CCLE_TPM_metabolic_genes <- as.data.frame(t(CCLE_TPM[,colnames(CCLE_TPM) %in% KEGG_metabolic_genes$Metabolic_gene]))

#some genes occur in multiple pathways.
gene_pathway_number <- num_of_pathways(pathway_file,rownames(CCLE_TPM_metabolic_genes))

#filtering growth conditions to only ones with at least 10 cell lines cultured in 
growth_conditions <- CCLE_annotation[,c(4,11,28)]; growth_conditions <- growth_conditions[growth_conditions$StrippedCellLineName %in% colnames(CCLE_TPM_metabolic_genes),]
count_growth_conditions <- as.data.frame(table(growth_conditions$GrowthPattern)); growth_conditions <- growth_conditions[!(growth_conditions$GrowthPattern == "Unknown"),]


#extracting cell lines based on available growth conditions
norm_tpm <- CCLE_TPM_metabolic_genes
norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% growth_conditions$StrippedCellLineName] #filtering data to only cell lines that we have culture data for
norm_tpm <- norm_tpm[,order(colnames(norm_tpm))]

#extracing growth conditions across cell lines
growth_conditions <- growth_conditions[order(growth_conditions$StrippedCellLineName),]
all_growth_conditions <- as.vector(growth_conditions$GrowthPattern)
unique_growth_conditions <- unique(all_growth_conditions)

##Calculate the pathway activities 
#mean ratio of genes in each pathway for each growth condition
set.seed(123)
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(unique_growth_conditions),dimnames = list(pathway_names,unique_growth_conditions))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(unique_growth_conditions),dimnames = list(pathway_names,unique_growth_conditions))

###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(unique_growth_conditions),dimnames = (list(pathway_names, unique_growth_conditions)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 10) next

  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachGrowthType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_growth_conditions, mean))

  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachGrowthType)[colAlls(mean_exp_eachGrowthType>0.001)]


  if(length(keep)<5) next
  
  #using the lowest value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))

  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachGrowthType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_growth_conditions, mean))
  ratio_exp_eachGrowthType <- t(mean_exp_eachGrowthType) / colMeans(mean_exp_eachGrowthType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachGrowthType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachGrowthType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 5) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachGrowthType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_growth_conditions, mean))
  ratio_exp_eachGrowthType <- t(mean_exp_eachGrowthType) / colMeans(mean_exp_eachGrowthType)
  mean_exp_pathway <- apply(ratio_exp_eachGrowthType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[unique_growth_conditions]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[unique_growth_conditions]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(unique_growth_conditions,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_growth_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachGrowthType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_growth_types_list <- lapply(times,function(x) sample(all_growth_conditions)) 
  names(shuffle_growth_types_list) <- times
  mean_exp_eachGrowthType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachGrowthType_list <- lapply(times,function(x) mean_exp_eachGrowthType_list[[x]] / rowMeans(mean_exp_eachGrowthType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(unique_growth_conditions),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- unique_growth_conditions
  for(c in unique_growth_conditions){
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


#heatmap of pathway activity scores
dat <- mean_expression_shuffle

dat[is.na(dat)] <- 1

mybreaks <- c(
  seq(0.6, .97, length.out=33),
  seq(0.98, 1.02, length.out=1),
  seq(1.03, 1.5,length.out=33))

color <- colorRampPalette(c("blue","white","red"))(67)

#adding number of cell lines for each growth condition to label
count_growth_conditions <- count_growth_conditions[!(count_growth_conditions$Var1 == "Unknown"),]
count_growth_conditions <- count_growth_conditions[order(factor(count_growth_conditions$Var1, levels = c("Adherent", "Mixed", "Suspension", "Organoid"))), ]
count_growth_conditions$labels <- paste(count_growth_conditions$Var1, "(n=",count_growth_conditions$Freq,")")

colnames(dat) <- count_growth_conditions$labels

p <- pheatmap(dat,cluster_cols = T,cluster_rows = T,color=color,breaks=mybreaks,vannotation_row = NA, annotation_colors = mat_colors)
ggsave(file.path(outDir, "KEGGpathway_activity_heatmap_growth_conditions.pdf"),p,width=10,height=15,units="in",device="pdf",useDingbats=FALSE)




#exporting results from pathway activity calculations
write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle_growth_condition.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle_growth_condition.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue_growth_condition.txt"),row.names=T,col.names=T,quote=F,sep="\t")


#using PCA to determine the most variable metabolic pathways across all cancer types
PCA_pathway <- prcomp(t(dat))

#calculating data variance via PCA
eigs <- PCA_pathway$sdev^2
data_variance_explained <- 100*(eigs / sum(eigs))
cum_var <- cumsum(data_variance_explained)

#plotting PCA data explained variance
data_variance_explained.df <- as.data.frame(data_variance_explained[1:2])
colnames(data_variance_explained.df) <- c("Data_Variance_Explained")
data_variance_explained.df$Variance <- c("Each_PC_Variance")

cum_var.df <- as.data.frame(cum_var[1:2])
colnames(cum_var.df) <- c("Data_Variance_Explained")
cum_var.df$Variance <- c("Cummulative_Variance")
variance_bind.df <- rbind(data_variance_explained.df,cum_var.df)
variance_bind.df$PC <- c(1:2)


p <- ggplot(variance_bind.df,aes(x = PC, y = Data_Variance_Explained, fill = Variance)) + geom_bar(stat = "identity",position = position_dodge(), width = 0.7) +  xlab("Prinicpal Components") +
  ylab("Data Variance (%)") +
  ggtitle("Variance across PCs")+
  theme_classic() + theme(axis.text.x=element_text(colour="black", size = 12), axis.text.y=element_text(colour="black", size = 12), axis.title = element_text(size = 13, colour = "black"), title = element_text(size = 15, colour = "black")) + ylim(0,100)

ggsave(file.path(outDir,"PCA_variance_explained_growth_conditions.pdf"),p,width = 8,height=6,units="in",device="pdf",useDingbats=FALSE)


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

#plotting PCA loadings (top 10 and bottom 10 values)
p<-ggplot(pathway_loading_score_condensed, aes(x=Loading, y=pathway)) +
  geom_bar(stat="identity")+theme_classic() +  theme(legend.position="none",
                                                     axis.text.x=element_text(colour="black", size = 12),
                                                     axis.text.y=element_text(colour="black", size = 12),
                                                     axis.line=element_line(linewidth=0.2,color="black"),
                                                     axis.ticks = element_line(colour = "black",size=0.2),
                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                     axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 15, colour = "black")) + ylab("") + xlab('Summed loadings |PC1-PC2|') + xlim(0,0.5)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity/summed_PC1_PC2_variance_growth_condition.pdf',plot = p,width = 12,height=12,units="in",device="pdf",useDingbats=FALSE)













