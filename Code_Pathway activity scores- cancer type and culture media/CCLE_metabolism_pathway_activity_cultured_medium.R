#following pathway activity score code was developed by Xiao et al (2019, Nature Commn, Metabolic landscape of the tumor microenvironment at single cell resolution)
#changed from single-cell data to growth media across cell lines


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


#output directory and metabolism gmt file pathway
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

#filtering media to only ones with at least 10 cell lines cultured in 
count_medium = as.data.frame(table(medium_metadata_unfiltered$FormulationID))
count_medium = count_medium[count_medium$Freq >= 10,] #filering medium with at least 10 cell lines cultured in 

medium_metadata_unfiltered <- medium_metadata_unfiltered[medium_metadata_unfiltered$FormulationID %in% c("RPMI + 10% FBS", "DMEM + 10% FBS", "EMEM + 10% FBS + 2mM Glutamine + 100uM NEAA",
                                                                                                         "RPMI + 10% FBS + 2mM Glutamine", "DMEM:F12 + 10% FBS", "McCoy's 5A + 10% FBS", 
                                                                                                         "MEM + 10% FBS", "F12 + 10% FBS", "RPMI + 20% FBS", "EMEM + 10% FBS",
                                                                                                         "DMEM + 10% FBS + 2mM Glutamine", "IMDM + 20% FBS + 4mM Glutamine + 1x Insulin-Transferrin-Selenium",
                                                                                                         "OPAC:Ad+++ (1:1)", "IMDM + 10% FBS", "RPMI + 10% FBS + 2mM Glutamine + 25mM HEPES + 25mM Sodium bicarbonate",
                                                                                                         "RPMI + 5% FBS", "L-15 + 10% FBS", "DMEM:F12 + 5% FBS + 2mM Glutamine + 5ug/ml Insulin + 10ug/ml Transferrin + 30nM Selenium + 10nM Hydrocortisone + 10nM Beta estradiol",
                                                                                                         "F12 + 15% FBS + 2mM Glutamine", "RPMI + 15% FBS", "IMDM + 10% FBS + 2mM Glutamine"),]
medium_metadata_unfiltered <- unique(medium_metadata_unfiltered)
medium_metadata_unfiltered <- medium_metadata_unfiltered[!medium_metadata_unfiltered$StrippedCellLineName == "JIMT1",] #removing cell line that is repeated twice in DepMap dataset

#extracting cell lines based on available media
norm_tpm <- CCLE_TPM_metabolic_genes
norm_tpm <- norm_tpm[,colnames(norm_tpm) %in% medium_metadata_unfiltered$StrippedCellLineName] #filtering data to only cell lines that we have culture data for
norm_tpm <- norm_tpm[,order(colnames(norm_tpm))]

#extracting culture mediums across cell lines
medium_metadata_unfiltered <- medium_metadata_unfiltered[medium_metadata_unfiltered$StrippedCellLineName %in% colnames(norm_tpm),]
medium_metadata_unfiltered <- medium_metadata_unfiltered[order(medium_metadata_unfiltered$StrippedCellLineName),]
all_media_types <- as.vector(medium_metadata_unfiltered$FormulationID)
media_types <- unique(all_media_types)

##Calculate the pathway activities 
#mean ratio of genes in each pathway for each media type
set.seed(123)
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(media_types),dimnames = list(pathway_names,media_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(media_types),dimnames = list(pathway_names,media_types))

###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(media_types),dimnames = (list(pathway_names, media_types)))

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 10) next

  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachMediaType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_media_types, mean))

  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachMediaType)[colAlls(mean_exp_eachMediaType>0.001)]


  if(length(keep)<5) next
  
  #using the lowest value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))

  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachMediaType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_media_types, mean))
  ratio_exp_eachMediaType <- t(mean_exp_eachMediaType) / colMeans(mean_exp_eachMediaType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachMediaType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachMediaType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 5) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachMediaType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_media_types, mean))
  ratio_exp_eachMediaType <- t(mean_exp_eachMediaType) / colMeans(mean_exp_eachMediaType)
  mean_exp_pathway <- apply(ratio_exp_eachMediaType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[media_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[media_types]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(media_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_media_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachMediaType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_media_types_list <- lapply(times,function(x) sample(all_media_types)) 
  names(shuffle_media_types_list) <- times
  mean_exp_eachMediaType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachMediaType_list <- lapply(times,function(x) mean_exp_eachMediaType_list[[x]] / rowMeans(mean_exp_eachMediaType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(media_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- media_types
  for(c in media_types){
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
  seq(1.03, 1.8,length.out=33))

color <- colorRampPalette(c("blue","white","red"))(67)

p <- pheatmap(dat,cluster_cols = T,cluster_rows = T,color=color,breaks=mybreaks,vannotation_row = NA, annotation_colors = mat_colors)
ggsave(file.path(outDir, "KEGGpathway_activity_heatmap_cultured_media.pdf"),p,width=10,height=19,units="in",device="pdf",useDingbats=FALSE)


#exporting results from pathway activity calculations
write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle_cell_media.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle_cell_media.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue_cell_media.txt"),row.names=T,col.names=T,quote=F,sep="\t")


#using PCA to determine the most variable metabolic pathways across all media types
PCA_pathway <- prcomp(t(dat))

#calculating data variance via PCA
eigs <- PCA_pathway$sdev^2
data_variance_explained <- 100*(eigs / sum(eigs))
cum_var <- cumsum(data_variance_explained)

#plotting PCA data explained variance
data_variance_explained.df <- as.data.frame(data_variance_explained[1:4])
colnames(data_variance_explained.df) <- c("Data_Variance_Explained")
data_variance_explained.df$Variance <- c("Each_PC_Variance")

cum_var.df <- as.data.frame(cum_var[1:4])
colnames(cum_var.df) <- c("Data_Variance_Explained")
cum_var.df$Variance <- c("Cummulative_Variance")
variance_bind.df <- rbind(data_variance_explained.df,cum_var.df)
variance_bind.df$PC <- c(1:4)

#plotting data variance captured by the first 4 PCs
p <- ggplot(variance_bind.df,aes(x = PC, y = Data_Variance_Explained, fill = Variance)) + geom_bar(stat = "identity",position = position_dodge(), width = 0.7) +  xlab("Prinicpal Components") +
  ylab("Data Variance (%)") +
  ggtitle("Variance across PCs")+
  theme_classic() + theme(axis.text.x=element_text(colour="black", size = 12), axis.text.y=element_text(colour="black", size = 12), axis.title = element_text(size = 13, colour = "black"), title = element_text(size = 15, colour = "black")) + ylim(0,100)

ggsave(file.path(outDir,"PCA_variance_explained_culture_media.pdf"),p,width = 8,height=6,units="in",device="pdf",useDingbats=FALSE)


#calculating absolutes summed loadings across across PCs that add up ~80% variance explained for all metabolic pathways
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
                                                     axis.line=element_line(linewidth=0.2,color="black"),
                                                     axis.ticks = element_line(colour = "black",size=0.2),
                                                     panel.border = element_blank(), panel.background = element_blank(),
                                                     axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 15, colour = "black")) + ylab("") + xlab('Summed loadings |PC1-PC4|') + xlim(0,2)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pathway_activity/summed_PC1_PC4_variance_culture_media.pdf',plot = p,width = 12,height=12,units="in",device="pdf",useDingbats=FALSE)

#boxplots to visualize metabolic pathway heterogeneity across different media types
#determining spread of metabolic pathways (TPM values) across media types
media_types <- unique(all_media_types)

TPM_non_normalized <- t((2^norm_tpm)-1)#undoing log2(TPM+1) to just get raw TPM values

mean_metabolic_TPM_media <- data.frame()
IQR_media_types_final.df <- data.frame()


for (p in names(pathways)) {
  genes <- pathways[[p]]
  
  #calculating spread of metabolic pathways across cell lines for metabolic pathways that have at least 3 genes in pathway
  if (length(genes) > 3)
  pathway_metabolic_TPM <- TPM_non_normalized[,colnames(TPM_non_normalized) %in% genes]
  pathway_metabolic_TPM_mean <- rowMeans(pathway_metabolic_TPM) #taking mean TPM across all pathway genes
  pathway_metabolic_TPM_statistics<- as.data.frame(pathway_metabolic_TPM_mean)
  colnames(pathway_metabolic_TPM_statistics) <- c("mean_TPM")
  pathway_metabolic_TPM_statistics$Media <- medium_metadata_unfiltered$FormulationID
  
  #calculating IQR of data spread across each media type  for each pathway
  IQR_media_types<- apply(pathway_metabolic_TPM_statistics, 2, function(x) by(pathway_metabolic_TPM_statistics$mean_TPM, pathway_metabolic_TPM_statistics$Media, quantile))
  
  IQR_media_types.df <- data.frame(Media = NA, Pathway = NA, IQR = NA)
  for (i in media_types) {
    IQR_media_types.df$Media <- i
    IQR_media_types.df$Pathway <-names(pathways[p])
    IQR_media_types.df$IQR <- IQR_media_types[["mean_TPM"]][[i]][["75%"]] - IQR_media_types[["mean_TPM"]][[i]][["25%"]]
    
    IQR_media_types_final.df <- rbind(IQR_media_types_final.df,IQR_media_types.df)
  }
  
  
  #assembling final metabolic score data frame
  pathway_metabolic_TPM_statistics$Pathway <- names(pathways[p])
  pathway_metabolic_TPM_statistics$cell_line <- rownames(pathway_metabolic_TPM)
  rownames(pathway_metabolic_TPM_statistics) <- NULL
  
  mean_metabolic_TPM_media <- rbind(mean_metabolic_TPM_media,pathway_metabolic_TPM_statistics)
  
}


#mean pathway TPM values across cancer cell lines and different metabolic pathways
for (i in c("Ascorbate and aldarate metabolism","Thiamine metabolism","Oxidative phosphorylation")) {
  state.df <- mean_metabolic_TPM_media[mean_metabolic_TPM_media$Pathway == i, ]
  IQR.df <- IQR_media_types_final.df[IQR_media_types_final.df$Pathway == i,]
  
  merged.df <- merge(state.df, IQR.df, by = "Media")
  
  #keeping same order of media types based on IQR descending order for ascorbate and aldarate metabolism
  if (i == c("Ascorbate and aldarate metabolism")) {
    merged.df <- merged.df[order(-merged.df$IQR),] #ordering IQR values
    
    ascorbate_order_cell_lines = merged.df$cell_line
    ascorbate_order_media_types_unique = unique(merged.df$Media)
  }
  
  merged.df <- merged.df[match(ascorbate_order_cell_lines, merged.df$cell_line), ]  
  merged.df$Media <- factor(merged.df$Media, levels = ascorbate_order_media_types_unique)
  
  
  g <- ggplot(merged.df,aes(x=Media,y=mean_TPM)) +
    geom_boxplot(fill = "lightblue", size=0.4,show.legend = F, outlier.size = 1) + labs(y=NULL,x=NULL) + 
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 12,angle=45,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 12),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",linewidth=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))+ ylab("mean metabolic pathway TPM") + ylim(-1,round(max(state.df$mean_TPM)+10)) + labs(title = i) + 
    annotate("text", x=1:length(unique(state.df$Media)), y=-1, label= paste(round(unique(merged.df$IQR), 2)), angle='30', color = "red") +
    ylim(-2,round(max(state.df$mean_TPM), 10) + 10)
  
  ggsave(file.path(outDir, (sprintf('%s_boxplot_state_scores_across_media.pdf',i))), plot = g,device = 'pdf',width =18, height = 18,dpi=300)
}











