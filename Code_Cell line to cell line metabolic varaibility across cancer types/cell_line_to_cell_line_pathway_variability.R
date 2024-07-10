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

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pancancer_pathway_variability")

#loading metabolic pathway file
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"


#loading z-scored TPM data 
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/CCLE_z_score_data.RData")

#extracing metabolic pathways and associated genes
pathways <-  gmtPathways(pathway_file)

KEGG_metabolic_genes = as.data.frame(unique(unlist(pathways)))

colnames(KEGG_metabolic_genes) = c("Metabolic_gene")


#extracting cell lines with z-scored metabolic genes
CCLE_TPM_metabolic_genes <- as.data.frame(t(CCLE_TPM_zscore[,colnames(CCLE_TPM_zscore) %in% KEGG_metabolic_genes$Metabolic_gene]))

CCLE_TPM_metabolic_genes <- as.data.frame(t(CCLE_TPM_metabolic_genes))

#getting unique cancer types
CCLE_TPM_metabolic_genes$cellType <- CCLE_cell_line_tumor_type$OncotreePrimaryDisease

celltypes <- unique(CCLE_TPM_metabolic_genes$cellType)


#cancer cell line to cell line variability, code adapted from Xiao et al (2019, Nature Commn, Metabolic landscape of the tumor microenvironment at single cell resolution)
#changed from single-cell data to cell line data
enrich_data_df <- data.frame(CELLTYPE=NULL,PATHWAY=NULL,NES=NULL,PVAL=NULL, ADJPVAL = NULL)
pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),types=character())

for (t in celltypes){
  t2 <- str_replace(t," ","")
  each_metabolic_tpm <- CCLE_TPM_metabolic_genes[CCLE_TPM_metabolic_genes$cellType==t,]
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
fwrite(enrich_data_df,"/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pancancer_pathway_variability/pre_rank_GSEA_KEGG_metabolic_pathways.csv")

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
ggsave(file.path(outDir,"pathway_variability_across_cancers.pdf"),p,width = 18, height=16,units="in",device="pdf",useDingbats=FALSE)



#determining percentage of cancers that exhibit significant cell line-cell line heterogeneity across metabolic pathways
pathway_counts <- rowCounts(dat >= -log10(0.05))
pathway_counts.df <- as.data.frame(pathway_counts)
pathway_counts.df$pathway_name <- rownames(dat)
pathway_counts.df <- pathway_counts.df[order(pathway_counts.df$pathway_counts, decreasing = FALSE),]
pathway_counts.df$pathway_counts <- (pathway_counts.df$pathway_counts/41)*100
pathway_counts.df$pathway_name <- factor(pathway_counts.df$pathway_name, levels = as.character(pathway_counts.df$pathway_name))



#bar plot of cancers that exhibit significant cell line-cell line heterogeneity across metabolic pathways
p<-ggplot(pathway_counts.df, aes(x = pathway_counts, y=pathway_name))  + 
geom_bar(stat="identity") +theme_classic()  + xlim(0,40) + theme(plot.title = element_text(size = 10),
axis.text.x=element_text(colour="black", size = 8),
axis.text.y=element_text(colour="black", size = 8),
axis.line=element_line(linewidth=0.2,color="black"),
axis.ticks = element_line(colour = "black",linewidth=0.2),
panel.border = element_blank(), panel.background = element_blank(),
axis.ticks.length= unit(.5, "mm"), axis.title = element_text(size = 10, colour = "black")) + ylab("") + 
xlab("% of cancers with cell-cell pathway variability") + ggtitle("Pathway variability across Pancancer") 

ggsave(file.path(outDir,"Percent_cancers_pathway_variability.pdf"),p,width = 8,height=6,units="in",device="pdf",useDingbats=FALSE)


#determining spread of metabolic pathways (TPM values) across cancers
cancer_types <- unique(celltypes)
all_cell_types <- as.vector(CCLE_cell_line_tumor_type$OncotreePrimaryDisease)

TPM_non_normalized <- (2^CCLE_TPM)-1 #undoing log2(TPM+1) to just get raw TPM values


mean_metabolic_TPM_cancers <- data.frame()
IQR_cancer_types_final.df <- data.frame()

cancer_types <- unique(celltypes)
all_cell_types <- as.vector(CCLE_cell_line_tumor_type$OncotreePrimaryDisease)

for (p in names(pathways)) {
  genes <- pathways[[p]]
  
  #calculating spread of metabolic pathways across cell lines for metabolic pathways that have at least 3 genes in pathway
  if (length(genes) > 3)
  pathway_metabolic_TPM <- TPM_non_normalized[,colnames(TPM_non_normalized) %in% genes]
  pathway_metabolic_TPM_mean <- rowMeans(pathway_metabolic_TPM) #taking mean TPM across all pathway genes
  pathway_metabolic_TPM_statistics<- as.data.frame(pathway_metabolic_TPM_mean)
  colnames(pathway_metabolic_TPM_statistics) <- c("mean_TPM")
  pathway_metabolic_TPM_statistics$Cancer_Type <- CCLE_cell_line_tumor_type$OncotreePrimaryDisease
  
  
  #calculating IQR of data spread across cell lines each cancer type for each pathway
  IQR_cancer_types<- apply(pathway_metabolic_TPM_statistics, 2, function(x) by(pathway_metabolic_TPM_statistics$mean_TPM, pathway_metabolic_TPM_statistics$Cancer_Type, quantile))

  
  IQR_cancer_types.df <- data.frame(Cancer_Type = NA, Pathway = NA, IQR = NA)
  for (i in cancer_types) {
    IQR_cancer_types.df$Cancer_Type <- i
    IQR_cancer_types.df$Pathway <-names(pathways[p])
    IQR_cancer_types.df$IQR <- IQR_cancer_types[["mean_TPM"]][[i]][["75%"]] - IQR_cancer_types[["mean_TPM"]][[i]][["25%"]]

    
    IQR_cancer_types_final.df <- rbind(IQR_cancer_types_final.df,IQR_cancer_types.df)
  }
  
  
  #assembling final metabolic score data frame
  pathway_metabolic_TPM_statistics$Pathway <- names(pathways[p])
  pathway_metabolic_TPM_statistics$cell_line <- rownames(pathway_metabolic_TPM)
  rownames(pathway_metabolic_TPM_statistics) <- NULL
  
  mean_metabolic_TPM_cancers <- rbind(mean_metabolic_TPM_cancers,pathway_metabolic_TPM_statistics)
  
}




#mean pathway TPM values across cancer cell lines and different metabolic pathways
for (i in c("Ascorbate and aldarate metabolism","Thiamine metabolism","Oxidative phosphorylation")) {
  state.df <- mean_metabolic_TPM_cancers[mean_metabolic_TPM_cancers$Pathway == i, ]
  IQR.df <- IQR_cancer_types_final.df[IQR_cancer_types_final.df$Pathway == i,]
  
  merged.df <- merge(state.df, IQR.df, by = "Cancer_Type")

  
  #keeping same order of cancer types based on IQR descending order for ascorbate and aldarate metabolism
  if (i == c("Ascorbate and aldarate metabolism")) {
    merged.df <- merged.df[order(-merged.df$IQR),] #ordering IQR values
    
    ascorbate_order_cell_lines = merged.df$cell_line
    ascorbate_order_cancer_types_unique = unique(merged.df$Cancer_Type)
  }
  
  merged.df <- merged.df[match(ascorbate_order_cell_lines, merged.df$cell_line), ]  
  merged.df$Cancer_Type <- factor(merged.df$Cancer_Type, levels = ascorbate_order_cancer_types_unique)
  
  
  g <- ggplot(merged.df,aes(x=Cancer_Type,y=mean_TPM)) +
    geom_boxplot(fill = "lightblue", size=0.4,show.legend = F, outlier.size = 1) + labs(y=NULL,x=NULL) + 
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 12,angle=45,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 12),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",linewidth=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))+ ylab("mean metabolic pathway TPM") + ylim(-1,round(max(state.df$mean_TPM))) + labs(title = i) + 
    annotate("text", x=1:length(unique(state.df$Cancer_Type)), y=-1, label= paste(round(unique(merged.df$IQR), 2)), angle='30', color = "red") +
          ylim(-2,round_any(max(state.df$mean_TPM), 10) + 10)

  ggsave(file.path(outDir, (sprintf('%s_boxplot_state_scores.pdf',i))), plot = g,device = 'pdf',width =18, height = 12,dpi=300)
}






#calculating z-scored metabolic state scores for further analysis
pathway_metabolic_zscore_mean_cancers <- data.frame()
for (p in pathway_counts.df$pathway_name) {
  genes <- pathways[[p]]
  
  
  if (length(genes) > 10)
  pathway_metabolic_zscore <- CCLE_TPM_zscore[,colnames(CCLE_TPM_zscore) %in% genes]
  pathway_metabolic_zscore_mean <- as.data.frame(rowMeans(pathway_metabolic_zscore)) #taking mean z-scores
  colnames(pathway_metabolic_zscore_mean) <- c("mean_zscore")
  pathway_metabolic_zscore_mean$Pathway <- names(pathways[p])
  pathway_metabolic_zscore_mean$cell_line <- rownames(pathway_metabolic_zscore_mean)
  pathway_metabolic_zscore_mean$Cancer_Type <- CCLE_cell_line_tumor_type$OncotreePrimaryDisease
  rownames(pathway_metabolic_zscore_mean) <- NULL
  
  pathway_metabolic_zscore_mean_cancers <- rbind(pathway_metabolic_zscore_mean_cancers,pathway_metabolic_zscore_mean)
}


#filtering out cell lines in which their cancer types have significant variability (adjusted p <= 0.05)
variable_pathways = pathway_counts.df$pathway_name
cell_line_metabolic_variability.df <- data.frame()

for (p in variable_pathways) {
  variable_metabolic_pathway_cancer_types <- dat[rownames(dat) == c(p),] >= -log10(0.05)
  metabolic_column_idx <- which(apply(variable_metabolic_pathway_cancer_types, 2, any)) #keeping only cancer types that have significant variables (i.e. TRUE statements)
  metabolic_column_names <- names(metabolic_column_idx)
  metabolic_column_names <- gsub(" ", "", metabolic_column_names)
  
  metabolic_CCLE_cancer_type_variability <-  pathway_metabolic_zscore_mean_cancers[pathway_metabolic_zscore_mean_cancers$Pathway == c(p),]
  metabolic_CCLE_cancer_type_variability$Cancer_Type <- gsub(" ", "", metabolic_CCLE_cancer_type_variability$Cancer_Type)
  
  metabolic_CCLE_cancer_type_variability <-metabolic_CCLE_cancer_type_variability[metabolic_CCLE_cancer_type_variability$Cancer_Type %in% metabolic_column_names,]
  
  cell_line_metabolic_variability.df <- rbind(cell_line_metabolic_variability.df,metabolic_CCLE_cancer_type_variability)
}


#saving state scores for other analysis
save(cell_line_metabolic_variability.df, file = "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Pancancer_pathway_variability/cell_line_to_cell_line_metabolic_variability_scores.RData")

