rm(list=ls()) #clear all
cat("\014") #clc

#importing libraries needed
library(viridis)
library(stringr)
library(plyr)
library(dplyr)
library(gtools)
library(matrixStats)
library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(tidyr)
library(devtools)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(randomcoloR)

#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis")

#loading data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE_mutation_data.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/organized_CRISPR_gene_effect.RData")
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CRISPR_depMap/Elastic_net_OXPHOS/cell_line_to_cell_line_metabolic_variability_scores.RData")

#loading significant pancancer KO targets
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/PLSR/OXPHOS/significant_OXPHOS_pan-cancer_targets.RData")


#including cell lines in mutation data that have knockout data
CRISPR_cell_lines <- rownames(CRISPR_CCLE_data)
CCLE_all_mutations.df = CCLE_mutation_cell_line[CCLE_mutation_cell_line$StrippedCellLineName %in% CRISPR_cell_lines,]

#top 25 common mutations (based on Mendiratta et al (2021, Nature Comm.))
top_cancer_mutations <- c("TP53", "PIK3CA", "LRP1B", "KRAS", "APC", "FAT4", "KMT2D", "KMT2C", "BRAF", "ARID1A", "FAT1", "PTEN", "ATM", "ZFHX3", "CREBBP",
                          "GRIN2A", "NF1", "PDE4DIP", "PTPRT", "TRRAP", "RNF213", "PREX2", "SPEN", "ERBB4", "KMT2A")

#combining gene and mutation type into one column
mutations_CCLE_top_pancancer_mutations.df <- CCLE_all_mutations.df[CCLE_all_mutations.df$HugoSymbol %in% top_cancer_mutations,]

#keeping samples in which mutations have functional consequence (i.e., likely GoF or LoF)
functional_mutations.df <- mutations_CCLE_top_pancancer_mutations.df[!(mutations_CCLE_top_pancancer_mutations.df$LikelyGoF == "False" & mutations_CCLE_top_pancancer_mutations.df$LikelyLoF == "False"),]

#mutation type frequency count across all cell lines with KO data
mutation_count <- data.frame(table(unique(functional_mutations.df[,c(4,11)])$HugoSymbol))


#plotting mutation types and frequency across cell lines
p <- ggplot(mutation_count, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(axis.title.x=element_text(colour="black", size = 14),
        axis.title.y=element_text(colour="black", size = 14),
        axis.text.y=element_text(colour="black", size = 12),
        axis.text.x=element_text(colour="black", size = 12, angle = 90)) + ylab("mutation frequency") + xlab(NULL) + 
        labs(title =  "Driver mutation frequency across DepMap cell lines (top 25 mutations)") + geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)
ggsave(file.path(outDir, 'driver_mutation_frequency.pdf'), plot = p,device = 'pdf',width =12, height = 8,dpi=300)


#heatmap of mutation/cancer type pair across all DepMap cell lines with KO data
cancer_type_mutation_count <- unique(functional_mutations.df[,c(4,11,12)]) %>% dplyr::group_by(OncotreePrimaryDisease, HugoSymbol) %>% dplyr::summarise(count = n())
cancer_type_mutation_count_wide <- as.data.frame(spread(cancer_type_mutation_count, key = HugoSymbol, value = count))
rownames(cancer_type_mutation_count_wide) <- cancer_type_mutation_count_wide$OncotreePrimaryDisease; cancer_type_mutation_count_wide <- cancer_type_mutation_count_wide[,-c(1)]
cancer_type_mutation_count_wide[is.na(cancer_type_mutation_count_wide)] <- 0

rownames(mutation_count) <- mutation_count$Var1

cancer_type_count <- unique(functional_mutations.df[,c(11,12)]) %>% dplyr::group_by(OncotreePrimaryDisease) %>% dplyr::summarise(count = n())

color <- colorRampPalette(c("white","blue"))(10)

p <- pheatmap(cancer_type_mutation_count_wide,cluster_cols = F,cluster_rows = F,color=color, display_numbers = cancer_type_mutation_count_wide, fontsize_number = 15, number_color = "black", 
              labels_col =paste0(mutation_count$Var1, " (n = ", mutation_count$Freq," cell lines)"),
              labels_row = paste0(cancer_type_count$OncotreePrimaryDisease, " (n = ", cancer_type_count$count," cell lines)"))
ggsave(file.path(outDir, 'number_cell_line_mutation_pair_heatmap_across_all_cell_lines.pdf'), plot = p,device = 'pdf',width =20, height = 20,dpi=300)


####
####

#filtering mutation data to only OXPHOS variable cancer types/cell lines 
#cleaning up CRISPR dataset
cell_lines_CRISPR <- data.frame(rownames(CRISPR_CCLE_data),CRISPR_CCLE_data)
rownames(cell_lines_CRISPR) <- NULL
colnames(cell_lines_CRISPR)[1] <- "cell_line"

#removing cell lines with significant missing data
genes_NA <- as.data.frame(colSums(is.na(cell_lines_CRISPR[,c(2:17932)]))) 
genes_NA$empty <- 0
genes_NA <-genes_NA[genes_NA$`colSums(is.na(cell_lines_CRISPR[, c(2:17932)]))` == 0,]
genes_CRISPR_removed_NA<- cell_lines_CRISPR[,c("cell_line", rownames(genes_NA))]


#filtering CRISPR data to only selected pathways' cell lines
metabolic_pathways <- "Oxidative phosphorylation"

Gene_dependency_scores_metabolic_pathways <- data.frame()
for (p in metabolic_pathways) {
  metabolic_pathway_cell_lines <- cell_line_metabolic_variability.df[cell_line_metabolic_variability.df$Pathway == p,]
  metabolic_score_Gene_dependency <- merge(metabolic_pathway_cell_lines,genes_CRISPR_removed_NA, by = "cell_line") 
  
  Gene_dependency_scores_metabolic_pathways <- rbind(Gene_dependency_scores_metabolic_pathways, metabolic_score_Gene_dependency)
}


#filtering datasets to remove samples with metabolic state scores close to 0 based on 33% and 67% distribution thresholds
Gene_dependency_scores_metabolic_pathways_filtered <- data.frame()
for (p in metabolic_pathways) {
  df <- Gene_dependency_scores_metabolic_pathways[Gene_dependency_scores_metabolic_pathways$Pathway == p, ]
  df_quants <- as.data.frame(t(quantile(df$mean_zscore, probs = c(0.33,0.67))))
  
  df_filtered <- df[df$mean_zscore < df_quants$`33%` | df$mean_zscore > df_quants$`67%`,]
  
  
  Gene_dependency_scores_metabolic_pathways_filtered <- rbind(Gene_dependency_scores_metabolic_pathways_filtered, df_filtered)
}


#combining mutation data with OXPHOS variable cell lines
functional_mutations_all.df <- CCLE_all_mutations.df[!(CCLE_all_mutations.df$LikelyGoF == "False" & CCLE_all_mutations.df$LikelyLoF == "False"),]

mutation_OXPHOS_variable_cell_lines <- unique(merge(functional_mutations_all.df[c(4,11)], Gene_dependency_scores_metabolic_pathways[,c("cell_line","Cancer_Type", "mean_zscore",as.character(unique(sig_candidates$Gene)))], by.x = "StrippedCellLineName", by.y = "cell_line"))


#classifying if cell lines are OXPHOS high, low or neither (i.e., 0) based on 33/67 percentile cutoffs defined above
mutation_OXPHOS_variable_cell_lines <- mutation_OXPHOS_variable_cell_lines %>% mutate(
  OXPHOS_state = ifelse(mean_zscore < df_quants$`33%`, "OXPHOS_low", 
                       ifelse(mean_zscore > df_quants$`67%`, "OXPHOS_high", 0)
  )
)

mutation_OXPHOS_variable_cell_lines <- mutation_OXPHOS_variable_cell_lines[,c(1:4,32, 5:31)]


#cell lines that differ between mutation and cancer type datasets
cell_lines_diff <- setdiff(Gene_dependency_scores_metabolic_pathways$cell_line, unique(mutation_OXPHOS_variable_cell_lines$StrippedCellLineName))


#filtering mutations to common mutations found in OXPHOS variable cell lines
frequent_mutations <- as.character(mutation_count[mutation_count$Freq >= 10, ]$Var1)
mutation_OXPHOS_variable_cell_lines_with_frequent_mutations <- mutation_OXPHOS_variable_cell_lines[mutation_OXPHOS_variable_cell_lines$HugoSymbol %in% frequent_mutations,]


#number of cell lines for each mutation
mutation_count_OXPHOS_variable_cell_lines <- data.frame(table(mutation_OXPHOS_variable_cell_lines_with_frequent_mutations[!mutation_OXPHOS_variable_cell_lines_with_frequent_mutations$OXPHOS_state == 0 & 
                                                                                                                            mutation_OXPHOS_variable_cell_lines_with_frequent_mutations$HugoSymbol %in% frequent_mutations,]$HugoSymbol))

p <- ggplot(mutation_count_OXPHOS_variable_cell_lines, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(axis.title.x=element_text(colour="black", size = 14),
        axis.title.y=element_text(colour="black", size = 14),
        axis.text.y=element_text(colour="black", size = 12),
        axis.text.x=element_text(colour="black", size = 12, angle = 90)) + ylab("mutation frequency") + xlab(NULL) + 
  labs(title =  "Driver mutation frequency across OXPHOS high and low cell lines") + geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)
ggsave(file.path(outDir, 'driver_mutation_frequency_OXPHOS_high_low_cell_lines.pdf'), plot = p,device = 'pdf',width =12, height = 8,dpi=300)


#filtering to mutations that have at last 15 cell lines (based on plot above)
mutation_candidates <- as.character(mutation_count_OXPHOS_variable_cell_lines[mutation_count_OXPHOS_variable_cell_lines$Freq >= 15, ]$Var1)

#converting hugosymbol data into columns with binary numbers
mutation_present_each_cell_line_binary <- data.frame(table(mutation_OXPHOS_variable_cell_lines[,c(1,2)]))
mutation_present_each_cell_line_binary <- spread(mutation_present_each_cell_line_binary, key = HugoSymbol, value = Freq)
mutation_present_each_cell_line_binary <- mutation_present_each_cell_line_binary[, names(mutation_present_each_cell_line_binary) %in% c(mutation_candidates, "StrippedCellLineName")]

#merging binary dataframe with mutation_count_OXPHOS_variable_cell_lines data.frame
mutation_OXPHOS_variable_cell_lines <- unique(merge(mutation_OXPHOS_variable_cell_lines[,-c(2)], mutation_present_each_cell_line_binary, by = "StrippedCellLineName"))


#adding in cell lines with missing mutation status
missing_cell_lines <- Gene_dependency_scores_metabolic_pathways[Gene_dependency_scores_metabolic_pathways$cell_line %in% cell_lines_diff ,c("cell_line","mean_zscore",as.character(unique(sig_candidates$Gene)))]

missing_cell_lines <- missing_cell_lines %>% mutate(
  OXPHOS_state = ifelse(mean_zscore < df_quants$`33%`, "OXPHOS_low", 
                        ifelse(mean_zscore > df_quants$`67%`, "OXPHOS_high", 0)
  )
)

missing_cell_lines <- merge(missing_cell_lines, unique(cell_line_metabolic_variability.df[,c(3,4)]), by = "cell_line")
missing_cell_lines <- missing_cell_lines[,c(1,31,2,30, 3:29)]
colnames(missing_cell_lines)[1] <- "StrippedCellLineName"

missing_cell_lines[names(mutation_OXPHOS_variable_cell_lines)[32:42]] <- NA


#adding cell lines with missing mutation data into dataframe
mutation_OXPHOS_variable_cell_lines_final.df <- rbind(mutation_OXPHOS_variable_cell_lines,missing_cell_lines)


#saving Rdata  of OXPHOS variable cancer types' cell lines OXPHOS scores, OXPHOS state category, KO effects for pan-cancer targets, and mutation status 
save(mutation_OXPHOS_variable_cell_lines_final.df, file = "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Mutation analysis/driver_mutations_OXPHOS_variable_cell_lines.RData")

#saving csv of OXPHOS variable cancer types' cell lines OXPHOS scores, OXPHOS state category, KO effects for pan-cancer targets, and mutation status 
fwrite(mutation_OXPHOS_variable_cell_lines_final.df,file.path(outDir,"OXPHOS_variable_cell_lines_mutation_KO_data.csv"))

###
###

#heatmap of mutation/cancer type pair results (no filtering of middle population with score close to 0)
mutations_OXPHOS_variable_cancers <- functional_mutations.df[functional_mutations.df$StrippedCellLineName %in% unique(Gene_dependency_scores_metabolic_pathways$cell_line) & 
                                                               functional_mutations.df$HugoSymbol %in% 
                                                               c("TP53", "PIK3CA", "LRP1B", "KRAS", "APC", "FAT4", "KMT2D", "KMT2C", "BRAF", "ARID1A", "FAT1", "PTEN", "ATM", "ZFHX3", "CREBBP",
                                                                 "GRIN2A", "NF1", "PDE4DIP", "PTPRT", "TRRAP", "RNF213", "PREX2", "SPEN", "ERBB4", "KMT2A"),]

mutations_OXPHOS_variable_cancers_count <- ddply(unique(mutations_OXPHOS_variable_cancers[,c(4,11,12)]), .(unique(mutations_OXPHOS_variable_cancers[,c(4,11,12)])$HugoSymbol, unique(mutations_OXPHOS_variable_cancers[,c(4,11,12)])$OncotreePrimaryDisease), nrow)

names(mutations_OXPHOS_variable_cancers_count) <- c("HugoSymbol", "OncotreePrimaryDisease", "Freq")

num_cell_lines_per_cancer <- data.frame(table(unique(CCLE_all_mutations.df[CCLE_all_mutations.df$StrippedCellLineName %in% unique(Gene_dependency_scores_metabolic_pathways$cell_line),c("StrippedCellLineName", "OncotreePrimaryDisease")])$OncotreePrimaryDisease))

num_cell_lines_per_mutation <- aggregate(mutations_OXPHOS_variable_cancers_count$Freq, by=list(mutations_OXPHOS_variable_cancers_count$HugoSymbol), FUN=sum)

color <- colorRampPalette(c("white","lightblue"))(20)
dat <- dcast(mutations_OXPHOS_variable_cancers_count, OncotreePrimaryDisease ~ HugoSymbol, value.var="Freq")
rownames(dat) <- dat$OncotreePrimaryDisease
dat <- dat[,-c(1)]
dat[is.na(dat)] <- 0


p <- pheatmap(dat,cluster_cols = F,cluster_rows = F,color=color, display_numbers = dat, fontsize_number = 15, number_color = "black", 
              labels_row=paste0(num_cell_lines_per_cancer$Var1, " (", num_cell_lines_per_cancer$Freq," cell lines)"),
              labels_col = paste0(num_cell_lines_per_mutation$Group.1, " (", num_cell_lines_per_mutation$x," cell lines)"))
ggsave(file.path(outDir, 'number_cell_line_mutation_pair_heatmap_across_OXPHOS_variable_cell_lines_no_cell_line_filtration.pdf'), plot = p,device = 'pdf',width =13, height = 10,dpi=300)

#barplot of mutation frequency across OXPHOS variable cell lines
p <- ggplot(num_cell_lines_per_mutation, aes(reorder(Group.1, -x), x)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(axis.title.x=element_text(colour="black", size = 14),
        axis.title.y=element_text(colour="black", size = 14),
        axis.text.y=element_text(colour="black", size = 12),
        axis.text.x=element_text(colour="black", size = 12, angle = 90)) + ylab("mutation frequency") + xlab(NULL) + 
  labs(title =  "Driver mutation frequency across DepMap cell lines") + geom_text(aes(label=x), position=position_dodge(width=0.9), vjust=-0.25)
ggsave(file.path(outDir, 'driver_mutation_frequency_OXPHOS_variable_cell_lines_no_cell_line_filtration.pdf'), plot = p,device = 'pdf',width =12, height = 8,dpi=300)




#heatmap of mutation/cancer type pair for OXPHOS high and low cell lines
OXPHOS_high_low_cell_lines <- na.omit(mutation_OXPHOS_variable_cell_lines_final.df[!mutation_OXPHOS_variable_cell_lines_final.df$OXPHOS_state == 0, ])
cancer_type_mutation_count <- OXPHOS_high_low_cell_lines[,c(2,32:42)] %>% dplyr::group_by(Cancer_Type) %>% dplyr::summarise(across(everything(), sum))
rownames(cancer_type_mutation_count) <- cancer_type_mutation_count$Cancer_Type

cancer_type_count <- data.frame(table(OXPHOS_high_low_cell_lines$Cancer_Type))
mutation_count <- data.frame(colSums(cancer_type_mutation_count[,c(2:12)]))

color <- colorRampPalette(c("white","lightblue"))(20)

p <- pheatmap(cancer_type_mutation_count[,-c(1)],cluster_cols = F,cluster_rows = F,color=color, display_numbers = cancer_type_mutation_count[,-c(1)], fontsize_number = 15, number_color = "black", 
              labels_row=paste0(cancer_type_count$Var1, " (", cancer_type_count$Freq," cell lines)"),
              labels_col = paste0(rownames(mutation_count), " (", mutation_count$colSums.cancer_type_mutation_count...c.2.12...," cell lines)"))
ggsave(file.path(outDir, 'number_cell_line_mutation_pair_heatmap_across_OXPHOS_high_low_cell_lines.pdf'), plot = p,device = 'pdf',width =12, height = 10,dpi=300)