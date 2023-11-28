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

#importing single cell R data
load("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/GBM_PTEN_malignant_cells.RData")

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data")

#importing gene datasets

pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"
pathways <-  gmtPathways(pathway_file)
KEGG_OXPHOS_genes <- pathways[["Oxidative phosphorylation"]]

ComplexI_genes <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/Harmonized single-cell landscape of glioblastoma/GSEA files/RESPIRATORY_CHAIN_COMPLEX_I.v2023.1.Hs.tsv",sep="\t")
ComplexI_genes <- unlist(strsplit(ComplexI_genes[17,2], ","))

ComplexII_genes <- c("SDHA", "SDHB", "SDHC", "SDHD") #GO:0045273

ComplexIII_genes <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/Harmonized single-cell landscape of glioblastoma/GSEA files/GOCC_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_III.v2023.1.Hs.tsv",sep="\t")
ComplexIII_genes <- unlist(strsplit(ComplexIII_genes[17,2], ","))

ComplexIV_genes <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/Harmonized single-cell landscape of glioblastoma/GSEA files/GOCC_RESPIRATORY_CHAIN_COMPLEX_IV.v2023.1.Hs.tsv",sep="\t")
ComplexIV_genes <- unlist(strsplit(ComplexIV_genes[17,2], ","))

ATPsynthase_genes <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/Harmonized single-cell landscape of glioblastoma/GSEA files/GOCC_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX.v2023.1.Hs.tsv",sep="\t")
ATPsynthase_genes <- unlist(strsplit(ATPsynthase_genes[17,2], ","))

Cytochrome_complex_genes <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CellxGene single cell data/Harmonized single-cell landscape of glioblastoma/GSEA files/GOCC_CYTOCHROME_COMPLEX.v2023.1.Hs.tsv",sep="\t")
Cytochrome_complex_genes <- unlist(strsplit(Cytochrome_complex_genes[17,2], ","))

#combining all OXPHOS and mitochondria genes
genes_all <- unique(c(KEGG_OXPHOS_genes,ComplexI_genes,ComplexIII_genes,ComplexIV_genes,ATPsynthase_genes,Cytochrome_complex_genes, ComplexII_genes))

#extracting only mito/OXPHOS genes
single_cell_OXPHOS_mito_genes.df <- glio_PTEN_subset_gene_exp.df[,colnames(glio_PTEN_subset_gene_exp.df) %in% genes_all]
single_cell_OXPHOS_mito_genes.df <- single_cell_OXPHOS_mito_genes.df %>% arrange(desc(rownames(single_cell_OXPHOS_mito_genes.df)))

#reordering PTEN mutation status df
PTEN_mutation_status.df <- PTEN_mutation_status.df %>% arrange(desc(SampleID))

remove(glio_PTEN_subset_gene_exp.df)
gc()

#compiling gene set list
gene_set_list <- list(KEGG_OXPHOS_genes = KEGG_OXPHOS_genes,ComplexII_genes = ComplexII_genes,
                      ComplexI_genes = ComplexI_genes,ComplexIII_genes = ComplexIII_genes,
                      ComplexIV_genes = ComplexIV_genes,ATPsynthase_genes = ATPsynthase_genes,
                      Cytochrome_complex_genes = Cytochrome_complex_genes)


#calculating mean expression for each gene set across all single cells
mean_gene_set_exp_per_sample.df <- data.frame()
for (s in names(gene_set_list)) {

  genes <- gene_set_list[[s]]
  df <- single_cell_OXPHOS_mito_genes.df[,colnames(single_cell_OXPHOS_mito_genes.df) %in% genes]
  rowmeans_sample <- as.numeric(rowMeans(df))

  mean_gene_set_exp_per_sample.df <- rbind(mean_gene_set_exp_per_sample.df,
                                           data.frame(Mean_exp = rowmeans_sample, SampleID = rownames(df), gene_set = s))

}

#merging mean gene expression df with mutation status df
mean_gene_set_exp_per_sample.df <- merge(mean_gene_set_exp_per_sample.df, PTEN_mutation_status.df,
                                         by= "SampleID")

#calculating median and IQR
IQR_spread <- mean_gene_set_exp_per_sample.df %>% group_by(gene_set, PTEN_status) %>%
  summarize(IQR = quantile(Mean_exp, probs = 0.75) - quantile(Mean_exp, probs = 0.25),
            median = quantile(Mean_exp, probs = 0.5),
            gene_set, PTEN_status)
IQR_spread <- unique(IQR_spread)


#merging IQR and gene expression data frames
mean_gene_set_exp_per_sample.df <- merge(mean_gene_set_exp_per_sample.df,IQR_spread, by = c("gene_set", "PTEN_status"))


#split violin function
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


#violin of mean gene expression in PTEN mutated vs non-mutated malignant cells
p <- ggplot(mean_gene_set_exp_per_sample.df, aes(x = gene_set, y=Mean_exp, fill = PTEN_status)) + facet_wrap(~gene_set, scales = "free") +
  geom_split_violin() + theme_classic() + xlab(NULL) + stat_summary(fun.y = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.75, linewidth = 0.3, position = position_dodge(width = .75), alpha = 0.7, show.legend = FALSE) + 
  ylab("mean gene expression") + scale_fill_manual(values=c("tomato", "skyblue1")) +
  theme(legend.position="right",panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=10),
        axis.line=element_line(linewidth=0.5,colour="black"),
        axis.ticks = element_line(colour = "black",linewidth=0.5),
        axis.text.x = element_blank(),
        axis.text.y=element_text(colour="black", size = 10),
        strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
        strip.text=element_text(size=8), text = element_text(size = 10),
        strip.text.x = element_text(size = 12)) +
  geom_text(x = 0.8, y = 1.1, size = 3, aes(label = paste0("IQR = ", signif(IQR,digits = 2))),
            data = unique(mean_gene_set_exp_per_sample.df[mean_gene_set_exp_per_sample.df$PTEN_status == "Mutated",c(1,2,6)])) +
  geom_text(x = 1.2, y = 0.9, size = 3, aes(label = paste0("IQR = ", signif(IQR,digits = 2))),
            data = unique(mean_gene_set_exp_per_sample.df[mean_gene_set_exp_per_sample.df$PTEN_status == "Not Mutated",c(1,2,6)]))

ggsave(file.path(outDir, 'violin_mean_gene_set_exp_PTEN_mutant_malignant_cells.pdf'), plot = p,device = 'pdf',width =8, height = 10,dpi=300)


#exporting mean gene expression scores for Matlab permutation analysis
fwrite(mean_gene_set_exp_per_sample.df,file.path(outDir,"mean_mito_OXPHOS_single_cell_exp_dataset.csv"))
