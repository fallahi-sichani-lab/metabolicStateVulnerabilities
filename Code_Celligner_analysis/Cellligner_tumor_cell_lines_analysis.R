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
library(umap)
library(dplyr)

#importing Celligner data
Celligner_data <- readRDS("/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/Celligner_aligned_data.rds")
rownames(Celligner_data) <- Celligner_data$X; Celligner_data <- Celligner_data[,-c(1)]

metadata <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/Celligner_info.csv")

#filtering to cancer types we wish to analyze and organizing data
filter_tumors_CL <- metadata %>%
  group_by(lineage, subtype, type) %>%
  count() 

metadata_filtered <- metadata %>%
  filter(subtype %in% c("uveal melanoma", "undifferentiated", "transitory", "neural crest-like",
                        "melanocytic", "Ewing sarcoma", "adenocarcinoma", "acute myeloid leukemia",
                        "melanoma", "colon adenocarcinoma", "ductal adenocarcinoma, exocrine",
                        "prostate adenocarcinoma", "glioblastoma", "glioblastoma multiforme",
                        "glioma", "osteosarcoma", "esophageal carcinoma", "squamous cell carcinoma",
                        "hepatocellular carcinoma", "lung adenocarcinoma") |
           lineage %in% c("breast") |
           grepl("non-small cell lung cancer", subtype) |
           grepl("cystadenocarcinoma", subtype) |
           grepl("renal cell carcinoma", subtype) |
           grepl("acute lymphoblastic leukemia", subtype) |
           lineage %in% c("thyroid") & grepl("carcinoma", subtype) |
           lineage %in% c("soft_tissue") & grepl("synovial", subtype) |
           lineage %in% c("soft_tissue") & grepl("alveolar", subtype) |
           lineage %in% c("soft_tissue") & grepl("embryonal", subtype)) %>%
  mutate(subtype = case_when(lineage == "skin" & subtype == "undifferentiated" ~ "melanoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "skin" & subtype == "melanocytic" ~ "melanoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "skin" & subtype == "transitory" ~ "melanoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "skin" & subtype == "neural crest-like" ~ "melanoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "colorectal" & subtype == "adenocarcinoma" ~ "colon adenocarcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "pancreas" & subtype == "ductal adenocarcinoma, exocrine" ~ "pancreas adenocarcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "pancreas" & subtype == "adenocarcinoma" ~ "pancreas adenocarcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "prostate" & subtype == "adenocarcinoma" ~ "prostate adenocarcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "central_nervous_system" & subtype == "glioblastoma multiforme" ~ "glioma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "central_nervous_system" & subtype == "glioblastoma" ~ "glioma", TRUE ~ subtype))  %>%
  mutate(subtype = case_when(lineage == "breast" ~ "breast", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "esophagus" & subtype == "squamous cell carcinoma" ~ "esophageal carcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "lung" & subtype == "lung adenocarcinoma" ~ "non-small cell lung cancer", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "lung" & subtype == "lung squamous cell carcinoma" ~ "non-small cell lung cancer", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "lung" & subtype %like% "non-small cell lung cancer" ~ "non-small cell lung cancer", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "ovary" & subtype %like% "cystadenocarcinoma" ~ "ovarian serous cystadenocarcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "kidney" & subtype %like% "renal cell carcinoma" ~ "renal cell carcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "blood" & subtype %like% "acute lymphoblastic leukemia" ~ "acute lymphoblastic leukemia", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "thyroid" & subtype %like% "carcinoma" ~ "thyroid carcinoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "soft_tissue" & subtype %like% "synovial" ~ "synovial sarcoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "soft_tissue" & subtype %like% "alveolar" ~ "alveolar rhabdomyosarcoma", TRUE ~ subtype)) %>%
  mutate(subtype = case_when(lineage == "soft_tissue" & subtype %like% "embryonal" ~ "embryonal rhabdomyosarcoma", TRUE ~ subtype)) %>%
  filter(!(subtype %in% c("adenocarcinoma", "undifferentiated", "squamous cell carcinoma"))) 

#output directory for plots
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021")


#ensemble IDs and gene IDs
hgnc <- read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/hgnc_complete_set_7.24.2018.txt")
hgnc <- hgnc[,c(2,20)]

#importing metabolic genes
pathway_file <- "/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/CCLE metabolic signature analysis/Single-Cell-Metabolic-Landscape-master/Data/KEGG_metabolism.gmt"
pathways <-  gmtPathways(pathway_file)
KEGG_metabolic_genes <- unique(unlist(pathways))

#extracting ENSGs of metabolic genes
hgnc_metabolic_genes <- hgnc[hgnc$symbol %in% KEGG_metabolic_genes,]

#replacing ENSG ids with Gene names
Celligner_metabolic_data <- Celligner_data[,colnames(Celligner_data) %in% hgnc_metabolic_genes$ensembl_gene_id]

hgnc_metabolic_genes <- hgnc_metabolic_genes[hgnc_metabolic_genes$ensembl_gene_id %in% colnames(Celligner_metabolic_data),]
hgnc_metabolic_genes <- hgnc_metabolic_genes[order(hgnc_metabolic_genes$ensembl_gene_id),]

Celligner_metabolic_data <- Celligner_metabolic_data[,order(colnames(Celligner_metabolic_data))]
colnames(Celligner_metabolic_data) <- hgnc_metabolic_genes$symbol

#filtering to only cell lines/cancer types of interest
Celligner_metabolic_data <- Celligner_metabolic_data[rownames(Celligner_metabolic_data) %in% unique(metadata_filtered$sampleID),]

#cancer type colors for UMAP
cancer_colors <- c(`colon adenocarcinoma`=  "#abd23f",
                   `acute myeloid leukemia` = "#e13978",
                   `uveal melanoma` = "#c091e3",                        
                   `melanoma` = "#e491c1",
                   `Ewing sarcoma` = "#45a132", 
                   `prostate adenocarcinoma` = "#3870c9",
                   `pancreas adenocarcinoma` = "#75dfbb",
                   `glioma` = "#e6c241",
                   `osteosarcoma` = "#d74829",
                   `acute lymphoblastic leukemia` = "#9f55bb",
                   `ovarian serous cystadenocarcinoma` = "#6c55e2",
                   `hepatocellular carcinoma` = "#6a6c2c",
                   `non-small cell lung cancer` = "#51d5e0",
                   `renal cell carcinoma` = "#d1d684",
                   `breast` = "#d8ab6a",
                   `esophageal carcinoma` = "#c44c90",
                   `alveolar rhabdomyosarcoma` = "#56e79d",
                   `embryonal rhabdomyosarcoma` = "#e13978",
                   `synovial sarcoma` = "#dfbc3a",
                   `thyroid carcinoma` = "#d04850") 




#Set up R UMAP, optimizing parameters
for (nn in c(20)){
  for (md in c(0.5)){
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
      
      umap_output <- umap(as.matrix(Celligner_metabolic_data),config = custom.config)
      
      
      #Collect Output
      output = list(umap = umap_output$layout,
                    data.zscore = Celligner_metabolic_data,
                    parameters = data.frame(nn=nn,md=md,ncomp=ncomp,metric=metric))
      
      
      UMAP_df = output[["data.zscore"]]
      UMAP_df$UMAP_1=output[["umap"]][,1]
      UMAP_df$UMAP_2=output[["umap"]][,2]
      UMAP_df$type = metadata_filtered$type
      UMAP_df$sampleID = metadata_filtered$sampleID
      UMAP_df$lineage = metadata_filtered$lineage
      UMAP_df$subtype = metadata_filtered$subtype
      
      
      
      plot_embedding_categorical <- function(category,.df){
        dpal_func = jcolors_contin("pal4",reverse = TRUE,bias = 0.5)
        dpal = rev(dpal_func(length(unique(.df[[category]]))+1))
        ggplot(data = .df,aes(color = type, size = type, x = `UMAP_1`,y = `UMAP_2`,fill = .data[[category]]))+
          geom_point(stroke = 0.2 ,shape = 21)+
          guides(fill= guide_legend(override.aes = list(size=3)))+
          xlab("UMAP 1") +
          ylab("UMAP 2") +
          theme_classic() +scale_fill_manual(values = cancer_colors) + 
          scale_color_manual(values=c(`CL`='black', `tumor`='white')) + 
          scale_size_manual(values=c(`CL`=3, `tumor`=2))
      }
      
      
      g = plot_embedding_categorical('subtype',UMAP_df)
      ggsave(sprintf('/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/cancer_types_tumors_cell_lines_UMAP-nn=%0.0f-md=%0.0e-%s.pdf',nn,md,metric),plot = g,device = 'pdf',width = 12, height = 7,dpi=300)
      
    }
  }
}


####
####

#functions to calculate % cell lines aligned to tumor types (heatmap)
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumors_samples,]), t(Celligner_aligned_data[cl_samples,]),
                      use='pairwise')
  
  
  return(tumor_CL_cor)
}


cell_line_tumor_class <- function(x, dist_mat, ann_mat, k=25, decreasing = T) {
  names(x) <- rownames(dist_mat)
  x <- sort(x, decreasing = decreasing)
  tumor_type <- dplyr::filter(ann_mat, sampleID %in% names(x[1:k]))$subtype %>%
    table() %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    head(1) %>%
    .[['.']] %>% as.character()
  
  return(tumor_type)
}

get_cell_line_tumor_class <- function(tumor_CL_cor, alignment) {
  cl_tumor_classes <- apply(tumor_CL_cor, 2, function(x) cell_line_tumor_class(x, tumor_CL_cor, alignment)) %>% 
    as.character()
  names(cl_tumor_classes) <- colnames(tumor_CL_cor)
  
  return(cl_tumor_classes)
}
  
  
cell_line_tumor_class_plot <- function(cl_tumor_classes, alignment, tumor_CL_cor, filename) {
  cl_tissue_type <- dplyr::filter(alignment, type=='CL')
  #cl_tissue_type[grep('rhabdomyosarcoma', cl_tissue_type$subtype),'tissue'] <- 'rhabdomyosarcoma'
  rownames(cl_tissue_type) <- cl_tissue_type$sampleID
  classification_freq <- table(cl_tumor_classes, cl_tissue_type[colnames(tumor_CL_cor),'subtype']) %>% as.data.frame()
  classification_freq <- reshape2::dcast(classification_freq, cl_tumor_classes ~ Var2, value.var = 'Freq') %>%
    tibble::column_to_rownames('cl_tumor_classes')
  print(setdiff(intersect(unique(dplyr::filter(alignment, type=='CL')$subtype),
                          unique(dplyr::filter(alignment, type=='tumor')$subtype)), 
                rownames(classification_freq)))
  
  esophagus_tumor <- rep(0, ncol(classification_freq))
  classification_freq <- rbind(classification_freq,`esophagus`= esophagus_tumor) 
  common_types <- intersect(rownames(classification_freq), colnames(classification_freq))
  
  prop_agree <- sum(diag(as.matrix(classification_freq[common_types, common_types])))/sum(as.matrix(classification_freq[common_types, common_types]))
  
  for(i in 1:ncol(classification_freq)) {
    classification_freq[,i] <- classification_freq[,i]/sum(classification_freq[,i])
  }
  
  
  agreement <- diag(as.matrix(classification_freq[common_types, common_types]))
  agreement_CL <- agreement
  names(agreement_CL) <- common_types
  agreement_tumor <- agreement
  names(agreement_tumor) <- common_types
  
  agreement_CL <- base::sort(agreement_CL, decreasing=T)
  agreement_tumor <- base::sort(agreement_tumor, decreasing=T)
  
  
  
  
  classification_freq <- classification_freq[names(agreement_tumor), names(agreement_CL)]
  rownames(classification_freq) <- gsub("_", " ", rownames(classification_freq))
  colnames(classification_freq) <- gsub("_", " ", colnames(classification_freq))
  
  p <- pheatmap::pheatmap(classification_freq, 
                     border_color = heatmap_params$square_border_color, 
                     na_col= heatmap_params$na_color, 
                     cluster_rows = F, 
                     cluster_cols = F, 
                     main="", 
                     fontsize = heatmap_params$title_font_size,
                     fontsize_col = heatmap_params$column_font_size,
                     fontsize_row = heatmap_params$row_font_size,
                     width = 3.5,
                     height = 3,
                     #fontface = heatmap_params$font_face,
                     angle_col=90, 
                     color= heatmap_params$color_vector, 
                     display_numbers = round(classification_freq, digits = 2),
                     number_color = "white",
                     fontsize_number = 12)
  
}


heatmap_params <- list(
  na_color = '#666666',
  title_font_size = 6,
  row_font_size = 6,
  column_font_size = 6,
  font_face = "plain",
  square_border_color = "white",
  color_palette = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab'),
  color_vector = c('#e3e3e3', rev(grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab')(100)))
)


#function to generate correlation histograms
cell_line_tumor_distance_distribution <- function(alignment, tumor_CL_cor) {
  alignment$compare_types <- alignment$subtype

  common_cancer_types <- intersect(dplyr::filter(alignment, type=='tumor')$compare_types, 
                                   dplyr::filter(alignment, type=='CL')$compare_types)
  
  tumor_names <- character()
  CL_names <- character()
  dist_list <- numeric()
  tissue_types <- character()
  for(cancer in common_cancer_types) {
    cur_tumors <- dplyr::filter(alignment, type=='tumor' & compare_types==cancer)$sampleID
    cur_CLs <- dplyr::filter(alignment, type=='CL' & compare_types==cancer)$sampleID
    cur_dist <- reshape2::melt(as.matrix(tumor_CL_cor[cur_tumors, cur_CLs]))
    tumor_names <- c(tumor_names, as.character(cur_dist$Var1))
    CL_names <- c(CL_names, as.character(cur_dist$Var2))
    dist_list <- c(dist_list, cur_dist$value)
    tissue_types <- c(tissue_types, rep(cancer, nrow(cur_dist)))
    
  }
  
  dist_df <- cbind.data.frame(tumor_names, CL_names, dist_list, tissue_types)
  dist_df$tissue_types <- gsub("_", " ", dist_df$tissue_types)
  mean_dist <- aggregate(dist_df$dist_list, list(dist_df$tissue_types), 
                         FUN = quantile, probs = 0.25) %>% dplyr::arrange(desc(x))
  
  mean_dist$Group.1 <- rev(mean_dist$Group.1)
  dist_df$tissue_types <- factor(dist_df$tissue_types, levels = mean_dist$Group.1)
  
  tumor_dist_spread <- ggplot2::ggplot(dplyr::filter(dist_df, tissue_types != 'all'),
                                       ggplot2::aes(x = dist_list, y = tissue_types, fill = tissue_types)) +
    ggridges::geom_density_ridges(alpha=0.8) +
    ggridges::theme_ridges() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", text=ggplot2::element_text(size=6),
                   axis.text = ggplot2::element_text(size=6)) +
    ggplot2::xlab("correlation between cell lines and tumors") +
    ggplot2::ylab('cancer type') + 
    ggplot2::scale_fill_manual(values = cancer_colors)
  
  return(tumor_dist_spread)
  
}


####
#####

#executing functions above to generate heatmap and ridge plots

tumor_CL_Corr <- calc_tumor_CL_cor(Celligner_metabolic_data,metadata_filtered)

cell_line_tumor_class <- get_cell_line_tumor_class(tumor_CL_Corr, UMAP_df)

p <-cell_line_tumor_class_plot(cell_line_tumor_class,UMAP_df,tumor_CL_Corr, filename)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/heatmap_percent_aligned_correct.pdf',
       plot = p,device = 'pdf',width = 8, height = 8,dpi=300)


p <- cell_line_tumor_distance_distribution(UMAP_df,tumor_CL_Corr)

ggsave('/Volumes/FallahiLab/Maize-Data/Data/Cara/NatureComm_Warren_2021/corr_cell_lines_tumors_ridge_plots.pdf',
       plot = p,device = 'pdf',width = 5, height =8,dpi=300)