rm(list=ls()) #clear all
cat("\014") #clc

#importing libraries needed
library(data.table)
library(ggplot2)



#output directory
outDir <- file.path("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Cellosaurus")

#input OXPHOS high/low cell lines from DepMap
cell_lines_DepMap <- read.csv("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Cellosaurus/OXPHOS_variable_cell_lines_mutation_KO_data.csv")

#input of cellosaurus data
cellosaurus_txt <- as.data.frame(read.delim("/Volumes/FallahiLab/Maize-Data/Data/Cara/CCLE data/Cellosaurus/cellosaurus.txt"))
colnames(cellosaurus_txt) <- "txt"


#keeping only entries that has cell line IDs and doubling information
cellosaurus_txt_v2 <- as.data.frame(cellosaurus_txt[cellosaurus_txt$txt %like% "ID   " | 
              cellosaurus_txt$txt %like% "CC   Doubling",])
colnames(cellosaurus_txt_v2) <- "text_info"


#identifying cell index of cell names
cell_line_idx <- which(cellosaurus_txt_v2$text_info %like% "ID   ")


#finding doubling times for each cell line

# Initialize an empty list to store the data frames
result_list <- list()

# Iterate over the indices of unique iterations
result_list <- parallel::mclapply(1:length(cell_line_idx), function(i) {
  
  c <- cell_line_idx[i]
  
  df <- as.data.frame(cellosaurus_txt_v2[c(c:(1+c)),])
  
  doubling_time_exist <- sum(df$`cellosaurus_txt_v2[c(c:(1 + c)), ]` %like% "CC   Doubling")
  
  if (doubling_time_exist == 1) {
    double_time <- as.character(df[df$`cellosaurus_txt_v2[c(c:(1 + c)), ]` %like% "CC   Doubling",])
  } else double_time <- NA
  
  cell_line_name <- cellosaurus_txt_v2[c,]
  
  cell_line_double_time <- data.frame(cell_line = cell_line_name, Doubling_Time = double_time)
 
  # Return the data frame for each iteration
  return(cell_line_double_time)
}, mc.cores = 10) #defining number of cores to use for parallel processing

# Combine the list of data frames into a single data frame
cell_line_double_time.df <- do.call(rbind, result_list)
cell_line_double_time.df <- na.omit(cell_line_double_time.df)

#modifying cell line and double time in data frame 
cell_line_double_time.df$cell_line <- gsub("^.{1,5}", "", cell_line_double_time.df$cell_line)
cell_line_double_time.df$cell_line <- gsub("-", "", cell_line_double_time.df$cell_line)
cell_line_double_time.df$cell_line <- gsub(" ", "", cell_line_double_time.df$cell_line)
cell_line_double_time.df$Doubling_Time  <- gsub(".*:\\s*", "", cell_line_double_time.df$Doubling_Time)
cell_line_double_time.df$Doubling_Time <- gsub("\\s*\\(.*", "", cell_line_double_time.df$Doubling_Time)

#keep cell lines with doubling time having hours/days
cell_line_double_time.df <- cell_line_double_time.df[cell_line_double_time.df$Doubling_Time %like% "Days" |
                                                       cell_line_double_time.df$Doubling_Time %like% "hour",]

#merging DepMap and cellosaurus cell lines
DepMap_cellosaurus.df <- merge(cell_lines_DepMap[,c(1:4)], cell_line_double_time.df, 
                               by.x = "StrippedCellLineName", by.y = "cell_line")

#cleaning up doubling times to get  numbers
DepMap_cellosaurus.df$Doubling_Time <- gsub("~", "", DepMap_cellosaurus.df$Doubling_Time)
DepMap_cellosaurus.df$Doubling_Time <- gsub("hours", "", DepMap_cellosaurus.df$Doubling_Time)
DepMap_cellosaurus.df$Doubling_Time <-gsub("\\+\\-.*", "", DepMap_cellosaurus.df$Doubling_Time)

#converting doubling times that has a range by calculating the mean
cell_lines_with_doubling_range_idx <- which(DepMap_cellosaurus.df$Doubling_Time %like% "-")

cell_lines_with_doubling_range <- DepMap_cellosaurus.df[cell_lines_with_doubling_range_idx,]

for (i in 1:nrow(cell_lines_with_doubling_range)) {
  numbers <- as.numeric(strsplit(cell_lines_with_doubling_range$Doubling_Time, "-")[[i]])

  average_value <- mean(numbers)
  
  cell_lines_with_doubling_range$Doubling_Time[i] <- average_value
}

#final data frame for analysis
DepMap_cellosaurus.df <- rbind(cell_lines_with_doubling_range,DepMap_cellosaurus.df[-cell_lines_with_doubling_range_idx,])
DepMap_cellosaurus.df$Doubling_Time <- as.numeric(DepMap_cellosaurus.df$Doubling_Time)


#t-test of doubling time of cell lines and boxplots
df <- DepMap_cellosaurus.df[!(DepMap_cellosaurus.df$OXPHOS_state == 0),]; df$OXPHOS_state <- factor(df$OXPHOS_state, levels = c("OXPHOS_low", "OXPHOS_high"))
stat.test <- t.test(df$Doubling_Time~df$OXPHOS_state)

#boxplot of doubling times in OXPHOS high vs. OXPHOS cell line
p <- ggplot(df,aes(x = OXPHOS_state, y = Doubling_Time, fill = OXPHOS_state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + labs(title = "Cell line doubling times (Cellosaurus)") + xlab("") + ylab("Doubling time (hr)") +
  theme_classic() + theme(legend.position="none",panel.grid.major = element_blank(),
                          panel.grid.minor= element_blank(),axis.text=element_text(colour = "black", size=12),
                          axis.line=element_line(linewidth=0.5,colour="black"),
                          axis.ticks = element_line(colour = "black",linewidth=0.5),
                          axis.text.x = element_text(colour="black", size = 12, angle = 30, hjust=1),
                          axis.text.y=element_text(colour="black", size = 12),
                          strip.background = element_rect(fill="white",linewidth=NA,colour = NULL),
                          strip.text=element_text(size=8), text = element_text(size = 12)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.7), size = 1.3, alpha = 0.6) + scale_fill_manual(values=c("skyblue1", "tomato")) + 
  annotate("text", x = 1.4, y = 100, label = paste("p =", format(stat.test[["p.value"]], scientific = TRUE, digits = 3))) + 
  annotate("text", x = 1, y = 120, label = paste("OXPHOS_low =", as.numeric(table(df$OXPHOS_state)[1]))) + 
  annotate("text", x = 2, y = 120, label = paste("OXPHOS_high =", as.numeric(table(df$OXPHOS_state)[2]))) + ylim(0,150)

ggsave(file.path(outDir, 'cell_line_doubling_times_across_OXPHOS_states.pdf'), plot = p,device = 'pdf',width = 5, height = 8,dpi=300)

