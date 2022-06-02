###### code to analyse one-to-one interaction data from Chiara ######

##### setup #####
setwd("~/Downloads/Chiara")
pacman::p_load(ggplot2, openxlsx, tidyverse, gplots)
data <- read.xlsx("Data_long_Chiara.xlsx")
data <- data %>% filter(!(Commensal == "COM" | Commensal == "ED1a"))

##### data manipulation ####
# calculating average of log2 inhibition 
data <- data %>% mutate(log2CR = rowMeans(log2(data[,2:4]), na.rm=TRUE))
data <- data[-c(2,3,4)]
data_wide <- data %>% spread(Commensal, log2CR)
data_wide <- data_wide %>% column_to_rownames(var = "Pathogen")

##### heatmap #####

library(pheatmap)

pdf("onetoone.pdf", width=12, height=10)
pheatmap(t(as.matrix(data_wide)),
          scale = "column",
          cluster_rows = TRUE, cluster_cols = TRUE,
          color = colorRampPalette(c("blue", "white", "red")) (50),
          angle_col = 45,
          border_color="black"
)
dev.off()

write_csv(data, "onetoone_data.csv")
