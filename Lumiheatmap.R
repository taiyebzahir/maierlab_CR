############### General requirements ##################################################

# A xlsx file named Tab1_v1.xlsx with Exp_Plate, Replicate, Community (see template)
# A xlsx file named Tab2_v2.xlsx with Plate_Design, Position,Condition, Drug_Conc, 
# Community and Border (see template)
# A folder called "Raw Data" with Tecan Files. All Raw files should only contain Sheets 
# with data, no empty sheets. All Sheet names need to end with the number of the Plate_Design
# R script called readTecanMeasurements.R


############### Workspace preparation #################################################

# set your working directory with Tab1.xlsx, Tab2.xlsx, folder with 
# Raw Data and species annotation

setwd("~/Downloads/All screening data/screen post drug_500ul/Lumi")


# load libraries
pacman::p_load(ggplot2, dplyr, openxlsx, ggpubr, tidyr, gridExtra, purrr, tidyverse, gplots) 

# source growth_curve_functions
source('./readTecanMeasurements.R')

# create folders for Results
dir.create("Results")

# load files, scripts and annotations
Raw_Data_files <- list.files('Raw Data')
Runs <- read.xlsx('Tab1_v1.xlsx', 1)
Runs <- Runs %>% modify_if(is.factor, as.character)  #line added later

Layout <- NULL
  for (i in 1:6) {
    Layout <- rbind(Layout, read.xlsx('Tab2_v2.xlsx', i))
  }
  
############### Reading and combining your data  ###########################################

Raw_Reads <- NULL

#loop through all files in the Raw Data folder and combine them

for (i in 1:length(Raw_Data_files)){
  Raw_Reads <- rbind(Raw_Reads, readTecanMeasurements(paste0('./Raw Data/',Raw_Data_files[i])))
}

Raw_Reads <- Raw_Reads %>% mutate(Exp_Plate = paste(Exp, Plate_Design, sep='_')) %>% 
  left_join(Runs, by='Exp') %>% left_join(Layout, by=c('identifier', 'Plate_Design'))

############### Data Analysis  #####################################################

# remove all Border values

Raw_Reads <- Raw_Reads %>% filter(!Border)

# normalise by "COM_control" control per column

Raw_Reads <- Raw_Reads %>% mutate(col = identifier %% 100)

Raw_Reads <- Raw_Reads %>% dplyr::group_by(Exp_Plate, col) %>% 
  mutate(signal_norm = (Luminescence/Luminescence[Condition=='COM_control'])) %>% ungroup()

write.csv(Raw_Reads, file = "./Results/Raw_Reads.csv")

####heatmap of luminescence#####

Raw_Reads <- read.csv("./Results/Raw_Reads.csv")

# making a dataframe with mean log2 fold change in luminescence signal
fold_change_growth <- NULL
fold_change_growth <- Raw_Reads %>% filter(Community == TRUE & Condition != 'DMSO' & 
      Condition!= 'COM_control') %>% group_by(identifier, Plate_Design) %>%
  mutate (log2FC_lumi = mean(log2(signal_norm))) %>% filter (Replicate == 1) %>% 
  subset(select = c(Condition,Conc,log2FC_lumi)) %>% 
  mutate(drug_conc = paste0(Condition, " ", Conc))
write.csv(fold_change_growth, "log2FC_lumiscreen.csv")

#plotting heat map of log_fold_change in S.Tm lumi normalized signal

growth_spread <- fold_change_growth[,-4] %>% spread(key = Conc, value = log2FC_lumi)
growth_spread <- growth_spread %>% remove_rownames %>% column_to_rownames(var="Condition")
growth_spread_40 <- NULL
growth_spread_160 <- NULL

for (row in 1:nrow(growth_spread)) {
  
    if(is.na(growth_spread[row,2.5])) {
    growth_spread_160 <- rbind(growth_spread_160, as.vector(growth_spread[row,]))
    }
  else {
    growth_spread_40 <- rbind(growth_spread_40, as.vector(growth_spread[row,]))
  }
}
growth_spread_160 <- growth_spread_160[c(-1,-2)] 
growth_spread_40 <- growth_spread_40[c(-6,-7)]

breaks_plot = seq(-5,5,length.out=31)
gradient1 = colorpanel( sum( breaks_plot[-1]<=0 ), "navy", "white" ) 
gradient2 = colorpanel( sum( breaks_plot[-1]>0 ), "white", "red" )
hm.colors = c(gradient1,gradient2)


pdf("CR_160.pdf", width=12, height=10)
heatmap.2(as.matrix(growth_spread_160),
          Colv = "NA", # this to make sure the columns are not ordered
          main = "Log2 fold change S.Tm growth",trace = "none",
          srtCol = 20, 
          breaks=breaks_plot,col=hm.colors,
          density.info="none",
          dendrogram = "none",
          colsep=1:ncol(growth_spread_160),
          rowsep=1:nrow(growth_spread_160),
          sepcolor="black",
          sepwidth=c(0.05,0.02),
          margins = c(5, 9)
          )
dev.off()

pdf("CR_40.pdf", width=11, height=4)
heatmap.2(as.matrix(growth_spread_40),
          Colv = "NA", # this to make sure the columns are not ordered
          trace = "none", scale = "none",
          srtCol = 20,
          breaks=breaks_plot,col=hm.colors,
          density.info="none",
          dendrogram = "none",
          colsep=1:ncol(growth_spread_40),
          rowsep=1:nrow(growth_spread_40),
          sepcolor="black",
          sepwidth=c(0.05,0.05),
          margins = c(5, 9)
)
dev.off()

