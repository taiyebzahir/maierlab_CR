############### General requirements #####################################################################

# A xlsx file named Tab1.xlsx with Exp_Plate, Replicate, Community (see template)
# A xlsx file named Tab2.xlsx with Plate_Design, Position,Condition, Drug_Conc, Community and Border (see template)
# A folder called "Raw Data" with Tecan Files. All Raw files should only contain Sheets with data, no empty sheets. All Sheet names need to end with the number of the Plate_Design
# R script called readTecanMeasurements.R


############### Workspace preparation #####################################################################

# set your working directory with Tab1.xlsx, Tab2.xlsx, folder with Raw Data and species annotation

setwd("~/Downloads/All screening data/screen post drug_500ul/OD")


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
  
############### Reading and combining your data  #####################################################################

Raw_Reads <- NULL

#loop through all files in the Raw Data folder and combine them

for (i in 1:length(Raw_Data_files)){
  Raw_Reads <- rbind(Raw_Reads, readTecanMeasurements(paste0('./Raw Data/',Raw_Data_files[i])))
}

Raw_Reads <- Raw_Reads %>% mutate(Exp_Plate = paste(Exp, Plate_Design, sep='_')) %>% 
  left_join(Runs, by='Exp') %>% left_join(Layout, by=c('identifier', 'Plate_Design'))

############### Data Analysis  #####################################################################

# substract baseline OD from all OD values (per plate)

Raw_Reads <- Raw_Reads %>% mutate(OD_corrected = Absorbance - 0.093)

# remove all Border values

Raw_Reads <- Raw_Reads %>% filter(!Border)

# normalise by "COM_control" control per column

Raw_Reads <- Raw_Reads %>% mutate(col = identifier %% 100)

Raw_Reads <- Raw_Reads %>% dplyr::group_by(Exp_Plate, col) %>% 
  mutate(OD_norm = (OD_corrected/OD_corrected[Condition=='COM_control'])) %>% ungroup()

write.csv(Raw_Reads, file = "./Results/Raw_Reads.csv")

####heatmap of OD#####

# making a dataframe with mean OD_normalised
overnight_meanOD <- NULL
overnight_meanOD <- Raw_Reads %>% 
  filter(Community == TRUE & Condition != 'DMSO' & Condition!= 'COM_control') %>% 
  group_by(identifier, Plate_Design) %>% mutate (mean_ODnorm = mean(OD_norm)) %>% ungroup() %>% 
  filter (Replicate == 1) %>% subset(select = c(Condition,Conc,mean_ODnorm)) 

write.csv(overnight_meanOD, "ONOD_normalised.csv")

#plotting heat map

OD_spread <- spread(overnight_meanOD, key = Conc, value = mean_ODnorm)
OD_spread <- OD_spread %>% remove_rownames %>% column_to_rownames(var="Condition")
OD_spread_40 <- NULL
OD_spread_160 <- NULL

for (row in 1:nrow(OD_spread)) {
  
    if(is.na(OD_spread[row,2.5])) {
    OD_spread_160 <- rbind(OD_spread_160, as.vector(OD_spread[row,]))
    }
  else {
    OD_spread_40 <- rbind(OD_spread_40, as.vector(OD_spread[row,]))
  }
}
OD_spread_160 <- OD_spread_160[c(-1,-2)]
OD_spread_40 <- OD_spread_40[c(-6,-7)]

pdf("ONOD_160.pdf", width=12, height=10)
heatmap.2(as.matrix(OD_spread_160),
          Colv = "NA", # this to make sure the columns are not ordered
          Rowv = TRUE, scale = "none",
          main = "Community OD",trace = "none",
          srtCol = 20,
          density.info="none",
          dendrogram = "none",
          colsep=1:ncol(OD_spread_160),
          rowsep=1:nrow(OD_spread_160),
          sepcolor="black",
          sepwidth=c(0.05,0.05),
          margins = c(5, 9))
dev.off()

pdf("ONOD_40.pdf", width=11, height=4)
heatmap.2(as.matrix(OD_spread_40),
          Colv = "NA", # this to make sure the columns are not ordered
          Rowv = TRUE, scale = "none",
          main = "Community OD",trace = "none",
          srtCol = 20,
          density.info="none",
          dendrogram = "none",
          colsep=1:ncol(OD_spread_40),
          rowsep=1:nrow(OD_spread_40),
          sepcolor="black",
          sepwidth=c(0.05,0.05),
          margins = c(5, 9))

dev.off()