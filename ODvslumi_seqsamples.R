### code to generate ODvslumi plot and show sequencing samples on the graph

setwd("~/Downloads/All screening data/screen post drug_500ul/OD")

# load libraries
pacman::p_load(ggplot2, dplyr, openxlsx, ggpubr, tidyr, gridExtra, purrr, tidyverse, gplots) 

##### ODVSlumi plot ####
# Is community-OD and SB-Lumi anticorrelated or are there exceptions? 

OD <- read.csv("ONOD_normalised.csv")
lumi <- read.csv("log2FC_lumiscreen.csv")
Data <- OD %>% left_join(lumi, by=c('Condition', 'Conc'))
data_com20 <- read.xlsx("COM20_diluted.xlsx")
data_com20 <- data_com20 %>% mutate(log2FC_lumi = log2(Lumi_norm))
library(ggpmisc)
library(plotly)
color_points <- colorpanel(7, "blue", "#ffe5bf", "red")
pdf("ODvslumi.pdf", width=9, height=6)
ggplot() + 
  geom_point(data = Data, aes(x = mean_ODnorm, y = log2FC_lumi, col = factor(Conc), 
                              text = drug_conc),
             alpha=0.7, size= 3) +
  scale_color_manual(values = color_points) +
  geom_point(data = data_com20, aes(x = OD_norm, y = log2FC_lumi), size = 2) +
  geom_smooth(aes(x = OD_norm, y = log2FC_lumi), data = data_com20,  
              method=lm, se=FALSE, color = "black") +
  stat_poly_eq(aes(x = OD_norm, y = log2FC_lumi, 
                   label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               data = data_com20, formula = y ~ x,
               parse = TRUE) +
  labs(x = "Relative community OD", 
       y = "Log2 fold change S.Tm growth") +
  theme_classic()

dev.off()

p <- ggplotly(scatterPlot, tooltip = "text")
htmlwidgets::saveWidget(as_widget(p), "ODvslumi.html")


####### seq samples ##########
ggplot(Data, aes(x = mean_OD, y = mean_log_signal)) + 
  geom_point(alpha = 0.5, aes(size = Group, colour = Group)) + 
  scale_size_manual(values = c(2, 5)) +
  scale_color_manual(values = c("A" = "steelblue",
                                "B" = "darkorange")) + 
  geom_smooth(method=lm, se=FALSE) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  labs(x = "Community OD", 
       y = "Log2 S.Tm luminescence") +
  theme_classic()

