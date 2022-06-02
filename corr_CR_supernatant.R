#### script compares supernatant data to onetoone interaction CR experiments 
### by Chiara

setwd("~/Downloads/Chiara")
pacman::p_load(ggplot2, openxlsx, tidyverse, gplots, plotly, ggpubr)
supCR_data <- read.xlsx("data_long_supernatant.xlsx", sheet = 2)
CR_data <- read.xlsx("Data_long_Chiara.xlsx")
CR_data <- CR_data %>% mutate(mean_CR = (Replicate1 + Replicate2 + Replicate3)/3 ) %>%
  select(-c(Replicate1,Replicate2,Replicate3))
supCR_data <- supCR_data %>% mutate(mean_supCR = (Replicate1 + Replicate2 + Replicate3)/3 ) %>%
  select(-c(Replicate1,Replicate2,Replicate3))
data_total <- merge(CR_data, supCR_data, by = c("Commensal", "Pathogen")) %>% 
  filter(!(Commensal == "ED1a" | Commensal == "Vp" | Commensal == "COM" | Commensal == "Er"))

# filtering ER Eubacterium rectale because of suspicion that it was not Er in the supernatant experiments


##### plotly #####

library(ggpmisc)

scatterPlot <- data_total %>% 
  ggplot(aes(x = mean_CR, y = mean_supCR)) + 
  geom_point(alpha=0.5, size= 3,  aes(text = Commensal)) + 
  labs(x = "Rel. pathogen growth in competition assay", 
       y = "Rel. pathogen growth in supernatant") + facet_wrap(~Pathogen) +
  stat_cor (method = "spearman", cor.coef.name = "rho") + 
  geom_smooth(method = lm, color = "blue") +
  theme_classic()

p <- ggplotly(scatterPlot, tooltip = "text")
htmlwidgets::saveWidget(as_widget(p), "CRvssupCR.html")

pdf("CRvssupCR.pdf", width=8, height=6)

data_total %>% 
  ggplot(aes(x = mean_CR, y = mean_supCR)) + 
  geom_point(alpha=0.5, size= 3) + 
  labs(x = "Rel. pathogen growth in competition assay", 
       y = "Rel. pathogen growth in supernatant") + facet_wrap(~Pathogen) +
  stat_cor (method = "spearman", cor.coef.name = "rho") + 
  geom_smooth(method = lm, color = "blue") +
  theme_classic()

dev.off()

