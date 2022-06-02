# script to merge metadata like conditions to miseq percent species and plot barplots of relative specie abundance #

##### load libraries import data ###
pacman::p_load(ggplot2, dplyr, openxlsx, tidyverse)

setwd("~/Downloads/lab_seq_analysis/sequencing")

#import sequencing samples metadata
samples <- read.xlsx('210818-Microbiology_Metadata_Sheet_1.5_Taiyeb.xlsx', 3)
#samples <- read.xlsx('Mouse_samples.xlsx')
colnames(samples) <- samples[1,]
samples <- samples[-1,]
samples <- samples %>% select(1,8,11,12)
samples <- samples %>% mutate(drug_conc = paste(`Condition: tag`, `Condition: tag2`, sep=' ')) 
samples <- samples [1:192,] # sample grouping is the replicate 

#import sequencing reads and tabulate in 1 table
source('sequencing_reads_import.R')
Reads_files <- list.files('./Screen_tables/Tables')  # species also below 1 %
#Reads_files <- list.files('./Mouse_tables/Tables')  # species also below 1 %
Seq_reads <- NULL
for (i in 1:length(Reads_files)){
  Seq_reads <- bind_rows(Seq_reads, sequencing_reads_import(paste0('./Screen_tables/Tables/',Reads_files[i])))
}
Seq_reads[is.na(Seq_reads)] <- 0
Seq_reads <- rownames_to_column(Seq_reads, "Sample Number")

##### merge the metadata with seq_reads #####
seq_data <- merge(Seq_reads,samples, by = "Sample Number")
names(seq_data)[names(seq_data)=="Sample Number"] <- "Sample_Number"
write.csv(seq_data, file = "seq_data_screen.csv", row.names = FALSE)

##### bar plots all drugs #####

colors_bar <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "#CAB2D6", # lt purple
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "brown"
)

seq_data <- seq_data %>% gather(Species, rel.abundance, NT5001:NT5076)
grouped_data <- seq_data %>% group_by(drug_conc)
grouped_data <- group_split(grouped_data)
for (p in grouped_data) {
  condition <- p[1,5]
  myplot <- ggplot(p, aes(fill=Species, y=rel.abundance, x=Sample)) + 
    geom_bar(position="stack", stat="identity") + ggtitle(condition) + theme_bw() + 
    scale_fill_manual(values = colors_bar)
    ggsave(filename=paste0('./abundance_plotsnew/',condition,".pdf"), 
           plot = myplot, device = cairo_pdf) 
}

##### barplots_mouse ####

day2 <- Seq_reads[c(1,5,9,13,17),]
day2 <- rownames_to_column(day2, "Sample")
p <- day2 %>% gather(Species, rel.abundance, NT5001:NT5025)
myplot <- ggplot(p, aes(fill=Species, y=rel.abundance, x=Sample)) + 
  geom_bar(position="stack", stat="identity") + theme_bw() + 
  scale_fill_manual(values = colors_bar)
ggsave(filename="plot2.pdf", plot = myplot, device = cairo_pdf)

day6 <- Seq_reads[c(2,6,10,14,18),]
day6 <- rownames_to_column(day6, "Sample")
p <- day6 %>% gather(Species, rel.abundance, NT5001:NT5025)
myplot <- ggplot(p, aes(fill=Species, y=rel.abundance, x=Sample)) + 
 geom_bar(position="stack", stat="identity") + theme_bw() +
  scale_fill_manual(values = colors_bar)
ggsave(filename="plot6.pdf", plot = myplot, device = cairo_pdf)

day57 <- Seq_reads[c(3,7,11,15,19),]
day57 <- rownames_to_column(day57, "Sample")
p <- day57 %>% gather(Species, rel.abundance, NT5001:NT5025)
myplot <- ggplot(p, aes(fill=Species, y=rel.abundance, x=Sample)) + 
 geom_bar(position="stack", stat="identity") + theme_bw() +
  scale_fill_manual(values = colors_bar)
ggsave(filename="plot57.pdf", plot = myplot, device = cairo_pdf)



#####PCoA plot#####

library(ape)
library(vegan)
distance <- vegdist(data, "bray")
res <- pcoa(distance)
biplot(res) 
