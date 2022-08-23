library("ggplot2")
library("tidyverse")
library(lubridate)
library("ggraph")
library(igraph)
library(RColorBrewer)
library(beeswarm)
library(cowplot)

setwd("~/Projects/Staph_military/Reanalysis/Transmission/")

dists <- read.csv("Staph_2015_distances.csv")
dists$Site1 <- substr(dists$ID1,7,7)
dists$Site2 <- substr(dists$ID2,7,7)
dists$Same_site <- dists$Site1 == dists$Site2
dists$Same_patient <- dists$Patient1 == dists$Patient2
dists$Same_Company <- dists$Company1 == dists$Company2

metadata <- read.csv("../Metadata/MLST_profiles_metadata.csv")
dates <- read.csv("../Metadata/Correct_dates_2015.csv")

within_host <- dists %>% filter(Same_patient == TRUE) 
between_host <- dists %>% filter(Same_patient == FALSE)

summary(within_host$Dist)
quantile(within_host$Dist, 0.95)

dists %>% ggplot(aes(x = Dist, fill = Same_patient, color = Same_patient,  bw = 0.5)) + geom_density(alpha = 0.5) + facet_grid(.~ST) + xlim(0,30) 

dists_25snps <- filter(dists, Dist < 26, Same_patient == FALSE)
network25 <- graph_from_data_frame(d=dists_25snps, directed=F) 
plot(network25, vertex.size = 5)
components_25 <- components(network25, mode = c("weak", "strong"))

# Within_plot

dists %>%  ggplot(aes(x = Same_patient ~ ST, y = Dist, fill = ST)) +
  geom_boxplot() + scale_x_discrete(labels=c("FALSE" = "Between hosts", "TRUE" = "Within host")) + scale_y_continuous(trans='log10') +
  xlab("") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2", "#d8b365")) + 
  ggtitle("Pairwise SNP distances by Sequence Type (ST)") + theme(plot.title = element_text(hjust = 0.5))

#scale_fill_brewer( type = "qual")


wilcox.test(Dist ~ Same_site, data=within_host)

Company_plot <- dists %>% ggplot(aes(x = Dist, fill = Same_Company, group = Same_Company)) +
  geom_density(alpha = 0.5) +  xlab("Pairwise SNP distance") + ylab("Relative Frequency") + 
  scale_fill_manual(values=c("#B9BFFF", "#F4B6B2"),labels=c("Different","Same")) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Company"))

Patient_plot <- dists %>% ggplot(aes(x = Dist, fill = Same_patient, group = Same_patient)) +
  geom_density(alpha = 0.5) +  xlab("Pairwise SNP distance") + ylab("Relative Frequency") + 
  scale_fill_manual(values=c("#B9BFFF", "#F4B6B2"),labels=c("Different","Same")) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Patient"))

SNPS_plot <- dists %>% ggplot(aes(x = Same_patient, y = Dist, fill = as.factor(ST))) +
  geom_boxplot() + scale_x_discrete(labels=c("FALSE" = "Between hosts", "TRUE" = "Within host")) + scale_y_continuous(trans='log10') +
  xlab("") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2", "#d8b365")) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="ST"))


Density_grid <- plot_grid(Company_plot,Patient_plot, nrow = 2, labels = c("B","C"))
plot_grid(SNPS_plot,Density_grid, labels = c("A",""))


