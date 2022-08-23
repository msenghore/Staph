library("ggplot2")
library("tidyverse")
library(lubridate)
library("ggraph")
library(igraph)
library(RColorBrewer)
library(beeswarm)

setwd("~/Projects/Staph_military/Reanalysis/Distance_analysis/")

dists <- read.csv("Gubbins_2015_distances_metadata.csv")
dists$Same_site <- dists$Site1 == dists$Site2

metadata <- read.csv("../Metadata/Staph_2015_mlst_lab_metadata.csv")
dates <- read.csv("../Metadata/Correct_dates_2015.csv")

within_host <- dists %>% filter(Same_patient == TRUE, ST1 == ST2) %>% filter(Dist < 1000)

between_host <- dists %>% filter(Same_patient == FALSE, ST1 == ST2) %>% filter(Dist < 1000)

wilcox.test(Dist ~ Same_site, data=within_host)
wilcox.test(Dist ~ ST1, data=(within_host %>% filter(ST1 %in% c(5,8))))
wilcox.test(Dist ~ ST1, data=(within_host %>% filter(ST1 %in% c(87,8))))
wilcox.test(Dist ~ ST1, data=(within_host %>% filter(ST1 %in% c(5,87))))

kruskal.test(Dist ~ ST1, data=(within_host %>% filter(ST1 %in% c(5,87,8))))
kruskal.test(Dist ~ ST1, data=(between_host %>% filter(ST1 %in% c(5,87,8))))

WH8 <- within_host %>% filter(ST1 == 8) %>% pull(Dist)
WH5 <- within_host %>% filter(ST1 == 5) %>% pull(Dist)
WH87 <- within_host %>% filter(ST1 == 87) %>% pull(Dist)

wilcox.test(Dist ~ Same_Company, data=between_host)

aov(Dist ~ ST1, data=(within_host %>% filter(ST1 %in% c(5,87,8)))) %>% summary()
aov(Dist ~ ST1, data=(between_host %>% filter(ST1 %in% c(5,87,8)))) %>% summary()

same_host_site <- dists %>% dplyr::filter(Same_patient == TRUE, Same_site == TRUE, Dist < 1000)

aggregate(between_host$Dist, list(between_host$Same_Company), FUN=mean) 

between_host %>% ggplot(aes(x = Same_Company, y = Dist, fill = Same_Company)) +
  geom_boxplot() +  scale_x_discrete(labels=c("FALSE" = "Different", "TRUE" = "Same")) + scale_y_continuous(trans='log10') +
  xlab("Company") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2")) + 
  theme(plot.title = element_text(hjust = 0.5))

bw_Host_Same_Company %>% ggplot(aes(x = Same_Company, y = Dist, fill = Same_platoon)) +
  geom_boxplot() +  scale_x_discrete(labels=c("FALSE" = "Different", "TRUE" = "Same")) + scale_y_continuous(trans='log10') +
  xlab("Platoon") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2")) + 
  theme(plot.title = element_text(hjust = 0.5))

Company_plot <- between_host %>% ggplot(aes(x = Dist, fill = Same_Company, group = Same_Company)) +
  geom_density(alpha = 0.7) +  xlab("") + ylab("Relative Frequency") + 
  scale_fill_manual(values=c("#B9BFFF", "#F4B6B2"), labels=c("Different","Same")) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Company"))

Platoon_plot <- bw_Host_Same_Company %>% ggplot(aes(x = Dist, fill = Same_platoon, group = Same_platoon)) +
  geom_density(alpha = 0.7) +  xlab("Pairwise SNP distance") + ylab("Relative Frequency") + 
  scale_fill_manual(values=c("#B9BFFF", "#F4B6B2"),labels=c("Different","Same")) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_legend(title="Platoon"))

Density_grid <- plot_grid(Company_plot,Platoon_plot, nrow = 2)

bw_Host_Same_Company <- between_host %>% filter(Same_Company == TRUE)

bw_Host_Same_Company$Same_platoon <- bw_Host_Same_Company$Platoon1 == bw_Host_Same_Company$Platoon2

aggregate(bw_Host_Same_Company$Dist, list(bw_Host_Same_Company$Same_platoon), FUN=summary)
wilcox.test(Dist ~ Same_platoon, data=(bw_Host_Same_Company))

aggregate(between_host$Dist, list(between_host$Same_Company), FUN=summary)
wilcox.test(Dist ~ Same_Company, data=(between_host))

summary(within_host$Dist)
quantile(within_host$Dist, 0.95)

### Plot box plot of within host diversity for ST 5,8 and 87

# Within_plot <- 
  
Dists_plot <-   dists %>% filter(ST1 %in% c("5","8","87"), ST1 == ST2) %>% ggplot(aes(x = Same_patient, y = Dist, fill = ST1)) +
    geom_boxplot() + scale_x_discrete(labels=c("FALSE" = "Between hosts", "TRUE" = "Within host")) + scale_y_continuous(trans='log10') +
    xlab("") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2", "#d8b365")) + 
    ggtitle("Pairwise SNP distances by Sequence Type (ST)") + theme(plot.title = element_text(hjust = 0.5))
  
  
  within_host %>% filter(ST1 %in% c("5","8","87"), ST1 == ST2) %>% ggplot(aes(x = Same_site, y = Dist, fill = ST1)) +
    geom_boxplot() + scale_x_discrete(labels=c("FALSE" = "Different Site", "TRUE" = "Same Site")) + scale_y_continuous(trans='log10') +
    xlab("") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2", "#d8b365")) + 
    ggtitle("Pairwise SNP distances by Sequence Type (ST)") + theme(plot.title = element_text(hjust = 0.5))
  #scale_fill_brewer( type = "qual")

  
  plot_grid(Dists_plot,Density_grid)
### Plot density plots for SNP distance

dists %>% filter (Dist < 26) %>% ggplot(aes(x = Dist, fill = Same_patient, color = Same_patient,  bw = 0.5)) + geom_density(alpha = 0.5) 
dists %>% filter (Dist < 101) %>% ggplot(aes(x = Dist, fill = Same_patient, color = Same_patient,  bw = 0.5)) + geom_density(alpha = 0.5) 


### Plot SNP_distances for comparisons from the same patiets

within_host %>% ggplot(aes(x = as.factor(ST1), y = Dist, color = ST1, fill = ST1)) + geom_boxplot() + ylab("SNP Distance") + xlab("Sequence type") 
within_host %>% ggplot(aes(x = as.factor(Company1), y = Dist, color = ST1, fill = ST1)) + geom_boxplot() + ylab("SNP Distance") + xlab("Company") 

within_host %>% filter(ST1 == ST2) %>% ggplot(aes(x = ST1, y = Dist, color = ST1, fill = ST1)) + geom_boxplot() + ylab("SNP Distance") + xlab("Sequence type") 


### Study within host diversity in ST8 and potentially, transmission networks
ST8_withinhost <- within_host %>% filter(ST1 == ST2, ST1 == 8)

##3 Run a logistic regression
SNP_rateST8 <- lm(Dist ~ as.numeric(Days) + Same_site, data = ST8_withinhost)
summary(SNP_rateST8) 
summary(ST8_withinhost$Dist)


### Plot a transmission network for everything under 9 snps

dists_9snps <- filter(dists, Dist < 11, Same_patient == FALSE)
network9 <- graph_from_data_frame(d=dists_9snps, directed=F) 
plot(network9, vertex.size = 5)
summary(network9)

dists_20snps <- filter(dists, Dist < 21, Same_patient == FALSE)
network20 <- graph_from_data_frame(d=dists_20snps, directed=F) 
plot(network20, vertex.size = 5)

components(network20, mode = c("weak", "strong"))

dists_20snps_ST5 <- filter(dists_20snps, ST1 == 5)
network20_ST5 <- graph_from_data_frame(d=dists_20snps_ST5, directed=F) 
plot(network20_ST5, vertex.size = 5)
components(network20_ST5, mode = c("weak", "strong"))


dists_15snps <- filter(dists, Dist < 16, Same_patient == FALSE)
network15 <- graph_from_data_frame(d=dists_15snps, directed=F) 
plot(network15, vertex.size = 5)
components(network15, mode = c("weak", "strong"))

dists_25snps <- filter(dists, Dist < 26, Same_patient == FALSE)
network25 <- graph_from_data_frame(d=dists_25snps, directed=F) 
plot(network25, vertex.size = 5)
components_25 <- components(network25, mode = c("weak", "strong"))



dists_25snps_ST5 <- filter(dists, Dist < 26, Same_patient == FALSE, ST1 == 5)
network25_ST5 <- graph_from_data_frame(d=dists_25snps_ST5, directed=F) 
plot(network25_ST5, vertex.size = 5)
components(network25_ST5, mode = c("weak", "strong"))

### Plot a network with metadata

## Use a threshold of 25 SNps to define clusters, summarize cluster

clusters <- components_25$membership %>% as.data.frame()
clusters$ID <- rownames(clusters)
clusters$cluster <- clusters$.
clusters <- dplyr::select(clusters, ID, cluster)
clusters$patient <- substr(clusters$ID,1,4)
clusters$visit <- substr(clusters$ID,6,6)
clusters$site <- substr(clusters$ID,7,7)

clusters <- merge(clusters,metadata,by="ID", all.x = TRUE)
clusters <- clusters %>% filter(ID != "CP000730.1")
clusters <- merge(clusters,dates,by="ID", all.x = TRUE)

### Generate tile plots for clusters
clusters %>% filter(cluster == 5) %>% ggplot(aes(x = visit, y = Site, fill = ST)) + 
  geom_tile() + facet_wrap( ~ patient, ncol = 5) + theme_bw()

clusters %>% filter(cluster == 1) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, ncol = 3) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 2) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 3) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 15) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 4) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 5) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 5) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 11) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 12) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 13) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 14) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 15) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 10) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")

clusters %>% filter(cluster == 16) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site)) +
  geom_tile(alpha = 0.5, width = 6) + facet_wrap(~patient, ncol = 4) + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")


cluster_summary <- dplyr::select(clusters, cluster) %>% unique()
rownames(cluster_summary) <- seq(1,17,1)

nares_25snps <- filter(dists, Dist < 26, Same_patient == FALSE, Site1 == "Nasal", Site2 == "Nasal")
network25_nares <- graph_from_data_frame(d=nares_25snps, directed=F) 
plot(network25_nares, vertex.size = 5)
components(network25_nares, mode = c("weak", "strong"))


### Distribution of body sites in ST8, population vs clusters

metadata$In_cluster <- metadata$ID %in% clusters$ID 
ST8meta <- metadata %>% filter(ST == 8)

table(ST8meta$Site, ST8meta$In_cluster )

ST5meta <- metadata %>% filter(ST == 5)
tab5 <- table(ST5meta$Site, ST5meta$In_cluster )
chisq.test(tab5)

### For each cluster get the number of isolates, number of patients, mean and median SNP Dist and ST

for(i in 1:17){
  qids <- ""
  qdists <- ""
  
  cluster_summary$N_isolates[i] <- clusters %>% filter(cluster == i) %>% nrow
  cluster_summary$N_patients[i] <- clusters %>% filter(cluster == i) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$ST[i] <- clusters %>% filter(cluster == i) %>% pull(ST) %>% unique()
  cluster_summary$N_wounded[i] <- clusters %>% filter(cluster == i, site.x == "C") %>% pull(patient) %>% unique() %>% length()
  #cluster_summary$Company[i] <- clusters %>% filter(cluster == i) %>% pull(Company) %>% unique()
  
  qids <- clusters %>% filter(cluster == i) %>% pull(ID)
  qdists <- dists_25snps %>% filter(ID1 %in% qids | ID2 %in% qids) %>% pull(Dist)
  
  cluster_summary$Mean_dist[i] <- mean(qdists)
  cluster_summary$Median_dist[i] <- median(qdists)
  cluster_summary$Ncol1[i] <- clusters %>% filter(cluster == i, visit == 1) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$Ncol2[i] <- clusters %>% filter(cluster == i, visit == 2) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$Ncol3[i] <- clusters %>% filter(cluster == i, visit == 3) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$Ncol4[i] <- clusters %>% filter(cluster == i, visit == 4) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$Ncol5[i] <- clusters %>% filter(cluster == i, visit == 5) %>% pull(patient) %>% unique() %>% length()
  cluster_summary$Ncol6[i] <- clusters %>% filter(cluster == i, visit == 6) %>% pull(patient) %>% unique() %>% length()
}


write.csv(cluster_summary, "Transmission_clusters_summary.csv", row.names = FALSE)


wilcox.test(N_isolates ~ ST, data = cluster_summary) %>% summary()
wilcox.test(Mean_dist ~ ST, data = cluster_summary) 

aggregate(cluster_summary$N_isolates, list(cluster_summary$ST), FUN=mean) 
aggregate(cluster_summary$N_isolates, list(cluster_summary$ST), FUN=median) 

aggregate(cluster_summary$N_patients, list(cluster_summary$ST), FUN=mean) 
aggregate(cluster_summary$N_patients, list(cluster_summary$ST), FUN=median) 

aggregate(cluster_summary$Mean_dist, list(cluster_summary$ST), FUN=mean)

## Tabulate proportion of sites in ST5 and ST8 clusters
tab <-  table(clusters$ST,clusters$Site)
prop.table(tab, margin = 1)
chisq.test(tab)

fulltab <-  table(metadata$ST,metadata$Site)
prop.table(fulltab, margin = 1)
chisq.test(fulltab)
fisher.test(fulltab)

boxplot(N_isolates ~ ST, data = cluster_summary)

#### Plot transmission network


links9 <- dists %>% filter(Dist < 10, Same_patient== FALSE,ID1 != "CP000730.1") %>% dplyr::select(ID1,ID2,Dist)

link9_IDs <- c(links9$ID1,links9$ID2) %>% unique()
links9_metadata <- metadata %>% filter(ID %in% link9_IDs)
links9_metadata$patient <- substr(links9_metadata$ID, 1,4)
### 1149.6C is duplicated index number 19
network9 <- graph_from_data_frame(links9, directed=F, vertices = links9_metadata) 

# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(4, "Set1") 

# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network9)$Platoon))]
V(network9)$shape <- ifelse(V(network9)$Company == "A150", "square", "circle")

s.factor(V(network9)$Company)

c2 = cluster_leading_eigen(network9) 
coords = layout_with_fr(network9)
# Make the plot
#plot(c2, network9, layout=coords , vertex.color=my_color, vertex.size = 7, vertex.shapes = my_shape, vertex.label = V(network9)$Subject, edge.weight = E(network9)$Dist)
plot(network9, vertex.color=my_color, vertex.size = 3, vertex.shapes = V(network9)$shape, vertex.label = V(network9)$Patient, edge.weight = E(network9)$Dist)
legend("bottomleft", legend=levels(as.factor(V(network9)$Platoon)) , title = "Platoon" , col = coul , bty = "n", pch=20 , pt.cex = 5, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
legend("bottomleft", legend=c("A150","C219"), title = "Company" ,bty = "n", pch=c(15,16) , pt.cex = 5, cex = 1.5, text.col="black" , horiz = FALSE, inset = c(0, 0.1))


components(network9, mode = c("weak", "strong"))  










### Identify patients where there seems to be an importation

Importers <- within_host %>% filter(Dist > 20) %>% pull(Patient1) %>% unique()
Importers <- data.frame(Importers)

for(i in 1:nrow(Importers)){
  query <- Importers$Importers[i]
  IDs <- within_host %>% dplyr::filter(Patient1 == query) %>% pull(ID1)
  Importers$index[i] <- IDs[1]
  
  Outliers <- within_host %>% filter(ID1 == IDs[1],Dist>20)
  Importers$Count[i] <- nrow(Outliers)
  Importers$Importees[i] <- (within_host %>% filter(ID1 == IDs[1],Dist>20) %>% pull(ID2))[1]
  
}


for(i in 1:nrow(Importers)){
  #ID <- Importers$Importers[i]
  Importers$Mean_Dist[i] <- within_host %>% filter(Patient1 == Importers$Importers[i]) %>% pull(Dist) %>% mean()
  Importers$CountST <- within_host %>% filter(Patient1 == Importers$Importers[i]) %>% pull(ST1) %>% unique() %>% length()
  }
### Re-run the parallels and importation steps on the entire within_host dataframe
### Steps:
### 1) Identify the index strain for each patient
### 2) Identify outliers as strains that are more the 12 SNPs apart
### 3) Count number of outliers
### 3) For each outlier, count the number of potential transmission pairs from other patients

Patients_with_multiple <- within_host %>% pull(Patient1) %>% unique()

Outlier_df <- data.frame(Patients_with_multiple)
Outlier_df$patient <- within_host %>% pull(Patient1) %>% unique()

for(i in 1: nrow(Outlier_df)){
  patient <- Outlier_df$Patients_with_multiple[i]
  IDs <- within_host %>% dplyr::filter(Patient1 == patient) %>% pull(ID1)
  Outlier_df$Index[i] <- IDs[1]
  Outlier_df$Outliers[i] <- within_host %>% dplyr::filter(ID1 == IDs[1], Dist > 12) %>% nrow()
  Outers <- within_host %>% dplyr::filter(ID1 == IDs[1], Dist > 12) %>% pull(ID2)
  Outlier_df$Outer1[i] <- Outers[1]
  Outlier_df$Outer2[i] <- Outers[2]
}   

### Pairs with more than two outliers
# 1148 and 2061 1148.1O both cases second and 3rd are identical

All_index_IDs <- c(Outlier_df$Index, Outlier_df$Outer1) 

All_index_IDs <- All_index_IDs[!is.na(All_index_IDs)]


## Sub_sample putative transmissions between indexes
Index_Transmissions <- dists %>% filter(Dist < 25, Same_patient == FALSE, ID1 %in% All_index_IDs, ID2 %in% All_index_IDs)

Index_t_tbl <- gather(Index_Transmissions, Demographic, Outcome, Same_patient:Same_platoon, factor_key=TRUE)

Index_t_tbl %>% filter(Demographic %in% c("Same_Company", "Same_platoon")) %>% ggplot(aes(x = Demographic, fill = Outcome)) + geom_bar()

table(Index_Transmissions$Company1)

Index_Transmissions_50 <- dists %>% filter(Dist < 50, Same_patient == FALSE, ID1 %in% All_index_IDs, ID2 %in% All_index_IDs)

Index_t_tbl_50 <- gather(Index_Transmissions_50, Demographic, Outcome, Same_patient:Same_platoon, factor_key=TRUE)

Index_t_tbl_50 %>% filter(Demographic %in% c("Same_Company", "Same_platoon")) %>% ggplot(aes(x = Demographic, fill = Outcome)) + geom_bar()



Index_Transmissions_50 %>% ggplot(aes(x = Dist, fill = Same_platoon)) + geom_bar()


Index_Transmissions_50 %>% ggplot(aes(x = Dist, fill = Same_platoon)) + geom_bar(stat="count", position ="fill", bins = 10)

Index_Transmissions %>% ggplot(aes(x = Dist, fill = Same_platoon, group = Same_platoon)) + geom_density(alpha = 0.5, bw = 1)


#### Identify potential donors for each Index strain

# 1) Import metadata and mlst Profies and merge

metadata <- read.csv("2015_Staph_metadata.csv") 
mlst <- read.csv("2015_MLST_profiles.csv")

metadata <- merge(metadata,mlst,by="ID")
metadata$Deci_date <- decimal_date(as.Date(metadata$Date))


metadata %>% filter(Subject %in% Patients_with_multiple) %>% ggplot(aes(x = Visit, y = Collection_Site, fill = ST)) + 
  geom_tile() + facet_wrap( ~ Subject, ncol = 8) + theme_bw()

metadata %>% filter(Subject %in% Patients_with_multiple) %>% ggplot(aes(x = Date, y = Collection_Site, fill = ST)) + 
  geom_tile() + facet_wrap( ~ Subject, ncol = 5) + theme_bw() + theme(axis.text.x = element_text(angle = 90))

c <- 0
for (a in Patients_with_multiple){
  c <- c+1
  #patient <- a
  #IDs <- within_host %>% dplyr::filter(Patient1 == patient) %>% pull(ID1)
  #index <- 
  #out_IDs <- 
  count_out <- length(out_IDs)
  Outlier_df$Patient[c] <- a
  Outlier_df$Index[c] <- IDs[1]
  Outlier_df$Count_outliers[c] <- within_host %>% dplyr::filter(ID1 == index, Dist > 12) %>% pull(ID2) %>% length()
}



#Check imported samples in transmissions
for (i in 1:nrow(Importers)){
  Index <- Importers$index[i]
  Import <- Importers$Importees[i]
  Importers$transmission_pairs[i] <- Transmissions %>% filter(ID2 == Import || ID2 == Import) %>% nrow()
}


# levels(trans_dists$ID1[i])[trans_dists$ID1[i]]





#### Plot transmission network


links9 <- Index_Transmissions %>% filter(Dist < 10) %>% dplyr::select(ID1,ID2,Dist)

link9_IDs <- c(links9$ID1,links9$ID2) %>% unique()
links9_metadata <- metadata %>% filter(ID %in% link9_IDs) 
### 1149.6C is duplicated index number 19
links9_metadata <- links9_metadata[-19,]
network9 <- graph_from_data_frame(links9, directed=F, vertices = links9_metadata) 

# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(5, "Set1") 

# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network9)$Platoon))]
V(network9)$shape <- ifelse(V(network9)$Company == "A150", "square", "circle")

  s.factor(V(network9)$Company)
  
  c2 = cluster_leading_eigen(network9) 
  coords = layout_with_fr(network9)
# Make the plot
#plot(c2, network9, layout=coords , vertex.color=my_color, vertex.size = 7, vertex.shapes = my_shape, vertex.label = V(network9)$Subject, edge.weight = E(network9)$Dist)
  plot(network9, vertex.color=my_color, vertex.size = 7, vertex.shapes = my_shape, vertex.label = V(network9)$ST, edge.weight = E(network9)$Dist)
  legend("bottomleft", legend=levels(as.factor(V(network9)$Platoon)) , title = "Platoon" , col = coul , bty = "n", pch=20 , pt.cex = 5, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
  legend("bottomleft", legend=c("A150","C219"), title = "Company" ,bty = "n", pch=c(15,16) , pt.cex = 5, cex = 1.5, text.col="black" , horiz = FALSE, inset = c(0, 0.1))

  
  components(network9, mode = c("weak", "strong"))  
  
  
  
## vertex.label = nodes[c(as.numeric(V(gg)$name)+1),]$name
plot(network, vertex.size = 5, vertex.shapes = Platoon)

links12 <- Index_Transmissions %>% filter(Dist < 13) %>% dplyr::select(ID1,ID2,Dist)
network12 <- graph_from_data_frame(d=links12, directed=F) 
plot(network, vertex.size = 5)

E(network)$weight <- links$Dist

%>% ggraph(layout = "fr")



actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(g)
