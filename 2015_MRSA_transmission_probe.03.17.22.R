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

metadata <- read.csv("../Metadata/MLST_profiles_metadata.csv")

within_host <- dists %>% filter(Same_patient == TRUE) %>% filter(Dist < 1000)

same_host_site <- dists %>% dplyr::filter(Same_patient == TRUE, Same_site == TRUE, Dist < 1000) 

summary(within_host$Dist)
quantile(within_host$Dist, 0.95)

### Plot box plot of within host diversity for ST 5,8 and 87

# Within_plot <- 
  
  dists %>% filter(ST1 %in% c("5","8","87"), ST1 == ST2) %>% ggplot(aes(x = Same_patient, y = Dist, fill = ST1)) +
    geom_boxplot() + scale_x_discrete(labels=c("FALSE" = "Between hosts", "TRUE" = "Within host")) + scale_y_continuous(trans='log10') +
    xlab("") + ylab("Pairwise SNP distance") + scale_fill_manual(values=c("#B9BFFF", "#F4B6B2", "#d8b365")) + 
    ggtitle("Pairwise SNP distances by Sequence Type (ST)") + theme(plot.title = element_text(hjust = 0.5))
  
  #scale_fill_brewer( type = "qual")



ggplot(data = within_host, aes( x = as.numeric(Days), y = Dist)) + geom_point(aes(fill = Same_site, color = Same_site)) + geom_smooth(method=lm)

lm(Dist ~ as.numeric(Days) + Same_site, within_host)


ggplot(data = same_host_site, aes( x = as.numeric(Days), y = Dist)) + geom_point() + geom_smooth(method=lm)

 
plot(within_host$Dist , within_host$Days)
 abline(lm(within_host$Dist ~ within_host$Days + within_host$Same_site), col = "red")

 
 ##3 Run a logistic regression
 SNP_rate <- lm(Dist ~ as.numeric(Days) + Same_site, data = within_host)
summary(SNP_rate) 

## Subset patients that had at least 5 isolates sequenced
Parallels <- c(2011,1142,1113,2017,2075,2189,2014,2040,2092,1041,2010,1004,1054,1148,1173,2130,2159,2166,2183,2209)
P_evo <- data.frame(Parallels)

## For each of these patients, identify the index case (first sample, subject ie to alphabetic order of specimen sites)
for(i in 1:nrow(P_evo)){
  query <- P_evo$Parallels[i]
  IDs <- within_host %>% dplyr::filter(Patient1 == query) %>% pull(ID1)
  P_evo$index[i] <- IDs[1]
  }

## Filter sites so the index strain is in site snp dists file
Parallel_dists <- filter(within_host, ID1 %in% P_evo$index)


Parallel_dists$
##Plot time series

Parallel_Evolution_plot <- Parallel_dists %>% ggplot(aes(x = as.numeric(Visit), y = Dist, shape = Site2, color = ST2, label = Site2)) + 
     geom_point() + facet_wrap(~ID1, scales = "free") + ylim(0,100) + xlab("Days since index strain") + ylab("SNP distance from index strain")

Parallel_dists %>% ggplot(aes(x = Days, y = Dist, shape = Site2, color = Visit2, group = Site2)) + 
  geom_point() + geom_smooth() + facet_wrap(~ID1, scales = "free", color = ST1) 


### Plot distance
### Counting transmission pairs

k_snps <- c(seq(1,10,1),seq(15,50,5),seq(60,100,10))

Ksnp_comparisons <- as.data.frame(k_snps)

for(i in 1: nrow(Ksnp_comparisons)){
  kmer <- Ksnp_comparisons$k_snps[i]
  Ksnp_comparisons$All_comparisons[i] <- dists %>% filter(Dist <= kmer) %>% nrow()
  Ksnp_comparisons$Within_host[i] <- dists %>% filter(Dist <= kmer, Same_patient == TRUE) %>% nrow()
  Ksnp_comparisons$All_pairs[i] <- dists %>% filter(Dist <= kmer, Same_patient == FALSE) %>% nrow()
  #Ksnp_comparisons$Unique_pairs[i] <- dists %>% filter(Dist <= kmer, Same_patient == FALSE) %>% filter(ID1 %in% All_index_IDs & ID2 %in% All_index_IDs) %>% nrow()
}

ksnp_tbl <- gather(Ksnp_comparisons,Category,Count, All_comparisons:All_pairs, factor_key=TRUE)

ksnp_tbl %>% ggplot(aes(x = k_snps, y = Count,)) + geom_smooth() + facet_wrap( ~ Category, ncol = 3, scales = "free")

ksnp_tbl %>% filter(k_snps < 26) %>% ggplot(aes(x = k_snps, y = Count,)) + geom_smooth() + facet_wrap( ~ Category, ncol = 3, scales = "free") +
  ylab("Cumulative number of comparisons") + xlab("SNP threshold")

### Plot density plots for SNP distance

dists %>% filter (Dist < 101) %>% ggplot(aes(x = Dist, fill = Same_patient, color = Same_patient,  bw = 2)) + geom_density(alpha = 0.5) 


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

dists_20snps <- filter(dists, Dist < 21, Same_patient == FALSE)
network20 <- graph_from_data_frame(d=dists_20snps, directed=F) 
plot(network20, vertex.size = 5)

### Plot a network with metadata
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
