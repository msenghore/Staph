library("ggplot2")
library("tidyverse")

setwd("~/Projects/Staph_military/2015_alignment/")

dists <- read.csv("MRSA_2015_snps_dist.csv")

within_host <- dists %>% dplyr::filter(Same_patient == TRUE) %>% dplyr::filter(SNP_dist < 1000)

same_host_site <- dists %>% dplyr::filter(Same_patient == TRUE, Same_site == TRUE, SNP_dist < 1000) 

summary(within_host$SNP_dist)
quantile(within_host$SNP_dist, 0.95)\  


ggplot(data = within_host, aes( x = as.numeric(Days), y = SNP_dist)) + geom_point(aes(fill = Same_site, color = Same_site)) + geom_smooth(method=lm)


ggplot(data = same_host_site, aes( x = as.numeric(Days), y = SNP_dist)) + geom_point() + geom_smooth(method=lm)

 
plot(within_host$SNP_Dist , within_host$Days)
 abline(lm(within_host$SNP_Dist ~ within_host$Days + within_host$Same_site), col = "red")

 
 ##3 Run a logistic regression
 SNP_rate <- lm(SNP_dist ~ as.numeric(Days) + Same_site, data = within_host)
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

##Plot time series

Parallel_dists %>% ggplot(aes(x = Days, y = SNP_dist, shape = Site2, color = Visit2, group = Site2, label = Site2)) + 
     geom_point() + geom_smooth() + facet_wrap(~ID1) + ylim(0,10) + geom_text()

Parallel_dists %>% ggplot(aes(x = Days, y = SNP_dist, shape = Site2, color = Visit2, group = Site2)) + 
  geom_point() + geom_smooth() + facet_wrap(~ID1, scales = "free") 

### Identify patients where there seems to be an importation

Importers <- within_host %>% filter(SNP_dist > 20) %>% pull(Patient1) %>% unique()
Importers <- data.frame(Importers)

for(i in 1:nrow(Importers)){
  query <- Importers$Importers[i]
  IDs <- within_host %>% dplyr::filter(Patient1 == query) %>% pull(ID1)
  Importers$index[i] <- IDs[1]
  
  Outliers <- within_host %>% filter(ID1 == IDs[1],SNP_dist>20)
  Importers$Count[i] <- nrow(Outliers)
  Importers$Importees[i] <- (within_host %>% filter(ID1 == IDs[1],SNP_dist>20) %>% pull(ID2))[1]
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
  Outlier_df$Outliers[i] <- within_host %>% dplyr::filter(ID1 == IDs[1], SNP_dist > 12) %>% nrow()
  Outers <- within_host %>% dplyr::filter(ID1 == IDs[1], SNP_dist > 12) %>% pull(ID2)
  Outlier_df$Outer1[i] <- Outers[1]
  Outlier_df$Outer2[i] <- Outers[2]
}   

### Pairs with more than two outliers
# 1148 and 2061 1148.1O both cases second and 3rd are identical

All_index_IDs <- c(Outlier_df$Index, Outlier_df$Outer1) 

All_index_IDs <- All_index_IDs[!is.na(All_index_IDs)]


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
  Outlier_df$Count_outliers[c] <- within_host %>% dplyr::filter(ID1 == index, SNP_dist > 12) %>% pull(ID2) %>% length()
}

## Sub_sample putative transmissions between patients
Transmissions <- dists %>% filter(SNP_dist < 25, Same_patient == FALSE)

#Check imported samples in transmissions
for (i in 1:nrow(Importers)){
  Index <- Importers$index[i]
  Import <- Importers$Importees[i]
  Importers$transmission_pairs[i] <- Transmissions %>% filter(ID2 == Import || ID2 == Import) %>% nrow()
}


# levels(trans_dists$ID1[i])[trans_dists$ID1[i]]





