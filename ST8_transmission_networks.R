library("ggplot2")
library("tidyverse")
library(lubridate)
library("ggraph")
library(igraph)
library(RColorBrewer)
library(beeswarm)
library(cowplot)

setwd("~/Projects/Staph_military/Reanalysis/Distance_analysis/ST8/")

dists <- read.csv("ST8_distances_metadata.csv")

metadata <- read.csv("../../Metadata/MLST_profiles_metadata.csv")

within_host <- dists %>% filter(Same_patient == TRUE) %>% filter(Dist < 1000)

same_host_site <- dists %>% dplyr::filter(Same_patient == TRUE, Same_site == TRUE, Dist < 1000) 

summary(within_host$Dist)
quantile(within_host$Dist, 0.95)


#### Build a transmission network for ST8
ST8_meta <- metadata %>% filter(ST =="8")
ST8_meta$patient <- substr(ST8_meta$ID,1,4)

 
## Pick an index strain for each patient and compar others to see if there is overlap
PatientsST8 <- c(dists$Patient1,dists$Patient2) %>% unique()

Outlier_df <- data.frame(PatientsST8)

for(i in 1: nrow(Outlier_df)){
    patient <- Outlier_df$PatientsST8[i]
  IDs <- within_host %>% dplyr::filter(Patient1 == patient) %>% pull(ID1)
  Outlier_df$Index[i] <- IDs[1]
  Outlier_df$Outliers[i] <- within_host %>% dplyr::filter(ID1 == IDs[1], Dist > 15) %>% nrow()
  Outers <- within_host %>% dplyr::filter(ID1 == IDs[1], Dist > 15) %>% pull(ID2)
  Outlier_df$Outer1[i] <- Outers[1]
  Outlier_df$Outer2[i] <- Outers[2]
}  

ST8_indexes <- c(Outlier_df$Index,Outlier_df$Outer1) 
ST8_indexes <- ST8_indexes[!is.na(ST8_indexes)]  


### Build a network for all ST isolates ###########
ST8_all_transmissions <- dists %>% filter(Dist < 25, Same_patient == FALSE) 
All_ST8_ID <- c(ST8_all_transmissions$ID1,ST8_all_transmissions$ID2) %>% unique()
ST8_all_links_metadata <- ST8_meta %>% filter(ID %in% All_ST8_ID)
ST8_all_links_metadata$patient  <- substr(ST8_all_links_metadata$ID, 1, 4)

full_networkST8 <- graph_from_data_frame(ST8_all_transmissions, directed=F, vertices = ST8_all_links_metadata)
plot(full_networkST8)
components(full_networkST8, mode = c("weak", "strong"))  

ST8_all_clusters <- components(full_networkST8, mode = c("weak", "strong"))  

ST8_all_cluster_index <- ST8_all_clusters$membership %>% as.data.frame()


###### Identify which network each isolate belongs to and plot a tile plot for each cluster
names <- rownames(ST8_all_cluster_index)
rownames(ST8_all_cluster_index) <- NULL
ST8_all_cluster_index <- cbind(names,ST8_all_cluster_index)
colnames(ST8_all_cluster_index) <- c("ID","Cluster")

### update the ST8 metadata table with the cluster information 

for(i in 1: nrow(ST8_meta)){
  ST8_meta$In_cluster[i] <- ST8_meta$ID[i] %in% names
  if(ST8_meta$In_cluster[i] == TRUE){
    q <- which(ST8_all_cluster_index$ID == ST8_meta$ID[i])
    ST8_meta$Cluster[i] <- ST8_all_cluster_index$Cluster[q]
  }
}

##### Create a matrix for all strains where they are paired with all possible indexes
source_df <- expand.grid(ST8_meta$ID,ST8_indexes)
colnames(source_df) <- c("Query", "Donor")

## Get cluster information for each Query and Donor
for(i in 1:nrow(source_df)){
  q <- which(ST8_meta$ID == source_df$Query[i])
  d <- which(ST8_meta$ID == source_df$Donor[i])
  source_df$Q_cluster[i] <- ST8_meta$Cluster[q]
  source_df$D_cluster[i] <- ST8_meta$Cluster[d]
}


source_df <- source_df %>% filter(Q_cluster == D_cluster, Q_cluster !="")

#### Extract on each donor from the index strains
for(i in 1:nrow(source_df)){
  donor <- levels(source_df$Donor[i])[source_df$Donor[i]]
  query <- levels(source_df$Query[i])[source_df$Query[i]]
  donor_cloud <- within_host %>% dplyr::filter(ID1 == donor | ID2 == donor, Dist < 16)
  cloud_IDs <- c(donor,donor_cloud$ID1,donor_cloud$ID2) %>% unique()
  cloud_visits <- c(donor_cloud$Visit1,donor_cloud$Visit2)
  comparisons <- dists %>% filter((ID1 == query & ID2 %in% cloud_IDs) | (ID2 == query & ID1 %in% cloud_IDs))
  
  source_df$Cloud_size[i] <- length(cloud_IDs)
  source_df$Cloud_duration[i] <- max(cloud_visits) - min(cloud_visits) + 1
  source_df$Cloud_min_snps[i] <- min(comparisons$Dist)
  source_df$Cloud_max_snps[i] <- max(comparisons$Dist)
  source_df$Cloud_mean_snps[i] <- mean(comparisons$Dist)
}

source_df$Self <- substr(source_df$Query,1,4) == substr(source_df$Donor,1,4)

### Loop through each ID and determine the likelest source based on mean snp diff
query_IDs <- source_df %>% pull(Query) %>% unique()

Query_hits <- data.frame(query_IDs)
tophits <- data.frame()
top2hits <- data.frame()
for (i in 1:nrow(Query_hits)){
  q <- levels(Query_hits$query_IDs[i])[Query_hits$query_IDs[i]]
  tmp <- source_df[order(source_df$Cloud_mean_snps) , ] %>% filter(Query == q)
  tophits<- rbind(tophits,tmp[1,])
  top2hits <- rbind(tophits,tmp[1,],tmp[2,])
}


ST8_link_IDs <- c(ST8_transmissions$ID1,ST8_transmissions$ID2) %>% unique()
ST8_links_metadata <- ST8_meta %>% filter(ID %in% ST8_link_IDs) 

network <- graph_from_data_frame(source_df, directed=F) 
plot(network,vertex.size = 7)

write.csv(source_df,"Potential_sources_for_each_strain.csv")


mtcars[order(mpg),]

1035_3I

within_host %>% dplyr::filter(ID1 == "1035_3I" | ID2 == "1035_3I", Dist > 15)
levels(trans_dists$ID1[i])[trans_dists$ID1[i]]


ST8_clusters_metadata <- merge(ST8_all_cluster_index,ST8_all_links_metadata, by = "ID" )

ST8_clusters_metadata  %>% ggplot(aes(x = Date, y = Site, fill = Site)) +
  geom_tile() + facet_grid(Cluster~patient) 



Plot_cluster_tile <- function(i){
  plotname <- paste("Plot for cluster",i,sep = " ")
  plot_tile <- ST8_clusters_metadata %>% filter(Cluster == i) %>% ggplot(aes(x = as.integer(visit), y = Site)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, ncol = 1) + 
  xlab("Visit Number") + ylab("") + theme_bw(base_size=16, base_family='Times New Roman') + xlim(1,6) +
    scale_x_continuous(breaks=seq(1,6,1))
  return(plot_tile)
}

pdf("plots.pdf", width=9, height=12)
cluster_patient_count <- for(a in 1:11){
  n_patients <- ST8_clusters_metadata %>% filter(Cluster == a) %>% pull(patient) %>% unique() %>% length()
  #tiles <- Plot_cluster_tile(a)
  #plot(tiles)
 return(n_patients)
  }
dev.off()

tile_plot1	<-	Plot_cluster_tile(1)
tile_plot2	<-	 Plot_cluster_tile(2)
tile_plot3	<-	 Plot_cluster_tile(3)
tile_plot4	<-	 Plot_cluster_tile(4)
tile_plot5	<-	 Plot_cluster_tile(5)
tile_plot6	<-	 Plot_cluster_tile(6)
tile_plot7	<-	 Plot_cluster_tile(7)
tile_plot8	<-	 Plot_cluster_tile(8)
tile_plot9	<-	 Plot_cluster_tile(9)
tile_plot10	<-	 Plot_cluster_tile(10)
tile_plot11	<-	 Plot_cluster_tile(11) 

ST8_clusters_metadata %>% filter(Cluster == 1) %>% ggplot(aes(x = as.integer(visit), y = Site)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, nrow = 1) + 
  xlab("Visit Number") + ylab("") + theme_bw(base_size=16, base_family='Times New Roman') + xlim(1,6) +
  scale_x_continuous(breaks=seq(1,6,1)) 


tile_plot4 <- scale_fill_manual(breaks = c("2", "1", "0.5"), 
                  values=c("red", "blue", "green"))

plot_grid(tile_plot2,tile_plot7,tile_plot11, ncol = 3) + ylab("Colonization Site") + xlab("Date of Isolation")
plot_grid(tile_plot6,tile_plot9, ncol = 2)

cluster_n <- c(8,2,7,5,4,3,2,6,3,3,2)

plot_grid(tile_plot1,tile_plot2, nrow = 2, rel_widths = c(8,2))




for (i in 1:11){
  fileout <- paste("tile_plot_v2_",i,".pdf", sep = "")
  #plot_the_tile <- paste("plot(tile_plot",i,")",sep = "")
  plotname
  size <- cluster_n[i]
  Plot_cluster_tile(i)
  # 2.1. Save the plot to a pdf
  if(size <4 ){ggsave(fileout, width = 3*size + 2, height = 3)}
  else if(4 <= size < 7){ggsave(fileout, width = 3*size +2, height = 6)} 
  else{ggsave(fileout, width = 3*size +2, height = 9)}
}


### Generate tile plots for isolates that have more than one colonizing ST
### ST8 and ST5, patients: 1004, 1148, 2013
### Other pairs: 1143, 2130, 
### Only two isolates: 2031, 2069
metadata$visit <- substr(metadata$ID,6,6)

metadata %>%  filter(patient %in% c("1004","1148","2013","1134","2130")) %>% 
  ggplot(aes(x = as.integer(visit), y = Site, fill = ST)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, nrow = 2) + 
  xlab("Visit Number") + ylab("") + theme_bw(base_size=16, base_family='Times New Roman') 

metadata %>%  filter(patient %in% c("1004","1148","2013")) %>% 
  ggplot(aes(x = as.integer(visit), y = Site, fill = ST)) +
  geom_tile(alpha = 0.5) + facet_wrap(~patient, nrow = 1) + 
  xlab("Visit Number") + ylab("") + theme_bw(base_size=16, base_family='Times New Roman') 


ST5_sites <- metadata %>% filter(ST ==5) %>% ggplot( aes(x=Site)) + geom_bar(aes( y = ..count../sum(..count..), fill = Site)) + ylab("Proportion") 
ST8_sites <- metadata %>% filter(ST ==8) %>% ggplot( aes(x=Site)) + geom_bar(aes( y = ..count../sum(..count..), fill = Site)) + ylab("Proportion") 
All_ST_sites <- metadata %>% ggplot( aes(x=Site)) + geom_bar(aes( y = ..count../sum(..count..), fill = Site)) + ylab("Proportion") + main("All isolates")

metadata %>% filter(ST %in% c(5,8)) %>% ggplot( aes(x=Site)) + geom_bar(aes( y = ..count../sum(..count..), fill = Site)) + ylab("Proportion") +
                      facet_wrap(~ST, nrow = 2)
                    

plot_grid(ST5_sites,ST8_sites, nrow = 2)
  opts(legend.position = 'none') +
  opts(axis.text.y = theme_blank(), axis.title.y = theme_blank()) + 
  opts(title = 'Male', plot.title = theme_text( size = 10) ) +  
  coord_flip()
   

  metadata$Oropharyinx <- metadata$Site == "Oropharynx"
tbl <- table(metadata$ST, metadata$Oropharyinx)
chisq.test(tbl)

metadata$ST5 <- metadata$ST == 5

tbl2 <- table(metadata$ST5, metadata$Oropharyinx)
chisq.test(tbl2)


###########################3 old code

ST8_transmissions <- dists %>% filter(ID1 %in% ST8_indexes , ID2 %in% ST8_indexes, Dist < 26)  

#links9 <- Index_Transmissions %>% filter(Dist < 10) %>% dplyr::select(ID1,ID2,Dist)

ST8_link_IDs <- c(ST8_transmissions$ID1,ST8_transmissions$ID2) %>% unique()
ST8_links_metadata <- ST8_meta %>% filter(ID %in% ST8_link_IDs) 

networkST8 <- graph_from_data_frame(ST8_transmissions, directed=F, vertices = ST8_links_metadata) 


# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(4, "Set1") 

# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(networkST8)$Platoon))]
V(networkST8)$shape <- ifelse(V(networkST8)$Company == "A150", "square", "circle")
as.factor(V(networkST8)$Company)

c2 = cluster_leading_eigen(networkST8) 
coords = layout_with_fr(networkST8)
# Make the plot
#plot(c2, networkST8, layout=coords , vertex.color=my_color, vertex.size = 7, vertex.shapes = my_shape, vertex.label = V(networkST8)$Subject, edge.weight = E(networkST8)$Dist)
plot(networkST8, vertex.color=my_color, vertex.size = 7, vertex.shapes = V(networkST8)$shape, vertex.label = V(networkST8)$patient, edge.weight = E(networkST8)$Dist, vertex.label.dist = 1)
legend("bottomleft", legend=levels(as.factor(V(networkST8)$Platoon)) , title = "Platoon" , col = coul , bty = "n", pch=20 , pt.cex = 2.5, cex = 1, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
legend("bottomleft", legend=c("A150","C219"), title = "Company" ,bty = "n", pch=c(15,16) , pt.cex = 2.5, cex = 1, text.col="black" , horiz = FALSE, inset = c(0, 0.1))


components(networkST8, mode = c("weak", "strong"))  

components(networkST8, mode = c("weak", "strong"))  

ST8_clusters <- components(networkST8, mode = c("weak", "strong"))  

ST8_cluster_index <- ST8_clusters$membership %>% as.data.frame()

names <- rownames(ST8_cluster_index)
rownames(ST8_cluster_index) <- NULL
ST8_clusters <- cbind(names,ST8_cluster_index)

colnames(ST8_clusters) <- c("ID","Cluster")
ST8_clusters$patient <- substr(ST8_clusters$ID,1,4)
clusters <- ST8_clusters %>% pull(Cluster) %>% unique()

for(i in clusters){
  plotname <- paste("Plot for cluster",i,sep = " ")
  patients <- ST8_clusters %>% filter(Cluster == i) %>% pull(patient)
  plot_tile <- ST8_meta %>% filter(patient %in% patients) %>% ggplot(aes(x = as.Date(Date), y = Site, fill = Site)) +
    geom_tile(alpha = 0.5) + facet_wrap(~patient, ncol = 3)
  plot(plot_tile)
}

patients1 <- ST8_clusters %>% filter(Cluster == 1) %>% pull(patient)
plot_tile1 <- ST8_meta %>% filter(patient %in% patients1) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 1") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile1)

patients2 <- ST8_clusters %>% filter(Cluster == 2) %>% pull(patient)
plot_tile2 <- ST8_meta %>% filter(patient %in% patients2) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 2") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile2)

patients3 <- ST8_clusters %>% filter(Cluster == 3) %>% pull(patient)
plot_tile3 <- ST8_meta %>% filter(patient %in% patients3) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 3") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile3)


patients4 <- ST8_clusters %>% filter(Cluster == 4) %>% pull(patient)
plot_tile4 <- ST8_meta %>% filter(patient %in% patients4) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 4") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile4)

patients5 <- ST8_clusters %>% filter(Cluster == 5) %>% pull(patient)
plot_tile5 <- ST8_meta %>% filter(patient %in% patients5) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 5") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile5)

patients6 <- ST8_clusters %>% filter(Cluster == 6) %>% pull(patient)
plot_tile6 <- ST8_meta %>% filter(patient %in% patients6) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 6") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile6)

patients7 <- ST8_clusters %>% filter(Cluster == 7) %>% pull(patient)
plot_tile7 <- ST8_meta %>% filter(patient %in% patients7) %>% ggplot(aes(x = (as.Date(Date, "%m/%d/%Y")), y = Site, fill = Site)) +
  geom_tile() + facet_wrap(~patient, ncol = 1) + ggtitle("Cluster 7") + theme(axis.text.x = element_text(angle = 90)) +
  scale_x_date(date_labels = "%b %d") + xlab("Date of sampling") + ylab("Colonization Site")
plot(plot_tile7)
