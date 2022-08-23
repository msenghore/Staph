library("tidyverse")
library(lubridate)

setwd("~/Projects/Staph_military/Reanalysis/Metadata/")

metadata <- read.csv("MLST_profiles_metadata.csv")
micro <- read.csv("2015_All_Staph_dates_microbiology.csv")
molecular <- read.csv("2015_molecular_typing.csv")

meta_IDs <- metadata %>% pull(ID)

meta_IDs[which(meta_IDs %in% micro$ID)]

metadata <- merge(metadata,micro, by="ID", all = T)

table(metadata$Species)

metadata <- metadata %>% filter(Species != "NA")

metadata <- merge(metadata,molecular, by="ID")



write.csv(metadata, "Staph_2015_mlst_lab_metadata.csv")
