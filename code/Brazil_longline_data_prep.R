### ##########################################
###
### Brazilian longline bycatch analysis
### DATA PREPARATION
### based on scripts from Namibia, Da Rocha et al. 2021: https://www.sciencedirect.com/science/article/abs/pii/S0006320720309733
###
### ##########################################

### adapted up by Steffen Oppel 31 July 2019
### data provided by Yann Rouxel and Dimas Gianuca in single xlsx spreadsheet



##############################################################
#### load ALL NECESSARY libraries
##############################################################
setwd("C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL")
library(tidyverse)
library(readxl)
library(lubridate)
library(dplyr)
library(data.table)
library(maps)


##############################################################
#### READ AND MANIPULATE ATF Fisheries observer data
##############################################################

#### read in original file and rename columns
## coordinates are in awful fucked up format with or without decimal point and with range of decimal places
## manually fixed two dates in Excel file (Set 742 and 918)
dd <- read_excel("data/Brazil_longline_bycatch_data.xlsx",   
                 sheet="DB_Observer_Full_ver00") %>%
  rename(Set=`Set#`,BYCATCH=`N birds`) %>%
  mutate(Longitude=ifelse(Longitude>-180,Longitude,Longitude/10000)) %>%
  mutate(Longitude=ifelse(Longitude>-180,Longitude,Longitude/10)) %>%
  mutate(Longitude=ifelse(Longitude>-10,Longitude*10,Longitude)) %>%
  mutate(Latitude=ifelse(Latitude>-90,Latitude,Latitude/10000)) %>%
  mutate(Latitude=ifelse(Latitude>-180,Latitude,Latitude/10)) %>%
  mutate(Latitude=ifelse(Latitude>-10,Latitude*10,Latitude)) %>%
  mutate(Moon.il=ifelse(Moon.il>1,Moon.il/10,Moon.il)) %>%   ### 5 sets have moon.ill >1 which should not be possible
  mutate(Date=as.Date(Date, origin = "1899-12-30"))

summary(dd)


##############################################################
#### CHECK THAT COORDINATES MAKE SENSE
##############################################################

map(database = "world") 
points(x = dd$Latitude, y = dd$Longitude, col="firebrick", pch=21)


##############################################################
#### SAVE PROCESSED DATA
##############################################################
saveRDS(dd,"data/Brazil_formatted_bycatch_data.rds")