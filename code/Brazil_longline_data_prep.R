### ##########################################
###
### Brazilian longline bycatch analysis
### DATA PREPARATION
### based on scripts from Namibia, Da Rocha et al. 2021: https://www.sciencedirect.com/science/article/abs/pii/S0006320720309733
###
### ##########################################

### adapted up by Steffen Oppel 31 July 2019
### data provided by Yann Rouxel and Dimas Gianuca in single xlsx spreadsheet

### updated on 16 Feb 2024 after Dimas Gianuca explained data structure
### BSL prior to 2011 are experimental and should not be used



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
library(tmap)
library(sf)


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
  mutate(Date=as.Date(Date, origin = "1899-12-30")) %>%
  mutate(BPUE=BYCATCH/(N_hooks/1000)) %>%
  filter(Latitude>-50)  ## filter out one suspicious location south-east of the Falklands

summary(dd)
dim(dd)

##############################################################
#### REMOVE THE EARLY DESIGN BSL
##############################################################

dd_red<- dd %>% filter(Year>2008)  ## remove data from before 2009
dim(dd_red)
table(dd_red$Toriline)
table(dd$Toriline)


##############################################################
#### CHECK THAT COORDINATES MAKE SENSE
##############################################################

#map(database = "world") 
#points(x = dd$Latitude, y = dd$Longitude, col="firebrick", pch=21)

dd_sf<-dd %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs=4326)  

tmap_mode("view")
  tm_shape(dd_sf)  +
  tm_symbols(col = 'BPUE', size = 0.1)


##############################################################
#### SAVE PROCESSED DATA
##############################################################
saveRDS(dd,"data/Brazil_formatted_bycatch_data.rds")
saveRDS(dd_red,"data/Brazil_formatted_bycatch_data2009_2018.rds")
