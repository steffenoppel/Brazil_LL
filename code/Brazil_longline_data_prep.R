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


##############################################################
#### LOAD ALL NECESSARY DATASETS
##############################################################






##############################################################
#### READ AND MANIPULATE ATF Fisheries observer data
##############################################################

#### 
dd <- read_excel("OBSERVER_Data_Final_OBSERVERS AND TRIPS_REVISED OCTOBER 2019.xlsx",   ## was "Copy of Copy of Namibia_OBSERVER_Data_Final_July2018.xlsx", 
                 sheet="LONGLINE_DATA")[,1:19]

# remove blank rows
dd <- subset(dd, OBSET_ID != "_")
dd <- subset(dd, !is.na(VESSEL_ID))
dd <- subset(dd, !is.na(Hooks_Observed))
dd <- subset(dd,Hooks_Observed>0)
dd <- subset(dd, !is.na(Birds_Obs_Caught))
dd <- subset(dd, !is.na(Total_Hooks_Set))

#Create seasons
dd$Month <- months(dd$Date_Start_Set)

dd$Season <- "Summer"
dd$Season[dd$Month=="April"]<- "Winter"
dd$Season[dd$Month=="May"]<- "Winter"
dd$Season[dd$Month=="June"]<- "Winter"
dd$Season[dd$Month=="July"]<- "Winter"
dd$Season[dd$Month=="August"]<- "Winter"
dd$Season[dd$Month=="September"]<- "Winter"

#Catch rate on hauls
dd$rate <- dd$Birds_Obs_Caught/dd$Hooks_Observed*1000

#Fix two sets where n observed sets is the same as total hooks set
for(i in 1:nrow(dd)){
  if(dd$Hooks_Observed[i]==dd$Total_Hooks_Set[i]){
    dd$Hooks_Observed[i] <- dd$Hooks_Observed[i]-1
  }
}

dim(dd)
head(dd)


# #First, for fisheries observer data, corrected observed hooks
# 
# dd$Hooks_Observed_Corr_Fac <- dd$Hooks_Observed*rnorm(373,0.627,0.021)
# dd$Hooks_Observed_BL <- dd$Total_Hooks_Set*rnorm(373,0.563,0.019)
# dd$rate <- dd$Birds_Obs_Caught/dd$Hooks_Observed_Corr_Fac*1000


## RETAIN ONLY COLUMNS NECESSARY FOR ANALYSIS AND FIX ERRORS
FishObsDat<-dd %>% dplyr::select(VESSEL_ID,Trip_ID,OBSET_ID,Total_Hooks_Set,Hooks_Observed,BSL,Birds_Obs_Caught,Season) %>%
  mutate(BSL=ifelse(BSL=="no","No",ifelse(BSL=="Single","Yes",BSL))) %>%
  mutate(Trip_ID=paste("FAO",Trip_ID, sep="")) %>%
  mutate(Hooks_Recovered=Hooks_Observed) %>%   ### ATF believes that the 'hooks_observed from FAO observers is actually the Hooks_Recovered
  mutate(Hooks_Observed=NA) %>%
  mutate(Year=2018) %>%
  mutate(Regulation="AFTER")
head(FishObsDat)
dim(FishObsDat)




##############################################################
#### READ AND MANIPULATE ATF Birdlife observers data
##############################################################

# read in data
# first 19 columns are data, others are reference/lookup
dat <- read_excel("Namibia_Demersal_Longline_Data_October 18.xlsx", ### changed from "Copy of Namibia_Demersal_Longline_Data_revised 240818.xlsx"
                  sheet="5_Haul")[,1:15]

corpse <- read_excel("Namibia_Demersal_Longline_Data_October 18.xlsx", ### changed from "Copy of Namibia_Demersal_Longline_Data_revised 240818.xlsx"
                     sheet="8_Bycatch")[,1:14]

gear <- read_excel("Namibia_Demersal_Longline_Data_October 18.xlsx", 
                   sheet="4_Set")[,1:35]

# remove blank rows
dat <- subset(dat, OBSET_ID != "_")
corpse <- subset(corpse, OBSET_ID != "_")

dat$mortality <- 0
# grab the bits of data & link by OBSPERIOD_ID

corpse.shot <- corpse %>%   
  group_by(OBSET_ID) %>%
  summarize(tot_mort = max(Bycatch_ID)) 

for (i in 1:nrow(dat)){
  for(j in 1:nrow(corpse.shot)){
    if(dat$OBSET_ID[i]==corpse.shot$OBSET_ID[j]) dat$mortality[i] <- corpse.shot$tot_mort[j]
  }
}

dat$Month <- months(dat$Date_Start_Haul)

dat$Season <- "Summer"
dat$Season[dat$Month=="April"]<- "Winter"
dat$Season[dat$Month=="May"]<- "Winter"
dat$Season[dat$Month=="June"]<- "Winter"
dat$Season[dat$Month=="July"]<- "Winter"
dat$Season[dat$Month=="August"]<- "Winter"
dat$Season[dat$Month=="September"]<- "Winter"

dat$Year <- 2017
dat$Year[(dat$Date_Start_Haul>(as.POSIXct("2017-06-30",format="%Y-%m-%d")))]<- 2018
dat$Year[(dat$Date_Start_Haul>(as.POSIXct("2018-06-30",format="%Y-%m-%d")))]<- 2019
dat <- dat[dat$Hooks_Observed<=dat$Hooks_Recovered,] #REMOVE one set with more observed hooks than set hooks

dim(dat)
head(dat)

### ADD INFORMATION ON BSL 

gear <- subset(gear, OBSET_ID != "_")
head(gear)
head(dat)

gear<- gear %>% select(OBSET_ID,BSL,Offal_discard)
dat <- dat %>% left_join(gear, by="OBSET_ID")



## RETAIN ONLY COLUMNS NECESSARY FOR ANALYSIS AND FIX ERRORS
ATFDat<-dat %>% dplyr::select(OBSTRIP_ID,OBSET_ID,Hooks_Recovered,Hooks_Observed,BSL,mortality,Season,Year) %>%
  mutate(Trip_ID=OBSTRIP_ID) %>%
  mutate(Total_Hooks_Set=NA) %>%
  rename(VESSEL_ID=OBSTRIP_ID, Birds_Obs_Caught=mortality) %>%
  mutate(Regulation="AFTER")
head(ATFDat)
dim(ATFDat)





##############################################################
#### READ AND MANIPULATE PRE-REGULATION DATA
##############################################################

d.BSL <- read.csv("namibian longline data with no BSL new seasons.csv",header=T)
head(d.BSL)
# d.BSL <- d.BSL[!is.na(d.BSL$TripNo),]
# dd4 <- dd3[,c(1,8,11:13,15)]
# dd4$Reg <- 1
# d.BSL1 <- d.BSL[,c(2,17,21,18,4,3)]
# d.BSL1$Reg <- 0
# names(d.BSL1) <- c("VESSEL_ID","Birds_Obs_Caught","Season","rate","Hooks_Observed_Corr_Fac","Year","Reg")
# dd5 <- rbind(dd4,d.BSL1)
# head(dd5)

preregdat<- d.BSL %>% dplyr::select(TripNo,Vessel,Year,HooksLine,SetNo,TotalMortality,Season) %>%
  mutate(Total_Hooks_Set=NA, BSL="No") %>%
  rename(VESSEL_ID=Vessel, Trip_ID=TripNo, OBSET_ID=SetNo,Birds_Obs_Caught=TotalMortality,Hooks_Observed=HooksLine) %>%
  mutate(Trip_ID=paste("PreREG",Trip_ID, sep="")) %>%
  mutate(OBSET_ID=paste(VESSEL_ID,OBSET_ID, sep="")) %>%
  mutate(Regulation="BEFORE")
  








##############################################################
#### JOIN ALL PRE- AND POSTREGULATION DATA
##############################################################
head(ATFDat)
head(FishObsDat)
head(preregdat)

#Join files into file called dd2
# names(dd)
# dd2 <- dd[c(3,7,11,12,15:25)]
# dat2 <- data.frame("VESSEL_ID"=dat$OBSTRIP_ID,"Date_Start_Set"=dat$Date_Start_Haul, 
#                    "Latitude_Start_Set"=dat$Latitude_Start_Haul,"Longitude_Start_Set"=dat$Longitude_Start_Haul,
#                    "Total_Hooks_Set"=dat$Hooks_Recovered,"Hooks_Observed"=dat$Hooks_Observed,
#                    "BSL"=dat$BSL,"Birds_Obs_Caught"=dat$mortality,"Offal_discard"=dat$Offal_discard,
#                    "Month"=dat$Month,"Season"=dat$Season,"Hooks_Observed_Corr_Fac"=dat$Hooks_Observed,
#                    "Hooks_Observed_BL"=dat$Hooks_Observed,"rate"=dat$rate,"Year"=dat$Year)
# dd3 <- rbind(dd2,dat2)


ALL_DATA<- bind_rows(FishObsDat,ATFDat,preregdat)
head(ALL_DATA)
dim(ALL_DATA)

#######################################################################################





##############################################################
#### LOOK AT THE DATA
##############################################################
unique(ALL_DATA$Year)
ALL_DATA %>% mutate(rate=Birds_Obs_Caught/Hooks_Observed) %>%
  group_by(Regulation,BSL) %>%
  summarise(Bycatch=mean(rate)*1000) 


fwrite(ALL_DATA,"ATF_Namibia_Longline_Bycatch_data2009_2019.csv")



##############################################################
#### FAO OBSERVERS DID NOT RECORD THE NUMBER OF ACTUALLY OBSERVED HOOKS
##############################################################

## Calculate the proportion of observed hooks from ATF
## assume that proportion for FAO was similar


## RATIO for ATF
ALL_DATA %>% dplyr::filter(!is.na(Hooks_Recovered)) %>%
  dplyr::filter(!is.na(Hooks_Observed)) %>%
  mutate(ratio=Hooks_Observed/Hooks_Recovered) %>%
  group_by(Year, Regulation) %>%
  summarise(mean=mean(ratio,na.rm=T), sd=sd(ratio,na.rm=T))


## THE FAO Observer ratio is much higher (and cannot be trusted)
ALL_DATA %>% dplyr::filter(!is.na(Total_Hooks_Set)) %>%
  mutate(ratio=Hooks_Recovered/Total_Hooks_Set) %>%
  group_by(Year, Regulation) %>%
  summarise(mean=mean(ratio,na.rm=T), sd=sd(ratio,na.rm=T))



