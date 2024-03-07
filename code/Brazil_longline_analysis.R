### ##########################################
###
### BRAZIL longline bycatch analysis
###
### ##########################################

### written by Steffen Oppel 12 February 2024
## data cleaning and preparation in "Brazil_longline_data_prep.R"
## model based on Gardner et al. 2008 and Field et al. 2019
## based on scripts from Namibia, Da Rocha et al. 2021: https://www.sciencedirect.com/science/article/abs/pii/S0006320720309733

## updated on 16 February 2024 to exclude data from before 2009
## updated 24 February to include more extensive data exploration - RF cannot find a signal, so added univariate exploration


##############################################################
#### load ALL NECESSARY libraries
##############################################################

library(tidyverse)
library(lubridate)
library(data.table)
library(runjags)
library(MCMCvis)
library(randomForest)
library(dplyr)
library(ggplot2)
filter<-dplyr::filter
select<-dplyr::select


##############################################################
#### LOAD PRE_ARRANGED DATA
##############################################################

setwd("C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL")
data<-readRDS("data/Brazil_formatted_bycatch_data2009_2018.rds") %>% mutate(BYCATCH_BIN=ifelse(BYCATCH==0,0,1))
alldata<-readRDS("data/Brazil_formatted_bycatch_data.rds") %>% mutate(BYCATCH_BIN=ifelse(BYCATCH==0,0,1))
head(data)
summary(data)



##############################################################
#### INSPECT DATA DISTRIBUTION
##############################################################
table(data$BYCATCH)

data %>% filter(BYCATCH>3)
alldata %>% filter(BYCATCH>3)

### ABUNDANCE OF BYCATCH REDUCED FROM 0.19 to 0.07 / 1000 hooks

data %>% group_by(Toriline ) %>%
  summarise(mean=mean(BPUE))

### OCCURRENCE OF BYCATCH REDUCED FROM 9% to 6% (very little)

data %>% mutate(BYCATCH_BIN=ifelse(BYCATCH==0,0,1)) %>%
  group_by(Toriline) %>%
  summarise(mean=mean(BYCATCH_BIN))

unique(data$Month)
unique(data$Year)


### CONFIRM THAT THERE IS NO BENEFIT TO KEEP ANY DATA FROM BEFORE 2009
### ABUNDANCE OF BYCATCH REDUCED FROM 0.14 to 0.09 - makes no difference to use old data

alldata %>% filter(!(Year<2009 & Toriline=="Yes")) %>%
  group_by(Toriline) %>%
  summarise(mean=mean(BPUE))

### OCCURRENCE OF BYCATCH REDUCED FROM 7% to 6% (even worse when just using reduced data)

alldata %>% filter(!(Year<2009 & Toriline=="Yes")) %>%
  mutate(BYCATCH_BIN=ifelse(BYCATCH==0,0,1)) %>%
  group_by(Toriline) %>%
  summarise(mean=mean(BYCATCH_BIN))



#### EXPLORE DISTRIBUTION OF VARIABLES WITH AND WITHOUT BYCATCH

data %>% mutate(BYCATCH_BIN=ifelse(BYCATCH==0,0,1)) %>%
  gather(key="Variable",value="Value", -Set,-Date,-BYCATCH_BIN,-BYCATCH,-BPUE,-Daylight_set,-N_hooks,-Night_set) %>%
  ggplot(aes(x=as.factor(BYCATCH_BIN), y=Value)) +
  facet_wrap(~Variable, scales="free_y", ncol=3) +
  geom_boxplot()


#### EXPLORE BYCATCH DATA DISTRIBUTION

data %>% #filter(BYCATCH>0) %>%
  ggplot(aes(x=BPUE, fill=Toriline)) +
  geom_histogram()


##############################################################
#### PRELIMINARY EXPLORATORY ANALYSIS TO DETERMINE WHICH VARIABLES ARE IMPORTANT
##############################################################
### very poor performance of these models
### will not be able to explain a lot of bycatch variation!!
### including BYCATCH_BIN in abundance model takes explained variance to >50% but even that is poor given that BYCATCH_BIN will explain all the 0s
RF<-randomForest(BPUE~N_hooks+Toriline+Year+Month+Latitude+Longitude+Season+Nightlight_set+Moon.il, data=data, mtry=3,ntree=1500, importance=T, replace=F)
varImpPlot(RF)
RF

RFbin<-randomForest(as.factor(BYCATCH_BIN)~N_hooks+Toriline+Year+Month+Latitude+Longitude+Season+Nightlight_set+Moon.il, data=data, mtry=3,ntree=1500, importance=T, replace=F)
varImpPlot(RFbin)
RFbin

RFabund<-randomForest(BPUE~N_hooks+Toriline+Year+Month+Latitude+Longitude+Season+Nightlight_set+Moon.il, data=data[data$BYCATCH_BIN==1,], mtry=3,ntree=1500, importance=T, replace=F)
varImpPlot(RF)
RFabund


par(mfrow=c(3,2))
partialPlot(RF,as.data.frame(data),"Longitude")
partialPlot(RF,as.data.frame(data),"Latitude")
partialPlot(RF,as.data.frame(data),"Season")
partialPlot(RF,as.data.frame(data),"Nightlight_set")
partialPlot(RF,as.data.frame(data),"Moon.il")
partialPlot(RF,as.data.frame(data),"Year")

#####  PLOT VARIABLE IMPORTANCE FOR BOTH MODELS ###########

VARbinary<-importance(RFbin, type=1)
VARnum<-importance(RF, type=1)


IMP<-bind_rows(VARbinary[,1],VARnum[,1]) %>%
  mutate(MAX = do.call(pmax, (.))) %>%
  mutate(model=c("occurrence","BPUE")) %>%
  gather(key=variable, value=MSE,-model,-MAX) %>%
  mutate(IMP=(MSE/MAX)*100) %>%
  arrange(model,desc(IMP))
    
    
ggplot(IMP, aes(x=variable, y=IMP)) +
  geom_bar(stat='identity', fill='lightblue') +
  coord_flip()+
  facet_wrap(~model, ncol=2) +
  ylab("Variable importance (%)") +
  xlab("Explanatory variable") +
  scale_y_continuous(limits=c(-20,115), breaks=seq(0,100,20), labels=seq(0,100,20))+
  scale_x_discrete(name="",limit = c("Toriline","Moon.il","N_hooks","Month","Latitude","Season","Nightlight_set","Longitude","Year"))  +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.x=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text=element_text(size=18, color="black"),
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


# ggsave("output/Brazil_LL_RF_variable_importance.jpg", width=9, height=9)




##############################################################
#### DATASET IS POORLY BALANCED - RESAMPLE
##############################################################
### retain all Toriline samples
### randomly sample same number of no Toriline samples
### repeat a 1000 times
### calculate bootstrapped confidence interval

table(data$Toriline)

BSL<- data %>% filter(Toriline=="Yes") 
none<- data %>% filter(Toriline=="No") 



## bootstrapping the BPUE samples (10000 random draws)
BSL.samples <- matrix(sample(BSL$BPUE, size = 10000 * nrow(BSL), replace = TRUE),10000, nrow(BSL))
BSL.statistics <- apply(BSL.samples, 1, mean)
BPUE_BSL<-data.frame(treatment="BSL",mean=mean(BSL.statistics),
                         lcl=quantile(BSL.statistics,0.025),ucl=quantile(BSL.statistics,0.975))

none.samples <- matrix(sample(none$BPUE, size = 10000 * nrow(BSL), replace = TRUE),10000, nrow(BSL))
none.statistics <- apply(none.samples, 1, mean)
BPUE_none<-data.frame(treatment="none",mean=mean(none.statistics),
                     lcl=quantile(none.statistics,0.025),ucl=quantile(none.statistics,0.975))


### calculating the BPUE differences

BSL_DIFF_BPUE<-data.frame(treatment="BSL_difference",mean=mean(BSL.statistics-none.statistics),
                           lcl=quantile(BSL.statistics-none.statistics,0.025),ucl=quantile(BSL.statistics-none.statistics,0.975))

### quantifying differences in per cent
BSL_DIFF_BPUE[,5:7]<-(BSL_DIFF_BPUE[,2:4]/BPUE_none[,2])*100


BSL_DIFF_BPUE


### COMBINE NUMBERS FOR OUTPUT FILE

DIFF_SUMMARY<-bind_rows(BPUE_BSL,BPUE_none,BSL_DIFF_BPUE) %>%
  rename(`prop.mean.diff(%)`=mean.1,`prop.lcl.diff(%)`=lcl.1,`prop.ucl.diff(%)`=ucl.1) %>%
  rename(`abs.mean.diff(N)`=mean,`abs.lcl.diff(N)`=lcl,`abs.ucl.diff(N)`=ucl)


fwrite(DIFF_SUMMARY,"output/Brazil_LL_Toriline_bootstrap_summary.csv")



### PLOT predicted OUTPUT FOR FISH AND SEABIRDS ###

bind_rows(BPUE_BSL,BPUE_none) %>%
  ggplot(aes(y=mean, x=treatment, col=treatment)) + geom_point(size=3)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.05)+
  scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1)) +
  xlab("") +
  ylab("Bycaught seabirds / 1000 hooks") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.position="none",
        strip.background = element_blank(),
        strip.placement = "outside", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("Brazil_LL_Toriline_bycatch_difference.jpg", width=8, height=7)





##############################################################
#### SPECIFY AND RUN JAGS MODELS 
##############################################################
### need to include lat, long, season, moon, night, hooks, toriline as fixed effects
### include random effects for month and year
### suggestions for how to include effort offset: https://stats.stackexchange.com/questions/446929/how-to-correctly-include-offset-in-bayesian-zero-inflated-poisson-model-in-winbu
### not implemented yet - just reduced variables

sink ("code/ATF_Brazil_LongLine_Bycatch_v2.jags")
cat(" 
model{
  
  # PRIORS FOR REGRESSION PARAMETERS
  
  intercept.occu ~ dnorm(0, 0.01)  ## intercept for occurrence of bycatch
  intercept.abund ~ dnorm(0, 0.01)  ## intercept for quantity of bycatch
  rho ~ dunif(0,50)  ## overdispersion parameter for negative binomial distribution

  tori.occu ~ dnorm(0, 0.01)
  tori.abund ~ dnorm(0, 0.01)
   
  breed.occu ~ dnorm(0, 0.01)
  breed.abund ~ dnorm(0, 0.01)
  
  effort.occu ~ dnorm(0, 0.01)
  effort.abund ~ dnorm(0, 0.01)
  

  # RANDOM YEAR EFFECTS FOR OCCURRENCE AND ABUNDANCE
  for(t in 1:nyears){
    occ.year[t]~dnorm(0,tau.occ.year)    ## trip-specific random effect for occurrence
    abund.year[t]~dnorm(0,tau.ab.year)    ## trip-specific random effect for abundance
    }
  tau.occ.year<-1/(sigma.occ.year*sigma.occ.year)
  sigma.occ.year~dunif(0,10)
  tau.ab.year<-1/(sigma.ab.year*sigma.ab.year)
  sigma.ab.year~dunif(0,10)
  

  # RANDOM MONTH EFFECTS FOR OCCURRENCE AND ABUNDANCE
  for(t in 1:12){
    occ.month[t]~dnorm(0,tau.occ.month)    ## month-specific random effect for occurrence
    abund.month[t]~dnorm(0,tau.ab.month)    ## month-specific random effect for abundance
  }
  tau.occ.month<-1/(sigma.occ.month*sigma.occ.month)
  sigma.occ.month~dunif(0,10)
  tau.ab.month<-1/(sigma.ab.month*sigma.ab.month)
  sigma.ab.month~dunif(0,10)
  


  #### LIKELIHOOD LOOP OVER  every set of longlines

  for(i in 1:N){
    
    # define the logistic regression model, where psi is the probability of bycatch occurring at all
    cloglog(psi[i]) <- intercept.occu +
                        effort.occu*N_hooks[i] +
                        tori.occu*tori[i] +
                        breed.occu*season[i] +
                        
                        occ.year[year[i]] +
                        occ.month[month[i]]
    z[i]~dbern(psi[i])
    
    # define the negative binomial regression model for abundance and multiply with bycatch probability
    mortality[i] ~ dnegbin(phi[i],rho)
    phi[i] <- rho/(rho+(z[i])*lambda[i]) - 1e-10*(1-z[i])
    log(lambda[i])<- intercept.abund +
                      effort.abund*N_hooks[i] +
                      tori.abund*tori[i] +
                      breed.abund*season[i]  +
                      abund.year[year[i]] +
                      abund.month[month[i]]
    
  } ## end loop over each observation
  
  
  # CONVERT TO ESTIMATES PER 1000 HOOKS BEFORE AND AFTER REGULATION
  #prereg <- (((1-exp(-exp(intercept.occu+log(1000)))))*(exp(intercept.abund)*1000))
  #postreg <- (((1-exp(-exp(intercept.occu+treat.occu+log(1000)))))*(exp(intercept.abund + treat.abund)*1000))

  
} ## end model
    ")
sink()



##############################################################
#### STANDARDIZE VARIABLES WITH LARGE NUMERICAL VALUES 
##############################################################

nhooks_scaled<-scale(data$N_hooks)
lat_scaled<-scale(data$Latitude)
long_scaled<-scale(data$Longitude)
lat_scaled<-scale(data$Latitude)
lat_scaled<-scale(data$Latitude)



##############################################################
#### PREPARE DATA FOR JAGS
##############################################################


jags.data <- list(mortality=data$BYCATCH,
                  N=dim(data)[1],
                  nyears=length(unique(data$Year)),
                  nmonths=length(unique(data$Month)),
                  year=as.numeric(as.factor(data$Year)),
                  month=as.numeric(data$Month),
                  tori=ifelse(data$Toriline=="No",0,1),
                  season=ifelse(data$Season=="Breading",1,0),
                  lat=as.numeric(lat_scaled),
                  long=as.numeric(long_scaled),
                  N_hooks=as.numeric(nhooks_scaled),
                  moon=data$Moon.il,
                  night=data$Nightlight_set
                  )


##############################################################
#### SET UP INITS AND PARAMETERS
##############################################################

inits = function(){list(tori.occu=rnorm(1,0,0.1),
                        tori.abund=rnorm(1,0,0.1),
                        night.occu=rnorm(1,0,0.1),
                        night.abund=rnorm(1,0,0.1)
                        )}

params <- c("tori.abund","moon.abund","night.abund","lat.abund","long.abund","breed.abund","effort.abund",
            "tori.occu","moon.occu","night.occu","lat.occu","long.occu","breed.occu","effort.occu")

n.chains = 4 
n.burnin = 500
n.iter = 1000
n.thin = 5
n.adapt = 500
n.sample = 200



##############################################################
#### RUN JAGS MODEL [takes 35 minutes]
##############################################################


# BRA_LL <- jagsUI(data=jags.data,
#                   model = "C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL/code/ATF_Brazil_LongLine_Bycatch_v1.jags",
#                   inits=inits,
#                   parameters.to.save =params,
#                   n.chains=n.chains,
#                   n.thin = n.thin,
#                   n.iter = n.iter,
#                   n.burnin = n.burnin, parallel=T, n.cores=4)



BRA_LL <- run.jags(data=jags.data, inits=inits, monitor=params,
                             model="C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL/code/ATF_Brazil_LongLine_Bycatch_v2.jags",
                             n.chains = n.chains, thin = n.thin, burnin = n.burnin, adapt = n.adapt,sample = n.sample, 
                             method = "rjparallel") 






#### MODEL ASSESSMENT ####
MCMCplot(BRA_LL$mcmc, params=c("tori.abund","moon.abund","night.abund","lat.abund","long.abund","breed.abund","effort.abund",
                                             "tori.occu","moon.occu","night.occu","lat.occu","long.occu","breed.occu","effort.occu"))
# ggsave("C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL/output/Fig_S1_parameter_estimates.jpg", height=11, width=8)
# ggsave("C:/STEFFEN/OneDrive - Vogelwarte/General - Little owls/MANUSCRIPTS/LittleOwlSurvival/Fig_S1_parameter_estimates.jpg", height=11, width=8)

MCMCtrace(BRA_LL$mcmc)
MCMCsummary(BRA_LL$mcmc)
MCMCdiag(BRA_LL$mcmc,
         round = 3,
         file_name = 'Brazil_LL_v1',
         dir = 'C:/STEFFEN/OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS/STEFFEN/RSPB/Marine/Bycatch/Brazil_LL/output',
         mkdir = 'Brazil_LL_v1',
         add_field = '1.0',
         add_field_names = 'Data version',
         save_obj = TRUE,
         obj_name = 'Brazil_LL_v1',
         add_obj = list(jags.data, sessionInfo()),
         add_obj_names = c('data-09Feb2024', 'session-info-09Feb2024'))








###############################################################################
####   CALCULATE CHANGE FROM POSTERIOR SAMPLES   #############################
###############################################################################


change<-(BRA_LL$sims.list$prereg-BRA_LL$sims.list$postreg)/BRA_LL$sims.list$prereg
quantile(change,0.5)
quantile(change,0.025)
quantile(change,0.975)





###############################################################################
####   CREATE OUTPUT TABLE FOR MANUSCRIPT   #############################
###############################################################################


LLsum<-data %>% mutate(Hooks_Observed=ifelse(is.na(Hooks_Observed),Hooks_Recovered*0.56,Hooks_Observed)) %>%
  mutate(rate=(Birds_Obs_Caught/Hooks_Observed)*1000) %>%
  group_by(Regulation) %>%
  summarise(byc_rate=mean(rate)) %>%
  bind_cols(as.data.frame(NamLLmodel$summary[4:3,c(1,3,7)])) %>%
  rename(lcl=`2.5%`,ucl=`97.5%`) %>%
  arrange(desc(Regulation))


#### CALCULATE THE CHANGE IN INTERACTION RATE ####
percchange<-function(x){((x[1]-x[2])/x[1])*100}
LLsum[3,2]<-apply(as.matrix(LLsum[,2]),2,percchange)
LLsum[3,3:5]<-c(mean(change,na.rm=T)*100,quantile(change,0.025)*100,quantile(change,0.975)*100)
LLsum[3,1]<-"CHANGE(%)"
LLsum

fwrite(LLsum,"Namibia_longline_bycatch_REG_comparison.csv")




###############################################################################
####   EVALUATE MODEL FIT WITH BAYESIAN P VALUE   #############################
###############################################################################
## does not fit - but not sure whether fit statistic is indeed calculated correctly
## abandoned from v5 onwards

# plot(NamLLmodel$sims.list$fit, NamLLmodel$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
# abline(0, 1, lwd = 2, col = "black")
# mean(NamLLmodel$sims.list$fit.new > NamLLmodel$sims.list$fit)
# mean(NamLLmodel$mean$fit) / mean(NamLLmodel$mean$fit.new)









##############################################################
#### EXTRAPOLATION OF FLEET-WIDE MORTALITY
##############################################################
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\Namibia\\Data")

## read in fleet-wide effort data
LLsum<-fread("Namibia_longline_bycatch_REG_comparison.csv")

LLeff <- read_excel("LL data combined 2016-2018.xlsx", 
                    sheet="Sheet1")
head(LLeff)




## manipulate and summarise fleet-wide seabird mortality based on fatal interactions

LLsummary<- LLeff %>%
  rename(effort=`Number of HOOKS_SET`) %>%
  mutate(Regulation=if_else(Year>2015,"AFTER","BEFORE")) %>%
  group_by(Regulation,Year) %>%
  summarise(tot_effort=sum(effort,na.rm=T)) %>%
  mutate(bycatch=tot_effort*LLsum$mean[match(Regulation,LLsum$Regulation)]/1000,
         bycatch.lcl=tot_effort*LLsum$lcl[match(Regulation,LLsum$Regulation)]/1000,
         bycatch.ucl=tot_effort*LLsum$ucl[match(Regulation,LLsum$Regulation)]/1000) %>%
  mutate(Fishery="Longline")




fwrite(LLsummary,"Namibia_fleetwide_longline_seabird_mortality.csv")






##############################################################
#### ASSESSMENT OF SAMPLING REPRESENTATIVITY
##############################################################
LLsummary<-LLsummary %>% select(Regulation,Year,tot_effort) %>%
  bind_rows(data.frame(Regulation="BEFORE",Year =as.numeric(c(2009,2010,2011,2012)), tot_effort=c(47481331,33071604,26131173,31225687)))


data %>% 
  mutate(Total_Hooks_Set=ifelse(is.na(Total_Hooks_Set),Hooks_Recovered/0.913,Total_Hooks_Set)) %>%
  mutate(Total_Hooks_Set=ifelse(is.na(Total_Hooks_Set),Hooks_Observed/(0.913*0.58),Total_Hooks_Set)) %>%
  group_by(Year, Regulation) %>%
  summarise(obs=sum(Total_Hooks_Set)) %>%
  inner_join(LLsummary, by=c("Year","Regulation")) %>%
  ungroup() %>%
  group_by(Regulation) %>%
  summarise(obs=sum(obs),tot_effort=sum(tot_effort)) %>%
  mutate(obs_ratio=obs/tot_effort) %>%
  select(tot_effort,obs,obs_ratio)



##############################################################
#### TEST BASIC MODEL TO ESTIMATE NUMBER OF HOOKS OBSERVED 
##############################################################


sink ("ATF_obsHooks.jags")
cat(" 
    model{
    
    ####  ESTIMATE OBSERVED HOOKS (NO DATA FOR FAO OBSERVERS)
    ## obsHook are provided data for 2 data sets, but NA for the third data set - need to interpolate
    
    ### PRIORS FOR REGRESSION PARAMETERS

    logitprop.obs.mn ~ dnorm(0,0.001)

    # RANDOM TRIP EFFECTS FOR OCCURRENCE AND ABUNDANCE
    for(t in 1:ntrips){
      obs.trip[t]~dnorm(0,tau.obs.trip)    ## trip-specific random effect for occurrence
    }
    tau.obs.trip<-1/(sigma.obs.trip*sigma.obs.trip)
    sigma.obs.trip~dunif(0,10)
    
    for(i in 1:Nobs){

    obsHook[i] ~ dbin(prop.obs[i],retrHook[i])
    logit(prop.obs[i])<-logitprop.obs.mn + obs.trip[trip[i]]						### random effect for observer on each trip
    #logit(prop.obs[i]) ~ dnorm(logitprop.obs.mn, logitprop.obs.tau)   ### this results in a syntax error

    } ## end loop over each observation that has  
    
    ratio<-mean(prop.obs[])


    
    } ## end model
    ")
sink()





##############################################################
#### PREPARE DATA FOR JAGS
##############################################################


jags.data <- list(retrHook=round(data$Hooks_Recovered),   ## maybe need to scale as in the 1000s? model runs forever if we do, but only 5 min if we leave original N
                  obsHook=round(data$Hooks_Observed),
                  ntrips=length(unique(data$Trip_ID)),
                  trip=as.numeric(as.factor(data$Trip_ID)),
                  Nobs=dim(data[!is.na(data$Hooks_Recovered),])[1])

inits = function(){list(logitprop.obs.mn=rnorm(1,0,0.1),
                        logitprop.obs.tau=rnorm(1,0,0.1))}

params <- c("ratio")




##############################################################
#### RUN JAGS MODEL [takes>25 hours]
##############################################################
n.chains = 3 
n.burnin = 10000
n.iter = 50000
n.thin = 1


test.obs <- jagsUI(data=jags.data,
                  model = "C:/STEFFEN/RSPB/Marine/Bycatch/Namibia/ATF_obsHooks.jags",
                  inits=inits,
                  parameters.to.save =params,
                  n.chains=n.chains,
                  n.thin = n.thin,
                  n.iter = n.iter,
                  n.burnin = n.burnin, parallel=T, n.cores=8)
