##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS
#
##########################################################################
# based on Lieury et al. 2015
# modified by Steffen oppel 4 April 2018
# goal is to assess population viability with and without captive released juveniles on Balkans

# curtailed data from 2006 - 2017 because we have territory monitoring data for that period

# pasted together data prep and model from various files on 16 May 2018
# update 17 May 2018: trying to troubleshoot dbinom error
# updated 22 May 2018: breed.prop cannot be a stochastic node, because 'slicer' gets stuck - tried with multiple distributions, priors and initials, but never works
# included Adam's suggestions for dbinom and output treatment (i.e. abandoned the idea of solving for CAPT.ADD)
# updated 23 May 2018: changed priors for fecundity as output was essentially 0-2

#### MAJOR REVISION after Aqaba meeting in April 2019
# included data from 2019
# changed survival probability of released birds to be similar to second-year birds (not juveniles, due to delayed release)
# 

library(readxl)
library(jagsUI)
library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
library(markdown)
library(rmarkdown)
library(knitr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE POPULATION MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
head(breed)

trendinput<- occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  group_by(year) %>%
  summarise(N=sum(occupancy))


breedinput<- breed %>% filter(Year>2005) %>%
  filter(!is.na(breed_success)) %>%
  mutate(count=1) %>%
  mutate(fledglings=ifelse(is.na(fledglings),0,fledglings)) %>%
  group_by(Year) %>%
  summarise(R=sum(count), J=sum(fledglings))








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE JUVENILE TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE)
load("RODBC_EGVU_telemetry_input.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A TIME SERIES DATA FRAME ###
mindate<-min(locs$Date)
maxdate<-max(locs$Date)
timeseries<-data.frame(date=seq(mindate, maxdate, "1 month"), month=0, year=0)
timeseries$month<-as.numeric(format(timeseries$date,"%m"))
timeseries$year<-as.numeric(format(timeseries$date,"%Y"))
dim(timeseries)

### CREATE A BLANK CAPTURE HISTORY for 5 years (60 months) for juvenile birds ###

CH.telemetry<-birds %>%
  filter(Age=="juv") %>%
  select(Name,Tag_year, origin) %>% 
  filter(!Name=="Zighmund")             ### remove single bird that was never free flying
CH.telemetry[,4:63]<-0									### NEEDS MANUAL ADJUSTMENT IF REPEATED FOR LONGER TIME SERIES



### FILL CAPTURE HISTORY WITH LIVE DEAD INFORMATION ###
### THIS MAY REQUIRE A FIX IN THE FUTURE IF BIRDS THAT ARE STILL ALIVE HAVE NOT BEEN TAGGED > 4 YEARS AGO
### WOULD NEED FILLING WITH NA

for(n in CH.telemetry$Name){
  xl<-locs[locs$Bird_ID==n,]
  mindate<-as.Date(birds$Fledge_date[birds$Name==n])
  mindate<-dplyr::if_else(is.na(mindate),as.Date(min(xl$Date)),mindate)
  if(n %in% c("Elodie", "Odiseas")){mindate<-mindate-15}
  if(is.na(birds$Stop_date[birds$Name==n])){maxdate<-as.Date(max(xl$Date))}else{			### for birds that are still alive
    maxdate<-as.Date(birds$Stop_date[birds$Name==n])}
  birdseries<-data.frame(date=seq(mindate, maxdate, "1 month"), live=1)
  stopcol<-length(birdseries$live)+3
  stopcol<-ifelse(stopcol>63,63,stopcol)
  CH.telemetry[CH.telemetry$Name==n,4:stopcol]<-birdseries$live
}

CH.telemetry<-CH.telemetry[order(CH.telemetry$Tag_year),]

y.telemetry<-as.matrix(CH.telemetry[,4:63])


##### CODE FOR FIRST AND and second years #######
## MIGRATION ENDS IN OCTOBER, so only 3 occ (Aug, Sept, Oct) in high risk category

x.telemetry<-y.telemetry
for (l in 1: nrow(x.telemetry)){
  n1 <- min(which(x.telemetry[l,]==1)) 		## find the first marking occasion
  if(n1>1){x.telemetry[l,1:(n1-1)]<-NA}		## set everything before the first occasion to NA
  x.telemetry[l,(n1):dim(x.telemetry)[2]]<-1			## Set everything after the first occasion to age = 1
  nmax <- n1+3		 		## ensure that age=2 is within study period 						$$$ CHANGE HERE IF DIFF TIME INTERVAL REQUIRED $$$
  x.telemetry[l,4:24]<-2				## since all birds start at the same time, age switches to 2 after 6 months	$$$ CHANGE HERE IF DIFF TIME INTERVAL REQUIRED $$$
  x.telemetry[l,25:60]<-3				## since all birds start at the same time, age switches to 2 after 6 months	$$$ CHANGE HERE IF DIFF TIME INTERVAL REQUIRED $$$
}
x.telemetry<-as.matrix(x.telemetry)
str(x.telemetry)



#### Function to create a matrix with information about known latent state z
state.telemetry <- y.telemetry
known.telemetry<-y.telemetry[,1]															#### WE ASSUME THAT FATE OF ALL BIRDS IS KNOWN
ch.telemetry<-state.telemetry
for (i in 1:dim(ch.telemetry)[1]){
  n1 <- min(which(ch.telemetry[i,]==1))
  n2 <- max(which(ch.telemetry[i,]==1))
  state.telemetry[i,n1:n2] <- 1
  state.telemetry[i,1:n1] <- NA
  if(n2<dim(ch.telemetry)[2]){
    if(known.telemetry[i]==1){state.telemetry[i,(n2+1):dim(state.telemetry)[2]]<-0}else{state.telemetry[i,(n2+1):dim(state.telemetry)[2]]<-NA}
  }
  state.telemetry[state.telemetry==0] <- NA
}
z.telemetry<-state.telemetry
z.telemetry<-as.matrix(z.telemetry)


#### create vector of first marking (this is 1 by default)
get.first.telemetry<-function(x)min(which(x!=0))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)



#### Bundle data FOR JAGS MODEL RUN and save workspace
INPUT.telemetry <- list(y.telemetry = y.telemetry, f.telemetry = f.telemetry, nind.telemetry = dim(CH.telemetry)[1], n.occasions.telemetry = dim(y.telemetry)[2], z.telemetry = z.telemetry, x.telemetry = x.telemetry)
#rm(list=setdiff(ls(), "INPUT.telemetry"))



#### Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,],na.rm=T)==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}


# Function to create a matrix of initial values for latent state z
ch.init <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  return(ch)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE TERRITORY MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_surveys.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)

load("RODBC_EGVU_surveys.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM MACEDONIA AND MERGE WITH BG DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#EGVUm<-read_excel("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Raw_Data\\archived_databases\\ReturnNeophron\\C4_C5\\Macedonia_mgmt.xlsx", sheet = "Survival")
EGVUm<-read_excel("C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\archived_databases\\ReturnNeophron\\C4_C5\\Macedonia_mgmt.xlsx", sheet = "Survival")

head(EGVUm)
EGVUm$visit_duration<-as.numeric(difftime(EGVUm$time_end,EGVUm$time_start, units = "hours"))
EGVUm$visit_duration[is.na(EGVUm$visit_duration)]<-mean(EGVUm$visit_duration,na.rm=T)

### filter the territories that have very few visits
EGVUm$Year<-as.numeric(format(EGVUm$date, format="%Y"))
EGVUm$terryear<-paste(EGVUm$Year,EGVUm$territory_name, sep="_")
mgmtMAC<-aggregate(count~Country+territory_name+Year+terryear, data=EGVUm, FUN=sum)
shit<-subset(mgmtMAC, count<2)

### remove the worthless surveys from the main table
EGVUm<-EGVUm[!(EGVUm$terryear %in% shit$terryear),] %>%
  mutate(STAGE=NA, STATUS=NA,failure_reason=NA,n_chicks=0,n_other=0,nest_ID=NA) %>%  ## insert the columns needed to join data
  select(names(EGVU))

EGVU<-rbind(EGVU,EGVUm)
dim(EGVU)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMAT DATA AND FIX EFFORT TIMES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(EGVU)

### recalculate effort for visits that have observed >0 birds with 0 effort
updDUR<-EGVU %>% filter(visit_duration==0) %>%
  mutate(dur=as.numeric(difftime(time_end,time_start,units='hours'))) %>%
  mutate(duration=ifelse(dur==0,0.1,dur)) %>%
  select(survey_ID,territory_name,date,time_start,time_end,duration)
EGVU$visit_duration[match(updDUR$survey_ID,EGVU$survey_ID)]<-updDUR$duration

### FORMAT AND SELECT RELEVANT COLUMNS
EGVU<- EGVU %>%
  filter(!is.na(visit_duration)) %>%
  filter(!is.na(n_adults)) %>%
  filter(visit_duration!=0) %>% 
  select(territory_name,date,visit_duration,n_adults) %>%
  mutate(Year=year(date),JDAY=yday(date),WEEK=week(date), MONTH=month(date)) %>%
  mutate(PRIMOCC=paste(Year,MONTH, sep="_")) %>%     ## CREATE UNIQUE PRIMARY OCCASION (MONTH WITHIN YEARS)
  filter(MONTH %in% c(4:8)) %>%
  filter(WEEK %in% c(13:35))

dim(EGVU)
head(EGVU)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE UNOCCUPIED TERRITORIES ACROSS ALL YEARS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
maxocc<-aggregate(n_adults~territory_name, EGVU, FUN=max)
remove<-maxocc[maxocc$n_adults==0,]
names(remove)[2]<-"REMOVE"
EGVU<-merge(EGVU,remove, by=c("territory_name"), all.x=T)
dim(EGVU)
EGVU<-EGVU[is.na(EGVU$REMOVE),]
dim(EGVU)
rm(maxocc, remove)
EGVU$REMOVE<-NULL



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A VISIT NUMBER FOR EACH TERRITORY, YEAR AND WEEK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EGVU$visit<-1

for (t in unique(EGVU$territory_name)){
  tx<-subset(EGVU, territory_name==t)
  for (y in unique(tx$PRIMOCC)){
    ty<-subset(tx, PRIMOCC==y)
    
    ty<-ty[order(ty$JDAY, decreasing=F),]
    
    EGVU$visit[EGVU$territory_name==t & EGVU$PRIMOCC==y] <- as.numeric(as.factor(ty$WEEK))		#rank(,ties.method= "min")
  }}
max(EGVU$visit)
EGVU$visit[EGVU$visit>5]<-5


### remove counts of >2 Adults and set to 2
hist(EGVU$n_adults)
EGVU$n_adults<-ifelse(EGVU$n_adults>2,2,EGVU$n_adults)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A CONTINUOUS NUMBER FOR YEAR AND TERRITORY (for random effects)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(EGVU)
EGVU$yearnum<-as.numeric(EGVU$Year)-2005
territories<-unique(EGVU$territory_name)
EGVU$terrnum<-match(EGVU$territory_name,territories)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP STATE MATRIX FOR KNOWN STATE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial values for the state variable: observed occurrence
# need to get ind in same order as y, otherwise init error because z and y order doesn't match

## REPLACED 'WEEK' with 'MONTH' throughout on 10 April 2018

enc.terrvis <- EGVU %>% filter (MONTH %in% c(4,5,6,7,8)) %>%
  select(terrnum, n_adults, PRIMOCC) %>%
  group_by(terrnum, PRIMOCC) %>%
  summarise(OBS=max(n_adults, na.rm=T)) %>%
  spread(key=PRIMOCC, value=OBS, fill = 0) %>%
  arrange(terrnum)

z.terrvis<-as.matrix(enc.terrvis[,c(2:dim(enc.terrvis)[2])])


# Create vector with last occasion of observation 
get.last.terrvis <- function(x) max(which(x==1))
f.terrvis <- apply(z.terrvis, 1, get.last.terrvis)

# Create vector with last occasion of observation of two birds
get.last2.terrvis <- function(x) max(which(x==2))
f2.terrvis <- apply(z.terrvis, 1, get.last2.terrvis)

# Ensure that state=1 until last observation
for (l in 1:dim(z.terrvis)[1]){
  if(f.terrvis[l]!="-Inf"){z.terrvis[l,(1:f.terrvis[l])]<-1}	# for all birds that were observed at least once set all occ to 1 before last observation
}

# Ensure that state=2 until last observation of two birds
for (l in 1:dim(z.terrvis)[1]){
  if(f2.terrvis[l]!="-Inf"){z.terrvis[l,(1:f2.terrvis[l])]<-2}	# for all birds that were observed at least once set all occ to 2 before last observation
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE DATA FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FOUR ARRAYS NEEDED - FOR OBSERVATIONS AND FOR EFFORT, FOR PARAMETER AND FOR INTERVAL
### EACH ARRAY HAS 3 dimensions for the number of sites (=territories), number of primary occasions and number of visits (=1 per day)


# PREPARING THE REQUIRED INPUT DATA:
R.terrvis<-length(unique(EGVU$terrnum))		# Number of individuals (= Territory-year combinations)
J.terrvis<-5						# Number of replicate surveys per month, changed from 7 before (WEEK)
K.terrvis<-dim(z.terrvis)[2]	# Number of primary occasions

# Set up some required arrays
site.terrvis <- 1:R.terrvis					# Sites
primoccs.terrvis <- 1:K.terrvis					# primary occasions (months across years)

# SET UP THE OBSERVATION DATA
y.terrvis <- array(NA, dim = c(R.terrvis, J.terrvis, K.terrvis))	# Detection histories
obseff.terrvis <- array(NA, dim = c(R.terrvis, J.terrvis, K.terrvis))	# observation effort


# CREATE A FULL DATA FRAME WITH ONE VALUE PER TERRITORY, PRIMARY AND SECONDARY OCCASION
## needs manual adjustment if number of primary occasions changes
fullEV.terrvis<-data.frame(terrnum=rep(unique(EGVU$terrnum),K.terrvis*J.terrvis),
                           PRIMOCC=rep(names(enc.terrvis)[2:(K.terrvis+1)],each=R.terrvis*J.terrvis),
                           visit=rep(rep(seq(1,J.terrvis,1),each=R.terrvis), K.terrvis))
fullEV.terrvis<-merge(fullEV.terrvis, EGVU, by=c("terrnum", "PRIMOCC", "visit"), all.x=T)
fullEV.terrvis$visit_duration[is.na(fullEV.terrvis$visit_duration)]<-0
fullEV.terrvis$primnum<-match(fullEV.terrvis$PRIMOCC,names(enc.terrvis))-1
head(fullEV.terrvis)



# FILL IN ARRAYS FOR OBSERVATIONS AND COVARIATES
for(k in 1:K.terrvis){
  
  dat.terrvis <- fullEV.terrvis %>% filter (primnum==k) %>%				## filter(WEEK==k+12) 
    select(terrnum, n_adults, visit) %>%
    group_by(terrnum, visit) %>%
    summarise(OBS=max(n_adults)) %>%
    spread(key=visit, value=OBS, fill = NA) %>%
    arrange(terrnum)
  
  effdat.terrvis <- fullEV.terrvis %>% filter (primnum==k) %>%
    select(terrnum, visit_duration, visit) %>%
    #select(terrnum, JDAY, visit) %>%
    group_by(terrnum, visit) %>%
    summarise(EFF=sum(visit_duration)) %>%
    #summarise(EFF=length(unique(JDAY))) %>%   ## this counts NA as 1 which is wrong
    spread(key=visit, value=EFF, fill = 0) %>%
    arrange(terrnum)
  
  y.terrvis[,,k] <-as.matrix(dat.terrvis[,2:6])
  obseff.terrvis[,,k] <-as.matrix(effdat.terrvis[,2:6])
}


# Standardize observation effort covariates
mean.eff.terrvis <- mean(obseff.terrvis, na.rm = TRUE)
sd.eff.terrvis <- sd(obseff.terrvis[!is.na(obseff.terrvis)])
obseff.terrvis <- (obseff.terrvis-mean.eff.terrvis)/sd.eff.terrvis     # Standardise observation effort



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP VECTOR OF TIME INTERVALS AND PARAMETER VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
primlookup<-data.frame(PRIMOCC=names(enc.terrvis)[2:(K.terrvis+1)],primnum=primoccs.terrvis) %>%
  mutate(YEAR=substring(PRIMOCC,1,4), MONTH=substring(PRIMOCC,6,6), DAY=15) %>%
  mutate(middate=ymd(paste(YEAR,MONTH,DAY)))
primlookup$nextdate<-c(primlookup$middate[2:dim(primlookup)[1]],NA)
primlookup<-primlookup %>%
  mutate(INTERVAL=interval(middate, nextdate) %/% months(1)) %>%
  mutate(PARM=ifelse(INTERVAL>5,2,1))

timeintervals.terrvis <- primlookup$INTERVAL[1:(K.terrvis-1)]	# time intervals between primary occasions
survparm.terrvis <- primlookup$PARM[1:(K.terrvis-1)]	# parameter for survival - breeding or non-breeding
yearindex.terrvis<-as.numeric(primlookup$YEAR)-2005








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INTEGRATED POPULATION MODEL FOR BALKAN EGYPTIAN VULTURES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2019_v1.jags")
cat("
model {
#-------------------------------------------------
# integrated population model for the balkan population
# - age structured model with 6 age classes: 
# - age-specific probability to recruit at ages 4 or 5
# - age-specific survival derived from tracking data
# - adult survival based on territory occupancy
# - pre breeding census, female-based assuming equal sex ratio & survival
# - productivity supplemented by captive-bred juveniles (0-10 per year)
#-------------------------------------------------


#-------------------------------------------------  
# 1. PRIORS FOR ALL DATA SETS
#-------------------------------------------------

# Priors and constraints FOR FECUNDITY
    ## this was originally on log scale with a density dependent offset
    # for (tpf in 1:T.count){  
    #   eta.fec[tpf]  ~ dnorm(0,0.01)T(-5,5)    ## inserted error as year-specific
    #   fecSE[tpf]<- mu.fec + eta.fec[tpf] ## removed the density dependence effect here and kept in normal scale, not in log scale as in original
    # } #t
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    

# Priors and constraints FOR POPULATION COUNTS OBSERVATION

    sigma.obs.count ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
    tau.obs.count<-pow(sigma.obs.count,-2)


# Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY

    for (iptel in 1:nind.telemetry){
      for (tptel in f.telemetry[iptel]:(n.occasions.telemetry-1)){
        phi.telemetry[iptel,tptel] <- beta.telemetry[x.telemetry[iptel,tptel]]
        p.telemetry[iptel,tptel] <- mean.p.telemetry
      } #tptel
    } #iptel

    for (utel in 1:3){					## change to 3 if > 2 years are modelled
      beta.telemetry[utel] ~ dunif(0, 1)              # Priors for age-specific survival
    }
    mean.p.telemetry ~ dunif(0.95, 1)                  # Prior for mean resighting (very high, only when tag malfunctions)



# Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    
    ## Priors for detection probability
    # Det prob varies for each site, det prob on logit scale, around 0.05 - 0.95 on normal scale
    for (ipterr in 1:nsite.terrvis){
      lmu.p.terrvis[ipterr] ~ dunif(-3, 2)				
    }							                        ## close k loop for det prob prior
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival for breeding and non-breeding season varies among years
    ## monthly BREEDING SEASON survival must be greater than 88% on logit scale
    ## monthly NON-BREEDING SEASON survival must be greater than 88% on logit scale

    for (nypterr in 1:nyears.terrvis){
      lm.phi.terrvis[1,nypterr] ~ dunif(2, 5)				      
      lm.phi.terrvis[2,nypterr] ~ dunif(2, 5)    
    }    
    

    ### RANDOM OBSERVATION EFFECT FOR EACH TERRITORY AND YEAR

    for (nyRpterr in 1:nyears.terrvis){
      rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    

    ### RANDOM SURVIVAL EFFECT FOR EACH TERRITORY AND YEAR

    for (ipterr2 in 1:nsite.terrvis){
      rand.surv.terrvis[ipterr2] ~ dnorm(0, tau.surv.terrvis)
    }
    tau.surv.terrvis <- 1 / (sd.surv.terrvis * sd.surv.terrvis)
    sd.surv.terrvis ~ dunif(0, 3)
    

# Priors for population process and immigration
  
    # Probability to start breeding at age = 4 or 5 years
    #reed.prop4 ~ dunif(0.02,0.03)                          ### based on Lieury et al. 2016, removed and fixed to 0.024 because slicer got stuck
    #reed.prop5 ~ dunif(0.120,0.130)                          ### based on Lieury et al. 2016, removed and fixed to 0.124 because slicer got stuck


    # Initial population sizes for first year of monitoring
    JUV[1] ~ dunif(50, 70)
    N1[1] ~ dunif(15, 35)
    N2[1] ~ dunif(10, 32)
    N3[1] ~ dunif(5, 20)
    N4[1] ~ dunif(5, 17)
    N5[1] ~ dunif(5, 15)
    N6[1] ~ dunif(100, 140)

    # JUV[1] ~ dpois(65)
    # N1[1] ~ dpois(25)
    # N2[1] ~ dpois(18)
    # N3[1] ~ dpois(13)
    # N4[1] ~ dpois(11)
    # N5[1] ~ dpois(10)
    # N6[1] ~ dpois(120)

#-------------------------------------------------  
# 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
#-------------------------------------------------
    
# -------------------------------------------------        
# 2.1. System process: female based matrix model
# -------------------------------------------------
    
for (tt in 2:T.count){

    nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
    JUV[tt] ~ dpois(nestlings[tt])                                                              ### need a discrete number otherwise dbin will fail
    N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(JUV[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
    N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
    N4[tt] ~ dbin(ann.surv.terrvis[tt], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(ann.surv.terrvis[tt], round(N4[tt-1]))                                                       ### number of 5-year old survivors
    N6[tt] ~ dbin(ann.surv.terrvis[tt], round((N5[tt-1]+N6[tt-1])))                                             ### number of 6-year or older (adult) birds

    # N1[tt] <- ann.phi.juv.telemetry*JUV[tt-1]                                               ### number of 1-year old survivors - add CAPT.ADD in here
    # N2[tt] <- ann.phi.sec.telemetry*N1[tt-1]                                                      ### number of 2-year old survivors
    # N3[tt] <- ann.phi.third.telemetry*N2[tt-1]                                                    ### number of 3-year old survivors
    # N4[tt] <- ann.surv.terrvis[tt]*N3[tt-1]                                                       ### number of 4-year old survivors
    # N5[tt] <- ann.surv.terrvis[tt]*N4[tt-1]                                                       ### number of 5-year old survivors
    # N6[tt] <- ann.surv.terrvis[tt]* (N5[tt-1]+N6[tt-1])                                           ### number of 6-year or older (adult) birds

} # tt
    


# -------------------------------------------------        
# 2.2. Observation process for population counts: state-space model of annual counts
# -------------------------------------------------

    for (tlc in 1:T.count){
      Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
      y.count[tlc] ~ dnorm(Nterr[tlc], tau.obs.count)								# Distribution for random error in observed numbers (counts)
    #}


# -------------------------------------------------        
# 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
# -------------------------------------------------
    #for (t in 1:nyear.fec){
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data




# -------------------------------------------------        
# 2.4. Likelihood for juvenile survival from telemetry
# -------------------------------------------------
 
  for (iltel in 1:nind.telemetry){

    # Define latent state at first capture
      z.telemetry[iltel,f.telemetry[iltel]] <- 1

    for (tltel in (f.telemetry[iltel]+1):n.occasions.telemetry){

    # State process
      z.telemetry[iltel,tltel] ~ dbern(mu1.telemetry[iltel,tltel])
      mu1.telemetry[iltel,tltel] <- phi.telemetry[iltel,tltel-1] * z.telemetry[iltel,tltel-1]

    # Observation process
      y.telemetry[iltel,tltel] ~ dbern(mu2.telemetry[iltel,tltel])
      mu2.telemetry[iltel,tltel] <- p.telemetry[iltel,tltel-1] * z.telemetry[iltel,tltel]

    } #tltel
  } #iltel


# -------------------------------------------------        
# 2.5. Likelihood for adult survival from territory monitoring
# -------------------------------------------------

    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL

    for (ilterr in 1:nsite.terrvis){
        z.terrvis[ilterr,1] ~ dbin(1,psi1.terrvis[ilterr])
        phi.terrvis[ilterr,1] <- 1
        surv.terrvis[ilterr,1] <- 1
        linphi.terrvis[ilterr,1] <- 0
    
      for (klterr in 2:nprim.terrvis){
        linphi.terrvis[ilterr,klterr] <-lm.phi.terrvis[parm.terrvis[klterr-1],year.terrvis[klterr]] + rand.surv.terrvis[ilterr] 
        phi.terrvis[ilterr,klterr] <- 1 / (1 + exp(-linphi.terrvis[ilterr,klterr]))
        surv.terrvis[ilterr,klterr] <- pow(phi.terrvis[ilterr,klterr],intv.terrvis[klterr-1])
        z.terrvis[ilterr,klterr] ~ dbin(phi.terrvis[ilterr,klterr],z.terrvis[ilterr,klterr-1])
      } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    ## observation prob depends on number of visits per week

    for (iobs in 1:nsite.terrvis){
      for (jobs in 1:nrep.terrvis){
        for (kobs in 1:nprim.terrvis){

          linpred.terrvis[iobs,jobs,kobs] <- lmu.p.terrvis[iobs]+beta.obs.eff.terrvis*eff.terrvis[iobs,jobs,kobs] + rand.obs.terrvis[year.terrvis[kobs]]
          p.terrvis[iobs,jobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,jobs,kobs]))
          y.terrvis[iobs,jobs,kobs] ~ dbin(p.terrvis[iobs,jobs,kobs], z.terrvis[iobs,kobs])
    
    
        } 						# close kobs loop over primary periods
      } 						 # close jobs loop over week visits
    } 							# close iobs loop over sites
    


# -------------------------------------------------        
# 3. DERIVED PARAMETERS
# -------------------------------------------------

### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES

  ann.phi.juv.telemetry<-pow(beta.telemetry[1],3)*pow(beta.telemetry[2],9)						### adjust here if first year consists of two parameters
  ann.phi.sec.telemetry<-pow(beta.telemetry[2],12)
  ann.phi.third.telemetry<-pow(beta.telemetry[3],12)


### 3.2 TERRITORY MONITORING DERIVED SURVIVAL AND RESIGHTING ESTIMATES

    # Detection probability averaged across all sites and occasions
    mean.p.terrvis <-mean(p.terrvis[,,])
    
    # Survival probability averaged across season and year
    for (nys in 1:nyears.terrvis){
      breed.surv.terrvis[nys]<- 1 / (1 + exp(-lm.phi.terrvis[1,nys]))
      nonbreed.surv.terrvis[nys]<- 1 / (1 + exp(-lm.phi.terrvis[2,nys]))
      ann.surv.terrvis[nys]<-pow(breed.surv.terrvis[nys],4)*pow(nonbreed.surv.terrvis[nys],7)
    }

### 3.3 POPULATION GROWTH DERIVED FROM COUNTS    

    # Annual population growth rate
    for (ipop in 1:(T.count-1)){
      lambda.t[ipop] <- Nterr[ipop+1] / max(Nterr[ipop],1)           ### inserted max of 1 to avoid this becoming an invalid parent error
      loglambda.t[ipop]<-log(lambda.t[ipop])## for calculating geometric mean of overall population growth rate
    }
    
    #### OVERALL POPULATION GROWTH RATE  #########
    mean.lambda<-exp((1/(T.count-1))*sum(loglambda.t[1:(T.count-1)]))   # Geometric mean




# -------------------------------------------------        
# 4. PREDICTION INTO THE FUTURE
# -------------------------------------------------

    # CAPTIVE RELEASE OF JUVENILE BIRDS
    CAPT.ADD ~ dunif(0,15)

    #eta.fec.fut  <- mean(eta.fec[1:T.count])    ## draw random error for fecundity based on past error
    #fecSE.fut<- mu.fec + eta.fec.fut
    fut.survival <-max(ann.surv.terrvis[1:nyears.terrvis])

    ## POPULATION PROCESS

for (fut in (T.count+1):(T.count+PROJECTION)){
  
    #eta.fec[fut]  ~ dunif(min(eta.fec[1:T.count]),max(eta.fec[1:T.count]))    ## draw random error for fecundity based on past error
    #fecSE[fut]<- mu.fec + eta.fec[fut]
    
    nestlings[fut] <- mu.fec * 0.5 * Nterr[fut]                                                              ### number of local recruits
    JUV[fut] ~ dpois(nestlings[fut])                                                              ### need a discrete number otherwise dbin will fail

    N1[fut] <- ann.phi.juv.telemetry*(JUV[fut-1]+CAPT.ADD)                                          ### number of 1-year old survivors - add CAPT.ADD in here
    N2[fut] <- ann.phi.sec.telemetry*N1[fut-1]                                                      ### number of 2-year old survivors
    N3[fut] <- ann.phi.third.telemetry*N2[fut-1]                                                    ### number of 3-year old survivors

    #fut.survival[fut] ~ dunif(min(ann.surv.terrvis[1:nyears.terrvis]),max(ann.surv.terrvis[1:nyears.terrvis]))
    N4[fut] <- fut.survival*N3[fut-1]                                                       ### number of 4-year old survivors
    N5[fut] <- fut.survival*N4[fut-1]                                                       ### number of 5-year old survivors
    N6[fut] <- fut.survival* (N5[fut-1]+N6[fut-1])                                           ### number of 6-year or older (adult) birds
    Nterr[fut] <- N4[fut] * 0.024 + N5[fut] * 0.124 + N6[fut]
    } # fut

   for (fut2 in (T.count+1):(T.count+PROJECTION-1)){
      lambda.t[fut2] <- Nterr[fut2+1] / max(Nterr[fut2],1)
      loglambda.t[fut2]<-log(lambda.t[fut2])## for calculating geometric mean of overall population growth rate
    
      } # fut2


    #### FUTURE POPULATION GROWTH RATE  #########
    fut.lambda<-exp((1/(PROJECTION-1))*sum(loglambda.t[(T.count+1):(T.count+PROJECTION-1)]))   # Geometric mean
    



} # END of the JAGS model
",fill = TRUE)
sink()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DATA AND INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Bundle data across all of the data sources

INPUT <- list(y.terrvis = y.terrvis,
              nsite.terrvis = R.terrvis,
              nrep.terrvis = J.terrvis,
              nprim.terrvis = K.terrvis,
              nyears.terrvis=max(yearindex.terrvis),
              eff.terrvis=obseff.terrvis,
              psi1.terrvis=z.terrvis[,1],
              intv.terrvis=timeintervals.terrvis,
              parm.terrvis=survparm.terrvis,
              year.terrvis=yearindex.terrvis,
              
              y.count=trendinput$N,
              T.count=length(trendinput$N),		## year is standardized so that covariate values are not too far away from zero

              R.fec=breedinput$R,
              J.fec=breedinput$J,

              y.telemetry = y.telemetry,
              f.telemetry = f.telemetry,
              nind.telemetry = dim(CH.telemetry)[1],
              n.occasions.telemetry = dim(y.telemetry)[2],
              z.telemetry = z.telemetry,
              x.telemetry = x.telemetry,
              
              PROJECTION=10)



## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("Nterr", "lambda.t","ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry",      #"breed.prop4","breed.prop5",
            "ann.surv.terrvis","mean.p.terrvis","mean.lambda","CAPT.ADD","fut.lambda","mu.fec")
          


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### provided by Adam Butler on 10 April 2018

ifix <- function(w){w[ w == -Inf] <- NA ; w }
ymax.terrvis <- ifix(apply(y.terrvis, c(1,3), max, na.rm=TRUE))
f1.terrvis <- ifix(apply(ymax.terrvis, 1, function(z.terrvis){max(which(z.terrvis == 1))}))
f2.terrvis <- ifix(apply(ymax.terrvis, 1, function(z.terrvis){max(which(z.terrvis == 2))}))

## any(is.na(f1) & is.na(f2)) ## FALSE ## check - good !

z.obs.terrvis <- z.init.terrvis <- array(dim=dim(ymax.terrvis), data = NA)
nyr <- ncol(ymax.terrvis)

for(k in 1:nrow(ymax.terrvis)){
  
  ## ######################################
  ## Set observed values
  
  if(! is.na(f2.terrvis[k])){  
    z.obs.terrvis[k,1:(f2.terrvis[k])] <- 2
  }
  else{     
    z.obs.terrvis[k,1:(f1.terrvis[k])] <- 1
  }
  
  ## ######################################
  ## Set initial values for stochastic nodes
  
  test.a <- (! is.na(f2.terrvis[k])) & (f2.terrvis[k] < nyr) & ((is.na(f1.terrvis[k])) | (f1.terrvis[k] <= f2.terrvis[k]))
  if(test.a){
    z.init.terrvis[k,(f2.terrvis[k]+1):nyr] <- 0
  }
  
  ## ###########
  
  test.b <- (! is.na(f1.terrvis[k])) & (f1.terrvis[k] < nyr) & (is.na(f2.terrvis[k]) | (f1.terrvis[k] > f2.terrvis[k]))
  if(test.b){
    z.init.terrvis[k,(f1.terrvis[k]+1):nyr] <- 0
  }
  
  ## ###########
  
  test.c <- (! is.na(f2.terrvis[k])) & (! is.na(f1.terrvis[k])) & (f1.terrvis[k] > f2.terrvis[k])
  if(test.c){
    z.init.terrvis[k,(f2.terrvis[k]+1):f1.terrvis[k]] <- 1
  }
  
  ## ###########
}

## ######################################################


### MODEL WITH RANDOM EFFECTS
INPUT$z.terrvis<-z.obs.terrvis
initIPM <- function(){list(z.terrvis = z.init.terrvis,
                         lmu.p.terrvis=runif(dim(z.terrvis)[1],-3, 2),
                         lm.phi.terrvis=matrix(runif(2*max(yearindex.terrvis),2, 5),ncol=max(yearindex.terrvis)),
                         #breed.prop4= 0.024,  
                         #breed.prop5=0.124,
                         sigma.obs.count=runif(1,0,10),
                         mu.fec = runif(1,0,2),
                         z.telemetry = cjs.init.z(INPUT$y.telemetry, INPUT$f.telemetry),
                         beta.telemetry = runif(3, 0, 1),
                         mean.p.telemetry = runif(1, 0.95, 1))}  




# MCMC settings
nc <- 3
nt <- 5
ni <- 100000
nb <- 50000 


# RUN THE MODEL ALLOWING FOR RANDOM ERRORS

rm(NeoIPMi)
for (i in 1:50){

try(NeoIPMi <- jags(data=INPUT,
                inits=initIPM,
                parameters.to.save=paraIPM,
                model.file="S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\EGVU_IPM_2018_v3.jags",
                n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T), silent=T)
if("NeoIPMi" %in% ls())break
  
}


NeoIPMi





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT THE OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load("EGVU_IPM_output.RData")
out<-as.data.frame(NeoIPMi$summary)
out$parameter<-row.names(NeoIPMi$summary)
write.table(out,"EGVU_IPM_estimates_v1.csv", sep=",", row.names=F)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT POPULATION TREND 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EV.ssm<-out[(grep("Nterr",out$parameter)),c(1,3,7)] %>%
  mutate(Year=c(trendinput$year,seq(2018,2027,1)))
names(EV.ssm)[1:3]<-c('mean','lcl','ucl')

## DRAW FIGURE PLOTTING MODEL OUTPUT AND RAW INPUT DATA



#pdf("EV_population_trend_all_Balkan.pdf", width=9, height=6)
#postscript("Fig1_Balkan.eps", width=9, height=6)
#jpeg("Fig1_Balkan.jpg", width=9, height=6, units="in", res=600, quality=100)
#par(oma=c(0,0,0,0),mar=c(4.2,4.5,0,0.5), cex=1.2)
ggplot()+
  geom_line(data=EV.ssm, aes(x=Year, y=mean), size=1,col='darkred')+
  geom_ribbon(data=EV.ssm,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_point(data=EV.ssm, aes(x=Year-0.1, y=mean), size=1,pch=16)+
  geom_point(data=trendinput, aes(x=year+0.1, y=N), size=1,col='darkgrey')+
  
  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,130),breaks=seq(0,130,30), labels=as.character(seq(0,120,30)))+
  scale_x_continuous(name="Year", breaks=seq(2006,2027,2), labels=as.character(seq(2006,2027,2)))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT THE SAMPLES WHERE FUTURE LAMBDA IS POSITIVE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samplesout<-as.data.frame(rbind(NeoIPMi$samples[[1]],NeoIPMi$samples[[2]],NeoIPMi$samples[[3]],NeoIPMi$samples[[4]]))



#### CALCULATE PROPORTION OF SIMULATIONS WHERE LAMBDA >1 FOR EACH CAPT.ADD
plotdat<-data.frame()
for (i in 1:nc){
  samplesout<-as.data.frame(NeoIPMi$samples[[i]])
  x<-samplesout %>% select(fut.lambda,CAPT.ADD,mu.fec) %>%
  mutate(EGVU=round(CAPT.ADD,0)) %>%
  mutate(count=1) %>%
  mutate(stablepop=ifelse(fut.lambda>0.99999,1,0)) %>%
  group_by(EGVU) %>%
  summarise(prop_stable=sum(stablepop)/sum(count))
plotdat<-rbind(plotdat,x)
}
write.table(plotdat,"EGVU_IPM_samples_for_plot.csv", sep=",", row.names=F)


#plotdat %>% arrange(CAPT.ADD) %>%
#  mutate(cumprob=seq(1:length(CAPT.ADD))/length(CAPT.ADD)) %>%
  
  
  ggplot(plotdat)+
  geom_point(aes(x=EGVU,y=prop_stable,width=2))+
  geom_smooth(aes(x=EGVU,y=prop_stable), method='gam')+   
  ylab("Probability of stable or increasing population") +
  xlab("Number of captive-bred juveniles added per year") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PRODUCE OUTPUT REPORT
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  ### create HTML report for overall summary report
  Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")
  
  rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\EGVU_CaptRelease.Rmd',
                    output_file = "EGVU_pop_model_captive_breeding.html",
                    output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel')
  
  
  
  
  
  
  
  




####### abandoned histogram #######

samplesout %>% select(fut.lambda,CAPT.ADD,mu.fec) %>%
  #filter(mu.fec<1.10) %>%
  filter(fut.lambda>0.99999) %>%
  #summarise(NEEDED=min(CAPT.ADD))
  
  ggplot()+
  geom_histogram(aes(x=CAPT.ADD,y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=1)+    
  ylab("Probability of stable or increasing population") +
  xlab("Number of captive-bred juveniles added per year") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
  
save.image("EGVU_IPM_output.RData")
