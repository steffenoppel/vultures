##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS
#
##########################################################################
# based on Lieury et al. 2015
# modified by Steffen oppel 4 April 2018
# goal is to assess population viability with and without captive released juveniles on Balkans

# pasted together data prep and model from various files on 16 May 2018
# update 17 May 2018: trying to troubleshoot dbinom error
# updated 22 May 2018: breed.prop cannot be a stochastic node, because 'slicer' gets stuck - tried with multiple distributions, priors and initials, but never works
# included Adam's suggestions for dbinom and output treatment (i.e. abandoned the idea of solving for CAPT.ADD)
# updated 23 May 2018: changed priors for fecundity as output was essentially 0-2

#### MAJOR REVISION after Aqaba meeting in April 2019
# included data from 2019
# changed survival probability of released birds to be similar to second-year birds (not juveniles, due to delayed release)
# included scenarios for future with n captive released birds and improvement in survival

# update on 1 May 2019: suggestion from Volen to also include 9 'rescued' chicks taken from the wild but rehabilitated and released with higher survival probability
# update 3 May: future lambda lower than mean lambda - included binomial draws for survival for future projection

# update on 26 July 2019: re-run model and provided extinction probabilities for report

## MAJOR UPDATE 26 July 2019: include egg reduction to harvest second eggs
## calculate fecundity when only 1 egg available (reduction from normal fecundity)
## add harvested eggs as additional juveniles to be released
## updated 27 July 2019: fixed fecundity accordingly in scenarios 


library(readxl)
library(jagsUI)
library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE POPULATION MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
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


### MODIFY BREEDINPUT FOR EGG HARVEST ###

breedinput1EGG<- breed %>% filter(Year>2005) %>%
  filter(!is.na(breed_success)) %>%
  mutate(count=1) %>%
  mutate(fledglings=ifelse(is.na(fledglings),0,fledglings)) %>%
  mutate(fledglings=ifelse(fledglings==2,1,fledglings)) %>%     ### remove the second egg and fledgling
  group_by(Year) %>%
  summarise(R=sum(count), J=sum(fledglings))

breedinput$J-breedinput1EGG$J

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE JUVENILE TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE)
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
  filter(!Name=="Zighmund") %>%            ### remove single bird that was never free flying
  filter(origin=="wild")                    ### USE ONLY WILD JUVENILES

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

## for reporting - indicate months of death for juveniles
get.last.telemetry<-function(x)min(which(x==0))
hist(apply(y.telemetry,1,get.last.telemetry))


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


#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_surveys.r")), wait = TRUE, invisible = FALSE)
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

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2019_EggRem5_NoRescue.jags")
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
# - fecundity reduced by harvesting second eggs
#-------------------------------------------------


#-------------------------------------------------  
# 1. PRIORS FOR ALL DATA SETS
#-------------------------------------------------

# Priors and constraints FOR FECUNDITY
    mu.fec[1] ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    mu.fec[2] ~ dunif(0,1)       # Priors on fecundity after removing 2nd egg can range from 0-1 chick per pair (uninformative)    

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
  
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(50, 70)   ##changed from JUV
    N1[1] ~ dunif(15, 35)
    N2[1] ~ dunif(10, 32)
    N3[1] ~ dunif(5, 20)
    N4[1] ~ dunif(5, 17)
    N5[1] ~ dunif(5, 15)
    N6[1] ~ dunif(100, 140)


#-------------------------------------------------  
# 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
#-------------------------------------------------
    
# -------------------------------------------------        
# 2.1. System process: female based matrix model
# -------------------------------------------------
    
for (tt in 2:T.count){

    nestlings[tt] <- mu.fec[1] * 0.5 * Nterr[tt]                                                              ### number of local recruits
    #JUV[tt] ~ dpois(nestlings[tt])                                                              ### need a discrete number otherwise dbin will fail
    N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
    N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
    N4[tt] ~ dbin(ann.surv.terrvis[tt], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(ann.surv.terrvis[tt], round(N4[tt-1]))                                                ### number of 5-year old survivors
    N6[tt] ~ dbin(ann.surv.terrvis[tt], round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds

} # tt
    


# -------------------------------------------------        
# 2.2. Observation process for population counts: state-space model of annual counts
# -------------------------------------------------

    for (tlc in 1:T.count){
      Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
      y.count[tlc] ~ dnorm(Nterr[tlc], tau.obs.count)								# Distribution for random error in observed numbers (counts)
    #}


# -------------------------------------------------        
# 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods (or binom for single broods after egg harvest)
# -------------------------------------------------
    #for (t in 1:nyear.fec){
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec[1]
      J.fec.red[tlc] ~ dbin(mu.fec[2],R.fec[tlc])
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
### INCLUDE SCENARIOS FOR N CAPTIVE BIRDS AND SURVIVAL IMPROVEMENT
### INCLUDE THE RESCUE OF 9 chicks and release with increased survival
### FUTURE FECUNDITY IS BASED ON mu.fec[1] for no action scenario and mu.fec[2] when there is captive-released birds


    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){

      #fut.fec[ncr] <-mu.fec[min(capt.release[ncr]+1,2)]  ## this will be mu.fec[1] for 0 capt.release and mu.fec[2] when captive birds are released
      fut.fec[1:5] <-mu.fec[2]  ## this will be mu.fec[1] for 0 capt.release and mu.fec[2] when captive birds are released
      fut.fec[5:PROJECTION] <-mu.fec[1]  ## this will be mu.fec[1] for 0 capt.release and mu.fec[2] when captive birds are released

      # SPECIFY IMPROVEMENT OF SURVIVAL
      for (is in 1:scen.imp.surv){ 
    
        # improve.surv ~ dunif(1,1.10)
        fut.survival[ncr,is] <-min(imp.surv[is]*mean(ann.surv.terrvis[1:nyears.terrvis]),1) ### invalid parent error if survival>1

        ## POPULATION PROCESS
        ### DOES NOT WORK WITH SAME ARRAY OF POP TRAJECTORY AS 3 dimensions required
        ## need to copy previous array elements
        nestlings.f[ncr,is,1] <- round(nestlings[T.count])   ##JUV[T.count]
        N1.f[ncr,is,1] <- N1[T.count]
        N2.f[ncr,is,1] <- N2[T.count]
        N3.f[ncr,is,1] <- N3[T.count]
        N4.f[ncr,is,1] <- N4[T.count]
        N5.f[ncr,is,1] <- N5[T.count]
        N6.f[ncr,is,1] <- N6[T.count]
        Nterr.f[ncr,is,1] <- Nterr[T.count]
        rescued[ncr,is,1] <- 0                ## we need to fill in a value to make the matrix complete, this is not actually used in any calculation

          for (fut in 2:PROJECTION){

            ### probabilistic formulation
            #rescued[ncr,is,fut] ~ dpois(5) T(1,11)                                         ### number of chicks rescued and rehabilitated to improve survival FROM HARVESTING SECOND EGGS
            nestlings.f[ncr,is,fut] <- (fut.fec[ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
            N1.f[ncr,is,fut] ~ dbin(min((imp.surv[is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]+capt.release[ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
            N2.f[ncr,is,fut] ~ dbin(min((imp.surv[is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))  ### number of 2-year old survivors
            N3.f[ncr,is,fut] ~ dbin(min((imp.surv[is]*ann.phi.third.telemetry),1),round(N2.f[ncr,is,fut-1]))                                                    ### number of 3-year old survivors
            N4.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is],round(N3.f[ncr,is,fut-1]))                                                       ### number of 4-year old survivors
            N5.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is],round(N4.f[ncr,is,fut-1]))                                                       ### number of 5-year old survivors
            N6.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is],round((N5.f[ncr,is,fut-1]+N6.f[ncr,is,fut-1])))                                   ### number of 6-year or older (adult) birds
            Nterr.f[ncr,is,fut] <- round((N4.f[ncr,is,fut] * 0.024) + (N5.f[ncr,is,fut] * 0.124) + (N6.f[ncr,is,fut]))


          } # fut

          for (fut2 in 1:(PROJECTION-1)){
            lambda.t.f[ncr,is,fut2] <- Nterr.f[ncr,is,fut2+1] / max(Nterr.f[ncr,is,fut2],1)
            loglambda.t.f[ncr,is,fut2]<-log(lambda.t.f[ncr,is,fut2])## for calculating geometric mean of overall population growth rate
          } # fut2


          #### FUTURE POPULATION GROWTH RATE  #########
          fut.lambda[ncr,is]<-exp((1/(PROJECTION-1))*sum(loglambda.t.f[ncr,is,1:(PROJECTION-1)]))   # Geometric mean

      } # end scenario of imp.surv

    } # end scenario of n capt released

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
              J.fec.red=breedinput1EGG$J,   ## added to model fecundity when second egg is removed

              y.telemetry = y.telemetry,
              f.telemetry = f.telemetry,
              nind.telemetry = dim(CH.telemetry)[1],
              n.occasions.telemetry = dim(y.telemetry)[2],
              z.telemetry = z.telemetry,
              x.telemetry = x.telemetry,
              
              ### Future Projection and SCENARIOS FOR TRAJECTORY
              PROJECTION=30,
              scen.capt.release=4,
              scen.imp.surv=7,
              capt.release=c(0,2,4,6),
              imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12))



## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("lambda.t","ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry",      #"breed.prop4","breed.prop5",
            "ann.surv.terrvis","mean.p.terrvis","mean.lambda","fut.lambda","mu.fec","Nterr", "Nterr.f")
          


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
                         sigma.obs.count=runif(1,0,10),
                         mu.fec = runif(2,0,1),
                         z.telemetry = cjs.init.z(INPUT$y.telemetry, INPUT$f.telemetry),
                         beta.telemetry = runif(3, 0, 1),
                         mean.p.telemetry = runif(1, 0.95, 1))}  


# MCMC settings
nc <- 4
nt <- 4
ni <- 50000
nb <- 10000


# RUN THE MODEL ALLOWING FOR RANDOM ERRORS

# rm(NeoIPMi)
# for (i in 1:50){
#   try(
NeoIPMbasic <- jags(data=INPUT,
                inits=initIPM,
                parameters.to.save=paraIPM,
                model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_fut_scenarios.jags",
                n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
  
NeoIPMeggredNoRescue <- jags(data=INPUT,
                       inits=initIPM,
                       parameters.to.save=paraIPM,
                       model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_EggRem5_NoRescue.jags",
                       n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

NeoIPMeggredRescue <- jags(data=INPUT,
                    inits=initIPM,
                    parameters.to.save=paraIPM,
                    model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_fut_scenarios_v2.jags",
                    n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
  
# , silent=T)
# if("NeoIPMi" %in% ls())break
# }


NeoIPMi$samples

save.image("EGVU_IPM_output2019_v2.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT THE OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)
#load("EGVU_IPM_output2019_chicks.RData")
out<-as.data.frame(NeoIPMi$summary)
out$parameter<-row.names(NeoIPMi$summary)
#write.table(out,"EGVU_IPM_estimates_v3.csv", sep=",", row.names=F)

out<-as.data.frame(NeoIPMeggred$summary)
out$parameter<-row.names(NeoIPMeggred$summary)
#write.table(out,"EGVU_IPM_future_v2.csv", sep=",", row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR REPORT /MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.lambda','ann.phi.juv.telemetry','ann.phi.sec.telemetry','ann.phi.third.telemetry','mu.fec')) %>%
  select(parameter,c(5,3,7))


## retrieve the adult survival estimates averaged across all years
selcol<-grep("ann.surv.terrvis",dimnames(NeoIPMeggred$samples[[1]])[[2]])
ann.surv.terrvis<-numeric()
for (c in 1:nc){
  ann.surv.terrvis<-c(ann.surv.terrvis,as.numeric(NeoIPMeggred$samples[[c]][,selcol]))
}

TABLE1[6,]<-c("adult survival",quantile(ann.surv.terrvis,0.5),quantile(ann.surv.terrvis,0.025),quantile(ann.surv.terrvis,0.975))
#write.table(TABLE1,"clipboard", sep="\t", row.names=F)


## retrieve the future lambdas
selcol<-grep("fut.lambda",dimnames(NeoIPMeggred$samples[[1]])[[2]])
fut.lambda<-NeoIPMeggred$samples[[1]][,selcol]
for (c in 2:nc){
  fut.lambda<-rbind(fut.lambda,NeoIPMeggred$samples[[c]][,selcol])
}

## give the projections proper scenario labels
capt.release=c(0,2,4,6)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12)

FUTLAM<-as.data.frame(fut.lambda) %>% gather(key="parm",value="f.lam") %>%
  group_by(parm) %>%
  summarise(median=quantile(f.lam,0.5),lcl=quantile(f.lam,0.025),ucl=quantile(f.lam,0.975)) %>%
  mutate(capt.release=0, imp.surv=0) 
FUTLAM$capt.release[grep(",",FUTLAM$parm)]<-capt.release[as.numeric(substr(FUTLAM$parm[grep(",",FUTLAM$parm)],12,12))]
FUTLAM$imp.surv[grep(",",FUTLAM$parm)]<-imp.surv[as.numeric(substr(FUTLAM$parm[grep(",",FUTLAM$parm)],14,14))]
FUTLAM$imp.surv <- ifelse(FUTLAM$imp.surv>1,paste("+",as.integer((FUTLAM$imp.surv-1)*100),"%"),"none")
FUTLAM
#write.table(FUTLAM[,c(1,5,6,2,3,4)],"clipboard", sep="\t", row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT POPULATION TREND FOR DIFFERENT SCENARIOS OF IMPROVEMENT AND CAPT RELEASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the population projections
EV.fut<-out[(grep("Nterr",out$parameter)),c(12,1,3,7)] %>%
  mutate(Year=c(trendinput$year,rep(seq(2019,2048,1),each=4*7)))
names(EV.fut)[1:4]<-c('parm','mean','lcl','ucl')

## give the projections proper scenario labels
capt.release=c(0,2,4,6)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12)

EV.fut <- EV.fut %>% mutate(capt.release=0, imp.surv=0)
EV.fut$capt.release[grep(",",EV.fut$parm)]<-capt.release[as.numeric(substr(EV.fut$parm[grep(",",EV.fut$parm)],9,9))]
EV.fut$imp.surv[grep(",",EV.fut$parm)]<-imp.surv[as.numeric(substr(EV.fut$parm[grep(",",EV.fut$parm)],11,11))]
EV.fut$imp.surv <- ifelse(EV.fut$imp.surv>1,paste("+",as.integer((EV.fut$imp.surv-1)*100),"%"),"none")
              
EV.fut[20:30,]  

## add past trajectory for the four capt.release scenarios
EV.fut<-rbind(EV.fut,(EV.fut %>% filter(Year<2019) %>% mutate(capt.release=2)),
              (EV.fut %>% filter(Year<2019) %>% mutate(capt.release=4)),
              (EV.fut %>% filter(Year<2019) %>% mutate(capt.release=6)))

## create factors for plot labels and order them appropriately
EV.fut$capt.release <- ifelse(EV.fut$capt.release>1,paste("+",EV.fut$capt.release,"chicks/year"),"no captive releases")
EV.fut$capt.release <- factor(EV.fut$capt.release, levels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
EV.fut$imp.surv <- factor(EV.fut$imp.surv, levels = c("+ 12 %","+ 10 %","+ 8 %","+ 6 %","+ 4 %","+ 2 %","none"))



### produce plot with 4 panels and multiple lines per year

#pdf("EV_population_projection_C3.pdf", width=10, height=7)
#postscript("Fig1_Balkan.eps", width=9, height=6)
#jpeg("Fig1_Balkan.jpg", width=9, height=6, units="in", res=600, quality=100)
#par(oma=c(0,0,0,0),mar=c(4.2,4.5,0,0.5), cex=1.2)
ggplot()+
  geom_line(data=EV.fut, aes(x=Year, y=mean, color=imp.surv), size=1)+
  geom_ribbon(data=EV.fut,aes(x=Year, ymin=lcl,ymax=ucl, fill=imp.surv),alpha=0.2)+
  geom_point(data=trendinput, aes(x=year+0.1, y=N), size=1,col='darkgrey')+
  facet_wrap(~capt.release,ncol=2) +
  
  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,300),breaks=seq(0,300,50), labels=as.character(seq(0,300,50)))+
  scale_x_continuous(name="Year", breaks=seq(2006,2048,5), labels=as.character(seq(2006,2048,5)))+
  guides(color=guide_legend(title="Survival increase"),fill=guide_legend(title="Survival increase"))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

dev.off()



## specify what survival should be for stable population
mean(out[(grep("surv",out$parameter)),1])*1.06
mean(out[(grep("surv",out$parameter)),1])*1.08



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUANTIFY PROPORTION OF SIMULATIONS (=PROBABILITY) THAT FUTURE LAMBDA IS POSITIVE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samplesout<-as.data.frame(rbind(NeoIPMeggred$samples[[1]],NeoIPMeggred$samples[[2]],NeoIPMeggred$samples[[3]]))
head(samplesout)

## give the projections proper scenario labels
capt.release=c(0,2,4,6)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12)

allsamp <- samplesout %>% gather(key="parm", value="value") %>%
  filter(grepl("fut.lambda",parm)) %>%
  mutate(capt.release=0, imp.surv=0)

allsamp$capt.release[grep(",",allsamp$parm)]<-capt.release[as.numeric(substr(allsamp$parm[grep(",",allsamp$parm)],12,12))]
allsamp$imp.surv[grep(",",allsamp$parm)]<-imp.surv[as.numeric(substr(allsamp$parm[grep(",",allsamp$parm)],14,14))]
allsamp$imp.surv <- ifelse(allsamp$imp.surv>1,paste("+",as.integer((allsamp$imp.surv-1)*100),"%"),"none")

## create factors for plot labels and order them appropriately
allsamp$capt.release <- ifelse(allsamp$capt.release>1,paste("+",allsamp$capt.release,"chicks/year"),"no captive releases")
allsamp$capt.release <- factor(allsamp$capt.release, levels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
allsamp$imp.surv <- factor(allsamp$imp.surv, levels = c("+ 12 %","+ 10 %","+ 8 %","+ 6 %","+ 4 %","+ 2 %","none"))


ggplot(allsamp)+
  #geom_point(aes(x=imp.surv,y=value, color=capt.release),size=0.5)+
  geom_hline(aes(yintercept=1), color='red', size=2)+
  geom_violin(aes(x=imp.surv,y=value)) +
  
  facet_wrap(~capt.release,ncol=2) +
  
  ## format axis ticks
  guides(color=guide_legend(title="Survival increase"),fill=guide_legend(title="Survival increase"))+
  
  
  ylab("Future population trend") +
  xlab("Improvement in survival probability") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




### quantify what proportion of simulations result in future lambda >1

allsamp %>% select(capt.release,imp.surv,value) %>%
  mutate(n=1, inc=ifelse(value>0.9999,1,0)) %>%
  group_by(imp.surv,capt.release) %>%
  summarise(prop.increase=sum(inc)/sum(n))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUANTIFY PROPORTION OF SIMULATIONS (=PROBABILITY) THAT FUTURE LAMBDA IS POSITIVE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(samplesout)
allsamp <- samplesout %>% gather(key="parm", value="value") %>%
  filter(grepl("Nterr.f",parm)) %>%
  filter(grepl(",30]",parm)) %>%   ### only use final year
  mutate(capt.release=0, imp.surv=0)

allsamp$capt.release[grep(",",allsamp$parm)]<-capt.release[as.numeric(substr(allsamp$parm[grep(",",allsamp$parm)],9,9))]
allsamp$imp.surv[grep(",",allsamp$parm)]<-imp.surv[as.numeric(substr(allsamp$parm[grep(",",allsamp$parm)],11,11))]
allsamp$imp.surv <- ifelse(allsamp$imp.surv>1,paste("+",as.integer((allsamp$imp.surv-1)*100),"%"),"none")

## create factors for plot labels and order them appropriately
allsamp$capt.release <- ifelse(allsamp$capt.release>1,paste("+",allsamp$capt.release,"chicks/year"),"no captive releases")
allsamp$capt.release <- factor(allsamp$capt.release, levels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
allsamp$imp.surv <- factor(allsamp$imp.surv, levels = c("+ 12 %","+ 10 %","+ 8 %","+ 6 %","+ 4 %","+ 2 %","none"))



### PLOT EXTINCTION RISK
pdf("EV_extinction_risk2028.pdf", width=10, height=7)
ggplot(allsamp)+
  geom_hline(aes(yintercept=10), color='red', size=2)+
  geom_violin(aes(x=imp.surv,y=value)) +
  
  facet_wrap(~capt.release,ncol=2) +
  
  ## format axis ticks
  guides(color=guide_legend(title="Survival increase"),fill=guide_legend(title="Survival increase"))+
  
  
  ylab("Population size in 2028") +
  xlab("Improvement in survival probability") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=16, color="black", vjust=0.5), 
        axis.title=element_text(size=16), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
dev.off()

### quantify what proportion of simulations result in future pop <10 birds (=extinct)

ext.probs<-allsamp %>% select(capt.release,imp.surv,value) %>%
  mutate(n=1, inc=ifelse(value<10,1,0)) %>%
  group_by(imp.surv,capt.release) %>%
  summarise(ext.prob=sum(inc)/sum(n)) %>%
  arrange(desc(ext.prob))

fwrite(ext.probs,"EGVU_extinction.prob_2028.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT PROBABILITY OF EXTINCTION OVER TIME 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(samplesout)
extprop <- samplesout %>% gather(key="parm", value="value") %>%
  filter(grepl("Nterr.f",parm)) %>%
  mutate(capt.release=as.numeric(substr(parm,9,9)), imp.surv=as.numeric(substr(parm,11,11)))  %>%
  filter(imp.surv %in% c(1,2,3,4)) %>% ### reduce plotting options
  mutate(Year=ifelse(nchar(parm)==14,substr(parm,13,13),substr(parm,13,14))) %>%
  mutate(n=1, inc=ifelse(value<10,1,0)) %>%
  group_by(imp.surv,capt.release,Year) %>%
  summarise(ext.prob=sum(inc)/sum(n)) %>%
  mutate(Year=as.numeric(Year)+2018) 


## create factors for plot labels and order them appropriately
extprop$capt.release <- factor(extprop$capt.release, labels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
extprop$imp.surv <- factor(extprop$imp.surv, labels = c("no improvement","surv +2%", "surv +4%", "surv +6%"))



### produce plot with 4 panels and multiple lines per year

#pdf("EV_extinction_probability_C3.pdf", width=10, height=7)
#postscript("Fig1_Balkan.eps", width=9, height=6)
#jpeg("Fig1_Balkan.jpg", width=9, height=6, units="in", res=600, quality=100)
#par(oma=c(0,0,0,0),mar=c(4.2,4.5,0,0.5), cex=1.2)
ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=capt.release), size=1)+
  facet_wrap(~imp.surv,ncol=2) +
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,0.8),breaks=seq(0,0.8,0.2), labels=as.character(seq(0,80,20)))+
  scale_x_continuous(name="Year", breaks=seq(2019,2048,5), labels=as.character(seq(2019,2048,5)))+
  guides(color=guide_legend(title="N captive releases"),fill=guide_legend(title="N captive releases"))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

dev.off()



  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PRODUCE OUTPUT REPORT
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  library(markdown)
  library(rmarkdown)
  library(knitr)
  
  
  ### create HTML report for overall summary report
  Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")
  
  rmarkdown::render('S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\EGVU_CaptRelease.Rmd',
                    output_file = "EGVU_pop_model_captive_breeding.html",
                    output_dir = 'S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel')
  
  
  
  
  
  
  
  


