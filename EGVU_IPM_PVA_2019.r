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

## MAJOR UPDATE 18 September 2019
## revised sequence of released chicks

## MAJOR REVSION 13 NOVEMBER 2019
## incorporated suggestions discussed at EVC in Albufeira in October 2019
## increase number of released chicks to 15 per year
## include releases for 10,20, and 30 years
## include gradual improvement of survival (not immediate improvement) over next 5 years

### 22 NOVEMBER 2019: INCLUDED DATA FROM 2019

### 29 NOVEMBER 2019: re-instated chick removal model


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

system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
head(breed)

trendinput<- occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  group_by(year) %>%
  summarise(N=sum(occupancy), R=sum(breeding, na.rm=T), J=sum(fledglings, na.rm=T))


breedinput<- breed %>% filter(Year>2005) %>%
  rename(year=Year) %>%
  left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
  filter(!is.na(breed_success)) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  mutate(count=1) %>%
  mutate(fledglings=ifelse(is.na(fledglings),0,fledglings)) %>%
  group_by(year) %>%
  summarise(R=sum(count), J=sum(fledglings))


### MODIFY BREEDINPUT FOR EGG HARVEST ###

breedinput1EGG<- breed %>% filter(Year>2005) %>%
  rename(year=Year) %>%
  left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  filter(!is.na(breed_success)) %>%
  mutate(count=1) %>%
  mutate(fledglings=ifelse(is.na(fledglings),0,fledglings)) %>%
  mutate(fledglings=ifelse(fledglings==2,1,fledglings)) %>%     ### remove the second egg and fledgling
  group_by(year) %>%
  summarise(R=sum(count), J=sum(fledglings))

breedinput$J-breedinput1EGG$J

### numbers for manuscript
breed %>% filter(Year>2005) %>%
  rename(year=Year) %>%
  left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  filter(!is.na(breed_success)) %>%
  group_by(territory_NAME) %>%
  summarise(Nyears=length(year))

occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%
  group_by(Country) %>%
  summarise(N=length(unique(territory_NAME)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE JUVENILE TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE)
load("RODBC_EGVU_telemetry_input.RData")
head(birds)

birds$Fledge_date[is.na(birds$Fledge_date)]<-birds$Tag_date[is.na(birds$Fledge_date)] ### fill in 'fledge' which equals tag date for wild non-juveniles


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
  filter(Age %in% c("juv","2cal_year","3cal_year","4cal_year")) %>%                   ### changed in 2019 to include immatures caught in Ethiopia
  select(Name,Tag_year, Age,origin) %>% 
  filter(!Name=="Zighmund") %>%            ### remove single bird that was never free flying
  filter(origin=="wild")     %>%               ### USE ONLY WILD JUVENILES
  arrange(Tag_year)

CH.telemetry[,5:64]<-0									### NEEDS MANUAL ADJUSTMENT IF REPEATED FOR LONGER TIME SERIES
x.telemetry<-as.matrix(CH.telemetry[,5:64]) ### create matrix with age progression

dim(CH.telemetry)
range(CH.telemetry$Tag_year)

### FILL CAPTURE HISTORY WITH LIVE DEAD INFORMATION ###
### THIS MAY REQUIRE A FIX IN THE FUTURE IF BIRDS THAT ARE STILL ALIVE HAVE NOT BEEN TAGGED > 4 YEARS AGO
### WOULD NEED FILLING WITH NA


for(n in CH.telemetry$Name){
  
  ### extract locations and start and end dates for each bird
  xl<-locs[locs$Bird_ID==n,]
  mindate<-as.Date(birds$Fledge_date[birds$Name==n])
  mindate<-dplyr::if_else(is.na(mindate),as.Date(min(xl$Date)),mindate)
  if(n %in% c("Elodie", "Odiseas")){mindate<-mindate-15}
  if(is.na(birds$Stop_date[birds$Name==n])){maxdate<-as.Date(max(xl$Date))}else{			### for birds that are still alive
  maxdate<-as.Date(birds$Stop_date[birds$Name==n])}
  birdseries<-data.frame(date=seq(mindate, maxdate, "1 month"), live=1)
  stopcol<-length(birdseries$live)+4
  stopcol<-ifelse(stopcol>64,64,stopcol)
  CH.telemetry[CH.telemetry$Name==n,5:stopcol]<-birdseries$live[1:(stopcol-4)]
  
  ### adjust end for birds not tagged>5 years ago
  if(stopcol<64){
    if(birds$Status[birds$Name==n]=="Alive"){
      CH.telemetry[CH.telemetry$Name==n,(stopcol+1):64]<-NA    ### set to NA because these occasions are in the future
    }
  } 
  
  ## find age of bird ##
  xage<-CH.telemetry$Age[CH.telemetry$Name==n]
  xage<-ifelse(xage=="juv",1,as.numeric(str_extract_all(xage,"\\(?[0-9]+\\)?", simplify=TRUE)))

  
  ## create matrix of age progression
  ## MIGRATION ENDS IN OCTOBER, so only 3 occ (Aug, Sept, Oct) in high risk category - extended to 5 months in 2019
  ## survival assumed equal for first 5 months, then another 12 months, then until end
  agechange1<-ifelse(xage==1,5,12)   ## define when bird switches to next age category
  agechange2<-ifelse(xage==1,17,24)  ## define when bird switches to last (third age category)
  
  x.telemetry[CH.telemetry$Name==n,1:agechange1]<-as.numeric(min(xage,3))			## Set everything after the first occasion to age at marking, capped at 3
  x.telemetry[CH.telemetry$Name==n,(agechange1+1):agechange2]<-as.numeric(min(xage+1,3))		## change to next age category, capped at 3
  x.telemetry[CH.telemetry$Name==n,(agechange2+1):60]<-as.numeric(min(xage+2,3))				## change to third and last age category, capped at 3

}

y.telemetry<-as.matrix(CH.telemetry[,5:64])
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
head(EGVU)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM MACEDONIA AND MERGE WITH BG DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### OBSOLETE AS OF 15 NOV because data imported into project database

# #EGVUm<-read_excel("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Raw_Data\\archived_databases\\ReturnNeophron\\C4_C5\\Macedonia_mgmt.xlsx", sheet = "Survival")
# EGVUm<-read_excel("C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\archived_databases\\ReturnNeophron\\C4_C5\\Macedonia_mgmt.xlsx", sheet = "Survival")
# 
# head(EGVUm)
# EGVUm$visit_duration<-as.numeric(difftime(EGVUm$time_end,EGVUm$time_start, units = "hours"))
# EGVUm$visit_duration[is.na(EGVUm$visit_duration)]<-mean(EGVUm$visit_duration,na.rm=T)
# 
# ### filter the territories that have very few visits
# EGVUm$Year<-as.numeric(format(EGVUm$date, format="%Y"))
# EGVUm$terryear<-paste(EGVUm$Year,EGVUm$territory_name, sep="_")
# mgmtMAC<-aggregate(count~Country+territory_name+Year+terryear, data=EGVUm, FUN=sum)
# shit<-subset(mgmtMAC, count<2)
# 
# ### remove the worthless surveys from the main table
# EGVUm<-EGVUm[!(EGVUm$terryear %in% shit$terryear),] %>%
#   mutate(STAGE=NA, STATUS=NA,failure_reason=NA,n_chicks=0,n_other=0,nest_ID=NA,GlobalID=NA) %>%  ## insert the columns needed to join data
#   select(names(EGVU))
# 
# EGVU<-rbind(EGVU,EGVUm)
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


### report sample size in manuscript
dim(y.terrvis)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INTEGRATED POPULATION MODEL FOR BALKAN EGYPTIAN VULTURES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2019_ExtendedProjection.jags")
cat("
model {
#-------------------------------------------------
# integrated population model for the balkan population
# - age structured model with 6 age classes: 
# - age-specific probability to recruit at ages 4 or 5
# - age-specific survival derived from tracking data
# - adult survival based on territory occupancy
# - pre breeding census, female-based assuming equal sex ratio & survival
# - FUTURE PROJECTION COMBINES SCENARIOS: baseline scenario (=do nothing, no chick removal)
# - productivity supplemented by captive-bred juveniles (0-10 per year)
# - fecundity reduced by harvesting second eggs ONLY FOR SCENARIOS OF CAPTIVE RELEASES
# - rescue of wild chicks optional (currently disabled)
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
### INCLUDE THE RESCUE OF 9 chicks and release with increased survival



    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){

      ### FUTURE FECUNDITY IS BASED ON mu.fec[1] for no action scenario and mu.fec[2] when there is captive-released birds

      for (fut in 1:5){
        fut.fec[fut,ncr] <-mu.fec[min(capt.release[fut,ncr]+1,2)]  ## this will be mu.fec[1] for 0 capt.release and mu.fec[2] when captive birds are released for first 5 years all second chicks are taken into captivity
      }
      for (fut in 6:PROJECTION){
        fut.fec[fut,ncr] <-mu.fec[1]  ## no chicks taken after 5 years anymore
      }



      # SPECIFY IMPROVEMENT OF SURVIVAL
      for (is in 1:scen.imp.surv){ 
    

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
        #rescued[ncr,is,1] <- 0                ## we need to fill in a value to make the matrix complete, this is not actually used in any calculation
        N1nestlings.f[ncr,is,1] <- 0
        N1released.f[ncr,is,1] <- N1[T.count]

          for (fut in 2:PROJECTION){

            fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean(ann.surv.terrvis[1:nyears.terrvis]),1) ### invalid parent error if survival>1

            ### probabilistic formulation
            #rescued[ncr,is,fut] ~ dpois(5) T(1,11)                                         ### number of chicks rescued and rehabilitated to improve survival FROM HARVESTING SECOND EGGS
            nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
            N1nestlings.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
            N1released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(capt.release[fut,ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
            N1.f[ncr,is,fut] <-  N1nestlings.f[ncr,is,fut] + N1released.f[ncr,is,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
            N2.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))  ### number of 2-year old survivors
            N3.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.third.telemetry),1),round(N2.f[ncr,is,fut-1]))                                                    ### number of 3-year old survivors
            N4.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round(N3.f[ncr,is,fut-1]))                                                       ### number of 4-year old survivors
            N5.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round(N4.f[ncr,is,fut-1]))                                                       ### number of 5-year old survivors
            N6.f[ncr,is,fut] ~ dbin(fut.survival[ncr,is,fut],round((N5.f[ncr,is,fut-1]+N6.f[ncr,is,fut-1])))                                   ### number of 6-year or older (adult) birds
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
# CREATE MATRIX OF RELEASE AND SURVIVAL SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
capt.release=seq(0,15,1)
project.time=c(10,20,30)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.10)
lag.time=c(5,10,15)
PROJECTION.years<-seq(1,30,1)

capt.rel.mat<- expand.grid(PROJECTION.years,capt.release,project.time) %>%
  rename(Year=Var1,REL=Var2,DUR=Var3) %>%
  mutate(scenario=paste(REL,DUR, sep="_")) %>%
  mutate(REL=ifelse(Year>DUR,0,REL)) %>%
  select(-DUR) %>%
  spread(key=scenario, value=REL)


surv.inc.mat<- expand.grid(PROJECTION.years,imp.surv,lag.time) %>%
  rename(Year=Var1,SURV=Var2,LAG=Var3) %>%
  mutate(scenario=paste(SURV,LAG, sep="_")) %>%
  mutate(ann.offset=(SURV-1)/LAG) %>%
  mutate(SURV=ifelse(Year<LAG,1+(Year*ann.offset),SURV)) %>%
  select(-LAG,-ann.offset) %>%
  spread(key=scenario, value=SURV)




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
              PROJECTION=30,                ## used 10 and 50 years previously, now trying 30 years
              scen.capt.release=ncol(capt.rel.mat)-1,
              scen.imp.surv=ncol(surv.inc.mat)-1,
              #capt.release=seq(0,15,1),
              #imp.surv=c(1,1.02,1.04,1.06,1.08,1.10))
              capt.release=capt.rel.mat[,2:ncol(capt.rel.mat)],
              imp.surv=surv.inc.mat[,2:ncol(surv.inc.mat)])



## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("mu.fec","lambda.t","ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry",      #"breed.prop4","breed.prop5",
            "ann.surv.terrvis","mean.p.terrvis","mean.lambda","fut.lambda","Nterr", "Nterr.f")
          


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


### THIS MODEL SIMPLY EXPLORES THE NO CHANGE SCENARIO AND QUANTIFIES FUTURE POPULATION TREND IF NOTHING IS DONE
# NeoIPMbasic <- jags(data=INPUT,
#                 inits=initIPM,
#                 parameters.to.save=paraIPM,
#                 model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_Baseline.jags",
#                 n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

# ### THIS MODEL QUANTIFIES FUTURE POPULATION TREND IF FOR 5 YEARS ALL SECOND EGGS ARE REMOVED AND NO CHICKS ARE RELEASED  
NeoIPMeggredNoRescue <- jags(data=INPUT,
                       inits=initIPM,
                       parameters.to.save=paraIPM,
                       model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_EggRem5_NoRescue.jags",
                       n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

# ### THIS MODEL QUANTIFIES FUTURE POPULATION TREND IF FOR 5 YEARS ALL SECOND EGGS ARE REMOVED AND CHICKS ARE RELEASED EVERY YEAR
# NeoIPMeggredRescue <- jags(data=INPUT,
#                     inits=initIPM,
#                     parameters.to.save=paraIPM,
#                     model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_EggRem5_WithRescue.jags",
#                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

### THIS MODEL QUANTIFIES FUTURE POPULATION TREND FOR A RANGE OF SCENARIOS OF CAPTIVE RELEASES AND SURVIVAL IMPROVEMENT
NeoIPM.ALL <- jags(data=INPUT,
                   inits=initIPM,
                   parameters.to.save=paraIPM,
                   model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_ExtendedProjection.jags",    ## was EGVU_IPM_2019_COMBINED.jags
                   n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)

# ### THIS MODEL QUANTIFIES FUTURE POPULATION TREND FOR A RANGE OF SCENARIOS OF SURVIVAL IMPROVEMENT WITH CAPTIVE RELEASES ONLY IN THE FIRST 5-10 YEARS
# NeoIPM.RED <- jags(data=INPUT,
#                            inits=initIPM,
#                            parameters.to.save=paraIPM,
#                            model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2019_REDRELEASE.jags",
#                            n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
  

save.image("EGVU_IPM_output2019_v4.RData")



  


  
  
  
  


