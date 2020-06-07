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

### 27 APRIL 2020: Included revised and corrected data from 2019 and modified model to allow changes in adult survival
### saved as new file
### included new data on adult survival and changed priors

### 8 May 2020: completely revised IPM to insert model components from adult survival exploration
## switched to a robust-design annual survival CJS model with 5 visits per year (April-Aug)
## juvenile survival adapted from Buechley et al. 2020 telemetry survival model
## REMOVED EGG HARVEST SECTION

### 18 May 2020: huge discussion whether to use only core area from Rhodopes or not -> outsourced into new script
### 27 May 2020: included Albania and N Macedonia and adjusted model to have single constant adult survival

### FINALISED ON 5 JUNE 2020 by including random effect for annual adult survival
### MODEL CONVERGED WITH 30,000 iterations

### OUTPUT NOT SATISFACTORY on 7 JUNE 2020 - reverted survival model from 8 May to insert recruitment into two-phase model
### IPM v4 has 'EGVU_adult_surv_terrmon_CJS_recruit.jags model in it

library(readxl)
library(jagsUI)
library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
library(janitor)
filter<-dplyr::filter
select<-dplyr::select


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE POPULATION MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
head(breed)

trendinput<- occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  group_by(year) %>%
  summarise(N=sum(occupancy), R=sum(breeding, na.rm=T), J=sum(fledglings, na.rm=T))

# countrytrendinput<- occu %>% filter(year>2005) %>%
#   filter(!(Country %in% c("Niger"))) %>%    # introduced in 2020 to remove Niger data
#   group_by(year,Country) %>%
#   summarise(N=sum(occupancy)) %>%
#   spread(key=Country, value=N)

countrytrendinput<- fread("EVGU_countrytrendinput2019.csv")

## calculate mean proportion of each country to total count
country.props<-colMeans(prop.table(as.matrix(countrytrendinput[,-1]), margin = 1), na.rm=T)


breedinput<- breed %>% filter(Year>2005) %>%
  rename(year=Year) %>%
  left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
  filter(!is.na(breed_success)) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  mutate(count=1) %>%
  mutate(fledglings=ifelse(is.na(fledglings),0,fledglings)) %>%
  group_by(year) %>%
  summarise(R=sum(count), J=sum(fledglings))




### PLOT RAW COUNT DATA 

countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  ggplot()+
  geom_line(aes(x=year, y=N, col=Country), size=1)




## export data
#fwrite(trendinput,"EGVU_adult_counts.csv")
#fwrite(breedinput,"EGVU_breeding_summary.csv")


### numbers for manuscript ###
# breed %>% filter(Year>2005) %>%
#   rename(year=Year) %>%
#   left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
#   filter(!(Country %in% c("Niger"))) %>%    # introduced in 2020 to remove Niger data
#   filter(!is.na(breed_success)) %>%
#   group_by(territory_NAME) %>%
#   summarise(Nyears=length(year))

occu %>% filter(year>2005) %>%
  filter(!(Country %in% c("Niger"))) %>%    # introduced in 2020 to remove Niger data
  group_by(Country) %>%
  summarise(N=length(unique(territory_NAME))) %>%
  adorn_totals("row")

breed %>% filter(Year>2005) %>%
  rename(year=Year) %>%
  left_join(occu[,1:4], by=c("territory_NAME","year")) %>%
  filter(!(Country %in% c("Niger"))) %>%    # introduced in 2020 to remove Niger data
  group_by(Country) %>%
  summarise(N=length(unique(territory_NAME))) %>%
  adorn_totals("row")

countrytrendinput %>% adorn_totals("col")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE JUVENILE TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE, intern = T)
load("RODBC_EGVU_telemetry_input.RData")
head(locs)
head(birds)
birds$Fledge_date[is.na(birds$Fledge_date)]<-birds$Tag_date[is.na(birds$Fledge_date)] ### fill in 'fledge' which equals tag date for wild non-juveniles
unique(birds$Age)
birds<-birds %>% filter(Age %in% c("juv","2cal_year")) %>%
  filter(Tag_date<ymd("2019-10-01")) %>%
  filter(!Name=="Zighmund") %>%            ### remove single bird that was never free flying
  filter(!release_method %in% c("hacking","fostering")) %>%            ### remove hacked and fostered birds as we will not use that technique
  dplyr::select(Name,Age,Tag_year,origin,release_method,Status,Fledge_date,Stop_date,Reason_death) %>%
  arrange(Tag_year) %>%
  
  
  # True States (S) - these are often unknown and cannot be observed, we just need them to initialise the model (best guess)
  # 1 dead
  # 2 alive with functioning tag
  # 3 alive with defunct tag OR without tag (when tag was lost)
  
  mutate(TS= ifelse(Status=="Alive",2,1)) %>%
  
  # Observed States (O) - these are based on the actual transmission history
  # 1 Tag ok, bird moving
  # 2 Tag ok, bird not moving (dead)
  # 3 Tag failed, bird observed alive
  # 4 Dead bird recovered
  # 5 No signal (=not seen)
  mutate(OS= ifelse(Status=="Alive",1,
                    ifelse(Status=="Unknown",5,
                           ifelse(Reason_death %in% c("unknown","Unknown"),2,4))))


### CALCULATE DISTANCE BETWEEN SUCCESSIVE LOCATIONS ###
## this approach takes way too long to compute and was removed on 30 April 2020 ##
poss_dist <- possibly(geosphere::distm, otherwise = NA)
locs<- locs %>% filter(Bird_ID %in% birds$Name) %>%
  filter(!LocID %in% c(68850, 35863, 38421, 40832, 3238)) %>%						### manually enter non-sensical GPS locations to be excluded
  mutate(Time=format(Time, format="%H:%M:%S")) %>%
  mutate(DateTime=ymd_hms(paste(Date,Time))) %>%
  mutate(Month=month(DateTime)) %>%
  mutate(PRIMOCC=paste(year(Date),Month, sep="_")) %>%     ## CREATE UNIQUE PRIMARY OCCASION (MONTH WITHIN YEARS)
  arrange(Bird_ID,DateTime) %>%
  select(LocID,Bird_ID,Date,long,lat, PRIMOCC)
head(locs)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A TIME SERIES DATA FRAME ###
mindate<-min(locs$Date)
maxdate<-max(locs$Date)+days(15)
timeseries<-data.frame(date=seq(mindate, maxdate, "1 month")) %>%
  mutate(month=month(date),year=year(date)) %>%
  mutate(date=format(date, format="%m-%Y")) %>%
  mutate(season=ifelse(month %in% c(2,3,4,9,10), 'migration',ifelse(month %in% c(11,12,1),"winter","summer"))) %>%
  mutate(col=seq_along(date)+1)
dim(timeseries)


### CREATE BLANK MATRICES TO HOLD INFORMATION ABOUT TRUE AND OBSERVED STATES ###

EV.obs.matrix<-birds %>% select(Name) %>%
  arrange(Name)
EV.obs.matrix[,2:max(timeseries$col)]<-NA									

EV.state.matrix<-birds %>% select(Name) %>%
  arrange(Name)
EV.state.matrix[,2:max(timeseries$col)]<-NA

EV.age.matrix<-birds %>% select(Name) %>%
  arrange(Name)
EV.age.matrix[,2:max(timeseries$col)]<-NA

EV.mig.matrix<-birds %>% select(Name) %>%
  arrange(Name)
EV.mig.matrix[,2:max(timeseries$col)]<-0

### FILL MATRICES WITH STATE INFORMATION ###
for(n in EV.obs.matrix$Name){
  
  ### extract locations and start dates for each bird
  xl<-locs[locs$Bird_ID==n,]
  mindate<-as.Date(birds$Fledge_date[birds$Name==n])
  mindate<-dplyr::if_else(is.na(mindate),as.Date(min(xl$Date)),mindate)
  if(n %in% c("Elodie", "Odiseas")){mindate<-mindate-15}
  mindate<-format(mindate, format="%m-%Y")
  
  ### extract end dates for each bird
  if(is.na(birds$Stop_date[birds$Name==n])){maxdate<-as.Date(max(xl$Date))}else{			### for birds that are still alive
    maxdate<-as.Date(birds$Stop_date[birds$Name==n])}
  maxdate<-format(maxdate, format="%m-%Y")
  
  ### specify columns in matrix to populate
  startcol<-timeseries$col[timeseries$date==mindate]
  stopcol<-timeseries$col[timeseries$date==maxdate]
  stopcol<-ifelse(stopcol<=startcol,min(startcol+1,max(timeseries$col)),stopcol)
  
  ## ASSIGN OBSERVED STATE
  EV.obs.matrix[EV.obs.matrix$Name==n,startcol:(stopcol-1)]<-1
  if(startcol==stopcol){EV.obs.matrix[EV.obs.matrix$Name==n,2:(stopcol-1)]<-NA} ## for the few cases where stopcol-1 is actually before startcol
  EV.obs.matrix[EV.obs.matrix$Name==n,stopcol:max(timeseries$col)]<-birds$OS[birds$Name==n]   ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!  
  
  ## ASSIGN INITIAL TRUE STATE (to initialise z-matrix of model)
  EV.state.matrix[EV.state.matrix$Name==n,(startcol+1):(stopcol-1)]<-2      ## state at first capture is known, hence must be NA in z-matrix
  EV.state.matrix[EV.state.matrix$Name==n,stopcol:max(timeseries$col)]<-birds$TS[birds$Name==n]    ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!
  EV.state.matrix[EV.state.matrix$Name==n,2:startcol]<-NA ## error occurs if z at first occ is not NA, so we need to specify that for birds alive for <1 month because stopcol-1 = startcol
  
  ## ASSIGN AGES FOR EACH MONTH
  ## find age of bird at start of tracking ##
  xage<-birds$Age[birds$Name==n]
  xage<-ifelse(xage=="juv",1,as.numeric(str_extract_all(xage,"\\(?[0-9]+\\)?", simplify=TRUE)))
  xage<-ifelse(xage==1,1,ifelse(xage==2,9,18))
  agelength<-length(seq(startcol:max(timeseries$col)))
  EV.age.matrix[EV.age.matrix$Name==n,startcol:max(timeseries$col)]<-seq(xage,(xage+agelength-1))      ## insert age progression
  
  ## CALCULATE MIGRATION ACROSS THE SEA FROM LAT DISPLACEMENT FOR FIRST TIME MIGRANTS
  lat.mat<-xl %>% mutate(date=format(Date, format="%m-%Y")) %>%
    left_join(timeseries, by="date") %>% group_by(col) %>% summarise(latdiff=max(lat)-min(lat), latmax=max(lat)) %>%
    mutate(age=seq_along(col)+xage-1) %>%
    mutate(seacross=ifelse(age<18 & latdiff>1.5 & latmax>33,1,0)) %>%
    select(col,seacross) %>%
    spread(key=col,value=seacross, fill=0)
  EV.mig.matrix[EV.age.matrix$Name==n,as.numeric(names(lat.mat))]<-lat.mat      ## insert migration
  if(birds$Reason_death[birds$Name==n] %in% c("drowned","Natural barrier")){
    EV.mig.matrix[EV.age.matrix$Name==n,max(as.numeric(names(lat.mat))):max(timeseries$col)]<-1
  }
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE INPUT DATA FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ENSURE SORT ORDER FOR ALL MATRICES IS IDENTICAL
EV.obs.matrix <- EV.obs.matrix %>% arrange(Name)
EV.state.matrix <- EV.state.matrix %>% arrange(Name)
EV.age.matrix <- EV.age.matrix %>% arrange(Name)
EV.mig.matrix <- EV.mig.matrix %>% arrange(Name)
birds <- birds %>% arrange(Name)


#### Convert to numeric matrices that JAGS can loop over
y.telemetry<-as.matrix(EV.obs.matrix[,2:max(timeseries$col)])
z.telemetry<-as.matrix(EV.state.matrix[,2:max(timeseries$col)])
age.mat<-as.matrix(EV.age.matrix[,2:max(timeseries$col)])
mig.mat<-as.matrix(EV.mig.matrix[,2:max(timeseries$col)])

#### SCALE THE AGE SO THAT NUMERICAL OVERFLOW DOES NOT OCCUR
max(age.mat, na.rm=T)
age.mat<-ifelse(age.mat>48,48,age.mat)
agescale<-scale(1:48)


#### create vector of first marking and of last alive record
get.first.telemetry<-function(x)min(which(!is.na(x)))
get.last.telemetry<-function(x)max(which(!is.na(x) & x==1))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)
l.telemetry<-apply(y.telemetry,1,get.last.telemetry)

rm(locs)

#### simple summaries for manuscript
table(birds$origin)

birds %>% filter(origin=="captive")

birds %>% filter(origin=="wild") %>% filter(Age=="juv") %>%
  mutate(surv=as.numeric(difftime(Stop_date,Fledge_date,"months"))/30) %>%
  mutate(surv=ifelse(surv<8,0,surv)) %>%
  group_by(surv) %>%
  summarise(n=length(Name))

birds %>% filter(origin=="captive") %>%
  mutate(surv=as.numeric(difftime(Stop_date,Fledge_date,"months"))/30) %>%
  mutate(surv=ifelse(surv<8,0,surv)) %>%
  group_by(surv) %>%
  summarise(n=length(Name))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE TERRITORY MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_surveys.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
load("RODBC_EGVU_surveys.RData")
head(EGVU)
dim(EGVU)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMAT DATA AND FIX EFFORT TIMES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(EGVU)

### recalculate effort for visits that have observed >0 birds with 0 effort
updDUR<-EGVU %>% filter(visit_duration==0 | is.na(visit_duration)) %>%
  mutate(dur=as.numeric(difftime(time_end,time_start,units='hours'))) %>%
  mutate(duration=ifelse(dur==0,ifelse(n_adults>0,1,0.1),dur)) %>%
  select(survey_ID,territory_name,date,time_start,time_end,duration)
EGVU$visit_duration[match(updDUR$survey_ID,EGVU$survey_ID)]<-updDUR$duration


### FORMAT AND SELECT RELEVANT COLUMNS
EGVU<- EGVU %>%
  mutate(visit_duration=ifelse(is.na(visit_duration),n_adults,visit_duration)) %>%    ### fill in blanks assuming that effort was >0 when adults were seen
  filter(!is.na(n_adults)) %>%
  filter(visit_duration!=0) %>%
  mutate(n_adults=ifelse(n_adults>2,2,n_adults)) %>%
  select(Country,territory_name,date,visit_duration,n_adults) %>%
  mutate(Year=year(date),JDAY=yday(date),MONTH=month(date)) %>%
  mutate(PRIMOCC=paste(Year,MONTH, sep="_")) %>%     ## CREATE UNIQUE PRIMARY OCCASION (MONTH WITHIN YEARS)
  filter(Year < 2020) %>%
  filter(MONTH %in% c(4:8))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMOVE UNOCCUPIED TERRITORIES ACROSS ALL YEARS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
remove<-EGVU %>% group_by(territory_name) %>% summarise(N=max(n_adults)) %>% filter(N==0)
remove2<-EGVU %>% group_by(territory_name) %>% summarise(first=min(Year)) %>% filter(first>2018)
EGVU<- EGVU %>% filter(!(territory_name %in% remove$territory_name)) %>% filter(!(territory_name %in% remove2$territory_name))
EGVU<- droplevels(EGVU)

dim(EGVU)
head(EGVU)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE MONTHLY VISITS FOR EACH TERRITORY AND YEAR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EGVUsum <- EGVU %>% filter (MONTH %in% c(4,5,6,7,8)) %>%
  group_by(Country,territory_name,Year,MONTH,PRIMOCC) %>%
  summarise(N=max(n_adults, na.rm=T),effort=sum(visit_duration, na.rm=T)) %>%
  ungroup() %>%
  mutate(yearnum=as.numeric(Year)-2005) %>%
  mutate(terrnum=as.numeric(territory_name)) %>%
  arrange(territory_name)

head(EGVUsum)

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])>0) next
    n2 <- max(which(ch[i,]>0))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}



enc.terrvis <- EGVUsum %>%
  group_by(terrnum, Year) %>%
  summarise(occ=max(N,na.rm=T)) %>%
  spread(key=Year, value=occ, fill = 0) %>%
  arrange(terrnum)

param.terrvis <- EGVUsum %>%
  mutate(param=ifelse(Year<2013,1,ifelse(Country=="Greece",1,2))) %>%
  group_by(terrnum, Year) %>%
  summarise(par=mean(param,na.rm=T)) %>%
  spread(key=Year, value=par, fill = 1) %>%
  arrange(terrnum)

obs.terrvis <- EGVUsum %>%
  group_by(terrnum, Year) %>%
  summarise(effort=sum(effort,na.rm=T)) %>%
  spread(key=Year, value=effort, fill = 0) %>%
  arrange(terrnum)

z.terrvis<-as.matrix(enc.terrvis[,2:dim(enc.terrvis)[2]])
z.obsvis<-as.matrix(obs.terrvis[,2:dim(obs.terrvis)[2]])

# Create vector with first occasion of observation
get.first.terrvis <- function(x) min(which(x>0))
f.obsvis <- apply(z.obsvis, 1, get.first.terrvis)
first.obs <- vector()
for (l in 1:dim(z.terrvis)[1]){
  first.obs[l]<-z.terrvis[l,f.obsvis[l]]
}

enc.terrvis[(f.obsvis>12),]
enc.terrvis[(first.obs==0),]
obs.terrvis[(first.obs==0),]
get.first.terrvis(z.obsvis[32,])



get.last.terrvis <- function(x) max(which(x>0))
l.obsvis <- apply(as.matrix(obs.terrvis[32,2:dim(enc.terrvis)[2]]), 1, get.last.terrvis)

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


### PREPARE RANDOM YEAR OFFSET FOR ADULT SURVIVAL

rand.phi.offset<-countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  mutate(prop=ifelse(is.na(N),0,country.props[match(Country,names(country.props))])) %>%
  group_by(year) %>%
  summarise(N=sum(N,na.rm=T),prop=sum(prop)) %>%
  mutate(TOT=N/prop) %>%
  mutate(surv=dplyr::lead(TOT)/TOT) %>%
  mutate(surv=ifelse(is.na(surv),0.999999,surv)) %>%
  mutate(surv=ifelse(surv>0.99,0.999999,surv)) %>%
  mutate(ann.phi=log(surv/(1 - surv)),phi.mean=log(0.93/(1 - 0.93))) %>%
  mutate(rand.phi.offset=ann.phi-phi.mean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PREPARE DATA FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### FOUR ARRAYS NEEDED - FOR OBSERVATIONS AND FOR EFFORT, FOR PARAMETER AND FOR INTERVAL
### EACH ARRAY HAS 3 dimensions for the number of sites (=territories), number of primary occasions and number of visits (=1 per MONTH)
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

# PREPARING THE REQUIRED INPUT DATA:
R.terrvis<-length(unique(EGVUsum$terrnum))		# Number of individuals (= Territory-year combinations)
J.terrvis<-5						# Number of replicate surveys per year, one for each month from April to August
K.terrvis<-length(unique(EGVUsum$Year))	# Number of primary occasions

# Set up some required arrays
site.terrvis <- 1:R.terrvis					# Sites
primoccs.terrvis <- 1:K.terrvis					# primary occasions (months across years)

# SET UP THE OBSERVATION DATA
y.terrvis <- array(NA, dim = c(R.terrvis, J.terrvis, K.terrvis))	# Detection histories
obseff.terrvis <- array(NA, dim = c(R.terrvis, J.terrvis, K.terrvis))	# observation effort


# CREATE A FULL DATA FRAME WITH ONE VALUE PER TERRITORY, PRIMARY AND SECONDARY OCCASION
## needs manual adjustment if number of primary occasions changes
fullEV.terrvis<-data.frame(terrnum=rep(unique(EGVUsum$terrnum),K.terrvis*J.terrvis),
                           Year=rep(unique(EGVUsum$Year),each=R.terrvis*J.terrvis),
                           MONTH=rep(rep(seq(4,8,1),each=R.terrvis), K.terrvis)) %>%
  left_join(EGVUsum,by=c("terrnum", "Year", "MONTH")) %>%
  mutate(yearnum=as.numeric(Year)-2005) %>%
  arrange(terrnum,Year,MONTH)
head(fullEV.terrvis)



# FILL IN ARRAYS FOR OBSERVATIONS AND COVARIATES
for(k in 1:K.terrvis){
  
  dat.terrvis <- fullEV.terrvis %>% filter(yearnum==k) %>%				## filter(WEEK==k+12) 
    select(terrnum, N, MONTH) %>%
    spread(key=MONTH, value=N, fill = NA) %>%
    arrange(terrnum)
  
  effdat.terrvis <- fullEV.terrvis %>% filter(yearnum==k) %>%				## filter(WEEK==k+12) 
    select(terrnum, effort, MONTH) %>%
    spread(key=MONTH, value=effort, fill = 0) %>%
    arrange(terrnum)
  
  y.terrvis[,,k] <-as.matrix(dat.terrvis[,2:6])
  obseff.terrvis[,,k] <-as.matrix(effdat.terrvis[,2:6])
}


# Standardize observation effort covariates
mean.eff.terrvis <- mean(obseff.terrvis, na.rm = TRUE)
sd.eff.terrvis <- sd(obseff.terrvis[!is.na(obseff.terrvis)])
obseff.terrvis <- (obseff.terrvis-mean.eff.terrvis)/sd.eff.terrvis     # Standardise observation effort





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE MATRIX OF RELEASE AND SURVIVAL SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
capt.release=seq(0,15,1)
project.time=c(10,20,30)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.10)
lag.time=c(5,10,15)
PROJECTION.years<-seq(1,50,1)

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
              phase=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2), ### starts in 2006, ends 2019 - phase switch after 2012
              #phase=as.matrix(param.terrvis[,2:dim(param.terrvis)[2]]), ### allow phase to vary between GR and BG - this is horribly complicated and I abandoned it immediately
              rand.phi.offset=rand.phi.offset$rand.phi.offset, ### specification of random annual survival offset on logit scale
              eff.terrvis=obseff.terrvis,
              firstobs=first.obs,
              f.obsvis=f.obsvis, 
              
              y.count=trendinput$N,
              T.count=length(trendinput$N),		## year is standardized so that covariate values are not too far away from zero
              
              R.fec=breedinput$R,
              J.fec=breedinput$J,
              #J.fec.red=breedinput1EGG$J,   ## added to model fecundity when second egg is removed
              
              country.prop=country.props,  ## proportion of count data that come from countries ALB, BG, GR, NMD
              countries=4,
              
              y.telemetry = y.telemetry,
              f.telemetry = f.telemetry,
              l.telemetry = l.telemetry,
              age.telemetry = age.mat, ### matrix(agescale[age.mat], ncol=ncol(age.mat), nrow=nrow(age.mat)), ##scaling age would be a pain for internal transformation
              capt.telemetry = ifelse(birds$origin=="wild",0,1),
              mig.telemetry = mig.mat,
              nind.telemetry = dim(y.telemetry)[1],
              n.occasions.telemetry = dim(y.telemetry)[2],
              migprog.age = c(0,1,rep(0,34)), ## specifies when first migration occurs
              captprog.age =c(0,0,0,0,0,0,0,0,0,0,1,rep(0,25)),
              
              ### Future Projection and SCENARIOS FOR TRAJECTORY
              PROJECTION=30,                ## used 10 and 50 years previously, now trying 30 years
              scen.capt.release=ncol(capt.rel.mat)-1,
              scen.imp.surv=ncol(surv.inc.mat)-1,
              #capt.release=seq(0,15,1),
              #imp.surv=c(1,1.02,1.04,1.06,1.08,1.10))
              capt.release=capt.rel.mat[,2:ncol(capt.rel.mat)],
              imp.surv=surv.inc.mat[,2:ncol(surv.inc.mat)])



## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("mu.fec","lambda.t","b.phi.age","b.phi.capt","b.phi.mig","ann.phi.capt.rel.first.year",
           "ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry", "ann.phi.terrvis",     #"breed.prop4","breed.prop5",
           "mean.phi.terrvis","mean.lambda","fut.lambda","Nterr", "Nterr.f")


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


### CREATE INITIAL VALUES FOR MODEL
INPUT$z.terrvis<-z.obs.terrvis
initIPM <- function(){list(lmu.p.terrvis=runif(dim(z.terrvis)[1],-3, 2),
                           mean.phi.terrvis=runif(2,0.75, 1),
                           sigma.obs.count=runif(1,0,10),
                           mu.fec = runif(1,0,1),
                           z.telemetry = z.telemetry,
                           mean.phi.telemetry = runif(1, 0.9, 0.999), ### two intercepts for juvenile and adults
                           base.obs.telemetry = rnorm(1,0, 0.001),                # Prior for intercept of observation probability on logit scale
                           base.fail.telemetry = rnorm(1,0, 0.001),               # Prior for intercept of tag failure probability on logit scale
                           base.recover.telemetry = rnorm(1,0, 0.001))}   


# ### REDUCE WORKSPACE FOR RUNNING MODEL
# rm.list<-data.frame(object=as.character(ls()), size=0)
# for (obj in ls()) {rm.list$size[rm.list$object==obj]<-object.size(get(obj))}
# rm.list %>% arrange(size)
# rm(list=setdiff(ls(), c("INPUT","initIPM","paraIPM","z.init.terrvis","z.terrvis","yearindex.terrvis","cjs.init.z","capt.rel.mat","surv.inc.mat","z.telemetry","z.obs.terrvis")))
# gc()
# save.image("EGVU_IPM_input_May2020.RData")


# MCMC settings
nc <- 3
nt <- 4
ni <- 5000
nb <- 1000


# ### THIS MODEL QUANTIFIES FUTURE POPULATION TREND FOR A RANGE OF SCENARIOS OF CAPTIVE RELEASES AND SURVIVAL IMPROVEMENT
# NeoIPM.ALL <- autojags(data=INPUT,
#                        inits=initIPM,
#                        parameters.to.save=paraIPM,
#                        model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2020_v1.jags",    ## was EGVU_IPM_2019_COMBINED.jags
#                        n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)##n.iter=ni,
# 
# ## THIS MODEL ASSUMES THAT ALL WILD CHICKS ARE TAKEN INTO TEMPORARY CUSTODY
# NeoIPM.chickremoval <- autojags(data=INPUT,
#                        inits=initIPM,
#                        parameters.to.save=paraIPM,
#                        model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2020_v1_chickremoval.jags",    ## was EGVU_IPM_2019_COMBINED.jags
#                        n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)##n.iter=ni, 
# 
# ## THIS MODEL ASSUMES THAT ALL CAPTIVE-REARED CHICKS ARE TAKEN FROM THE WILD AND CALCULATES HOW MANY NEED TO BE BRED
# paraIPM<-c(paraIPM,"need.to.breed")
# NeoIPM.chicksupplement <- autojags(data=INPUT,
#                                 inits=initIPM,
#                                 parameters.to.save=paraIPM,
#                                 model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2020_v1_chicksupplementation.jags",    ## was EGVU_IPM_2019_COMBINED.jags
#                                 n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)##n.iter=ni, 
# 
# 
# 





### COUNTRY-SPECIFIC MODEL
INPUT$y.count=as.matrix(countrytrendinput[,2:5])
initIPM <- function(){list(lmu.p.terrvis=runif(dim(z.terrvis)[1],-3, 2),
                           mean.phi.terrvis=runif(1,0.75, 1),
                           sigma.obs.count=runif(4,0,10),
                           mu.fec = runif(1,0,1),
                           z.telemetry = z.telemetry,
                           mean.phi.telemetry = runif(1, 0.9, 0.999), ### two intercepts for juvenile and adults
                           base.obs.telemetry = rnorm(1,0, 0.001),                # Prior for intercept of observation probability on logit scale
                           base.fail.telemetry = rnorm(1,0, 0.001),               # Prior for intercept of tag failure probability on logit scale
                           base.recover.telemetry = rnorm(1,0, 0.001))}   



NeoIPM.ALL <- autojags(data=INPUT,
                       inits=initIPM,
                       parameters.to.save=paraIPM,
                       model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2020_v3.jags",    ## was EGVU_IPM_2019_COMBINED.jags
                       n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)##n.iter=ni,


save.image("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\EGVU_papers\\PVA_CaptiveRelease\\EGVU_IPM2020_output_v3_30y.RData")








### MODEL WITH SIMPLE OCCUPANCY (NO ROBUST DESIGN) AND RECRUITMENT
enchist.terrvis<-as.matrix(enc.terrvis[,2:dim(enc.terrvis)[2]])
effort.terrvis<-as.matrix(obs.terrvis[,2:dim(obs.terrvis)[2]])

# Create vector with first occasion of observation
get.first.terrvis <- function(x) min(which(x>0))
f.obsvis <- apply(effort.terrvis, 1, get.first.terrvis)
first.obs <- vector()
for (l in 1:dim(enchist.terrvis)[1]){
  first.obs[l]<-enchist.terrvis[l,f.obsvis[l]]
}

## scale observation effort

range(effort.terrvis)
eff.scale<-scale(0:1993)

INPUT$y.count=as.matrix(countrytrendinput[,2:5])
INPUT$y.terrvis = enchist.terrvis
INPUT$nsite.terrvis = nrow(enchist.terrvis)
INPUT$nprim.terrvis = ncol(enchist.terrvis)
INPUT$phase=c(1,1,1,1,1,1,2,2,2,2,2,2,2,2)  ## cutoff in 2012
INPUT$eff.terrvis=scale(effort.terrvis)
INPUT$firstobs=first.obs
INPUT$f.obsvis=f.obsvis


initIPM <- function(){list(mean.p.terrvis=runif(1,0.5,1),
                           mean.phi.terrvis=runif(2,0.75, 1),
                           sigma.obs.count=runif(4,0,10),
                           mu.fec = runif(1,0,1),
                           z.telemetry = z.telemetry,
                           mean.phi.telemetry = runif(1, 0.9, 0.999), ### two intercepts for juvenile and adults
                           base.obs.telemetry = rnorm(1,0, 0.001),                # Prior for intercept of observation probability on logit scale
                           base.fail.telemetry = rnorm(1,0, 0.001),               # Prior for intercept of tag failure probability on logit scale
                           base.recover.telemetry = rnorm(1,0, 0.001))}   



NeoIPM.ALL <- autojags(data=INPUT,
                       inits=initIPM,
                       parameters.to.save=paraIPM,
                       model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures\\EGVU_IPM_2020_v4.jags",    ## was EGVU_IPM_2019_COMBINED.jags
                       n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)##n.iter=ni,


save.image("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\EGVU_papers\\PVA_CaptiveRelease\\EGVU_IPM2020_output_v4_30y.RData")








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INTEGRATED POPULATION MODEL FOR BALKAN EGYPTIAN VULTURES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2020_v1.jags")
cat("
model {
#-------------------------------------------------
# integrated population model for the balkan population
# - age structured model with 6 age classes: 
# - age-specific probability to recruit at ages 4 or 5
# - age-specific survival derived from tracking data
# - adult survival based on territory occupancy; CHANGES in 2016
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
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)

# Priors and constraints FOR POPULATION COUNTS OBSERVATION
    sigma.obs.count ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
    tau.obs.count<-pow(sigma.obs.count,-2)


# Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
#### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(phi.telemetry[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
                                      b.phi.capt*(capt.telemetry[i]) +     ### survival dependent on captive-release (captive-raised or other)
                                      b.phi.mig*(mig.telemetry[i,t]) +     ### survival dependent on first migration across the sea
                                      b.phi.age*(age.telemetry[i,t])     ### survival dependent on age (juvenile or other)
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi.telemetry ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi.telemetry/(1 - mean.phi.telemetry))    # logit transformed survival intercept
    
    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for captive release on survival probability on logit scale
    
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(p.obs.telemetry[i,t]) <- base.obs.telemetry
        logit(tag.fail.telemetry[i,t]) <- base.fail.telemetry
        logit(p.found.dead.telemetry[i,t]) <- base.recover.telemetry
      } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs.telemetry ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale


# Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    

    ## Priors for detection probability
    # Det prob varies for each site, det prob on logit scale, around 0.05 - 0.95 on normal scale
    for (ipterr in 1:nsite.terrvis){
      lmu.p.terrvis[ipterr] ~ dunif(-3, 2)				
    }							                        ## close k loop for det prob prior
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival
    ## annual survival must be greater than 80% on logit scale
    
    for (nypterr in 1:2){   ## only 2 survival periods
      lm.phi.terrvis[nypterr] <- log(mean.phi.terrvis[nypterr]/(1 - mean.phi.terrvis[nypterr]))    # logit transformed survival intercept		      
      mean.phi.terrvis[nypterr] ~ dunif(0.75, 1)   # informative prior for annual survival probabilities BEFORE 2016
    }    
    
    
    ### RANDOM OBSERVATION EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
      rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    

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

    nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
    N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
    N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
    N4[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N4[tt-1]))                                                ### number of 5-year old survivors
    N6[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds

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
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data



# -------------------------------------------------        
# 2.4. Likelihood for juvenile survival from telemetry
# -------------------------------------------------
 
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind.telemetry){
    
      for (t in f.telemetry[i]:(n.occasions.telemetry-1)){
    
        # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
    
          ps[1,i,t,1]<-1    ## dead birds stay dead
          ps[1,i,t,2]<-0
          ps[1,i,t,3]<-0
    
          ps[2,i,t,1]<-(1-phi.telemetry[i,t])
          ps[2,i,t,2]<-phi.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
          ps[2,i,t,3]<-phi.telemetry[i,t] * tag.fail.telemetry[i,t]
    
          ps[3,i,t,1]<-(1-phi.telemetry[i,t])
          ps[3,i,t,2]<-0
          ps[3,i,t,3]<-phi.telemetry[i,t]
    
        # Define probabilities of O(t) [last dim] given S(t)  [first dim]
    
          po[1,i,t,1]<-0
          po[1,i,t,2]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t]) * (1-p.found.dead.telemetry[i,t])
          po[1,i,t,3]<-0
          po[1,i,t,4]<-p.found.dead.telemetry[i,t]
          po[1,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t] * (1-p.found.dead.telemetry[i,t])
    
          po[2,i,t,1]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
          po[2,i,t,2]<-0
          po[2,i,t,3]<-0
          po[2,i,t,4]<-0
          po[2,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t]
    
          po[3,i,t,1]<-0
          po[3,i,t,2]<-0
          po[3,i,t,3]<-0
          po[3,i,t,4]<-0
          po[3,i,t,5]<-1
    
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind.telemetry){
      # Define latent state at first capture
      z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
        for (t in (f.telemetry[i]+1):n.occasions.telemetry){
          # State process: draw S(t) given S(t-1)
          z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
          # Observation process: draw O(t) given S(t)
          y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
      } #t
    } #i


# -------------------------------------------------        
# 2.5. Likelihood for adult survival from territory monitoring
# -------------------------------------------------

    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
      #z.terrvis[ilterr,f.obsvis[ilterr]] <- firstobs[ilterr]
    
      for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
        z.terrvis[ilterr,klterr] ~ dbin(mean.phi.terrvis[phase[klterr]],z.terrvis[ilterr,klterr-1])
      } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    ## observation prob depends on monthly observation effort
    
    for (iobs in 1:nsite.terrvis){
      for (kobs in f.obsvis[iobs]:nprim.terrvis){
        for (jobs in 1:nrep.terrvis){
    
          linpred.terrvis[iobs,jobs,kobs] <- lmu.p.terrvis[iobs]+beta.obs.eff.terrvis*eff.terrvis[iobs,jobs,kobs] + rand.obs.terrvis[kobs]
          p.terrvis[iobs,jobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,jobs,kobs]))
          y.terrvis[iobs,jobs,kobs] ~ dbin(p.terrvis[iobs,jobs,kobs], z.terrvis[iobs,kobs])
    
    
        } 						# close kobs loop over monthly surveys
      } 						 # close jobs loop over year
    } 							# close iobs loop over sites
    


# -------------------------------------------------        
# 3. DERIVED PARAMETERS
# -------------------------------------------------

### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES

## for WILD BIRDS
for (ageprog in 1:36){
        logit(phi.wild.telemetry[ageprog]) <- lp.mean.telemetry +      ### intercept for mean survival 
                                      b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
                                      b.phi.age*ageprog     ### survival dependent on age (juvenile or other)
}

  ann.phi.juv.telemetry<-phi.wild.telemetry[1]*phi.wild.telemetry[2]*phi.wild.telemetry[3]*phi.wild.telemetry[4]*phi.wild.telemetry[5]*phi.wild.telemetry[6]*phi.wild.telemetry[7]*phi.wild.telemetry[8]*phi.wild.telemetry[9]*phi.wild.telemetry[10]*phi.wild.telemetry[11]*phi.wild.telemetry[12]					### multiply monthly survival from age 1-12
  ann.phi.sec.telemetry<-phi.wild.telemetry[13]*phi.wild.telemetry[14]*phi.wild.telemetry[15]*phi.wild.telemetry[16]*phi.wild.telemetry[17]*phi.wild.telemetry[18]*phi.wild.telemetry[19]*phi.wild.telemetry[20]*phi.wild.telemetry[21]*phi.wild.telemetry[22]*phi.wild.telemetry[23]*phi.wild.telemetry[24]					### multiply monthly survival from age 13-24
  ann.phi.third.telemetry<-phi.wild.telemetry[25]*phi.wild.telemetry[26]*phi.wild.telemetry[27]*phi.wild.telemetry[28]*phi.wild.telemetry[29]*phi.wild.telemetry[30]*phi.wild.telemetry[31]*phi.wild.telemetry[32]*phi.wild.telemetry[33]*phi.wild.telemetry[34]*phi.wild.telemetry[35]*phi.wild.telemetry[36]					### multiply monthly survival from age 25-36

## for CAPTIVE-REARED DELAYED RELEASE BIRDS
for (captageprog in 1:36){
      logit(phi.capt.telemetry[captageprog]) <- lp.mean.telemetry +      ### intercept for mean survival
                                            b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
                                            b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
                                            b.phi.age*captageprog     ### survival dependent on age (juvenile or other)
  }

  ann.phi.capt.rel.first.year<-phi.capt.telemetry[8]*phi.capt.telemetry[9]*phi.capt.telemetry[10]*phi.capt.telemetry[11]*phi.capt.telemetry[12]*phi.capt.telemetry[13]*phi.capt.telemetry[14]*phi.capt.telemetry[15]*phi.capt.telemetry[16]*phi.capt.telemetry[17]*phi.capt.telemetry[18]*phi.capt.telemetry[19]*phi.capt.telemetry[20]*phi.capt.telemetry[21]*phi.capt.telemetry[22]*phi.capt.telemetry[23]*phi.capt.telemetry[24]					### first year for delayed-release bird is longer, from month 8 to 24 

### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    

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

    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){

      ### FUTURE FECUNDITY IS BASED ON mu.fec - removed the chick removal scenario

      for (fut in 1:PROJECTION){
        fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
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
        N2wild.f[ncr,is,1] <- 0
        N2released.f[ncr,is,1] <- 0

          for (fut in 2:PROJECTION){

            fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean.phi.terrvis[2],1) ### invalid parent error if survival>1

            ### probabilistic formulation
            nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
            N1.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### +rescued[ncr,is,fut] number of 1-year old survivors
            N2wild.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))                ### number of 2-year old wild survivors
            N2released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.capt.rel.first.year),1),round(capt.release[fut,ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
            N2.f[ncr,is,fut] <-  N2wild.f[ncr,is,fut] + N2released.f[ncr,is,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
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
# SPECIFY INTEGRATED POPULATION MODEL WITH COUNTRY-SPECIFIC COUNT INPUT AND CONSTANT ADULT SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2020_v2.jags")
cat("
    model {
    #-------------------------------------------------
    # integrated population model for the balkan population
    # - age structured model with 6 age classes: 
    # - age-specific probability to recruit at ages 4 or 5
    # - age-specific survival derived from tracking data
    # - adult survival based on territory occupancy; CHANGES in 2016
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
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
      sigma.obs.count[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    # Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
    #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(phi.telemetry[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
        b.phi.capt*(capt.telemetry[i]) +     ### survival dependent on captive-release (captive-raised or other)
        b.phi.mig*(mig.telemetry[i,t]) +     ### survival dependent on first migration across the sea
        b.phi.age*(age.telemetry[i,t])     ### survival dependent on age (juvenile or other)
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi.telemetry ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi.telemetry/(1 - mean.phi.telemetry))    # logit transformed survival intercept
    
    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for captive release on survival probability on logit scale
    
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(p.obs.telemetry[i,t]) <- base.obs.telemetry
        logit(tag.fail.telemetry[i,t]) <- base.fail.telemetry
        logit(p.found.dead.telemetry[i,t]) <- base.recover.telemetry
      } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs.telemetry ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    
    
    # Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    
    
    ## Priors for detection probability
    # Det prob varies for each site, det prob on logit scale, around 0.05 - 0.95 on normal scale
    for (ipterr in 1:nsite.terrvis){
      lmu.p.terrvis[ipterr] ~ dunif(-3, 2)				
    }							                        ## close k loop for det prob prior
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival
    ## annual survival must be greater than 80% on logit scale
    
    #for (nypterr in 1:2){   ## only 2 survival periods
      lm.phi.terrvis <- log(mean.phi.terrvis/(1 - mean.phi.terrvis))    # logit transformed survival intercept		      
      mean.phi.terrvis ~ dunif(0.75, 1)   # informative prior for annual survival probabilities BEFORE 2016
    #}    
    
    
    ### RANDOM OBSERVATION EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
      rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    
    
    # Priors for population process and immigration
    
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(70, 120)   ##changed from JUV
    N1[1] ~ dunif(15, 45)
    N2[1] ~ dunif(5, 20)
    N3[1] ~ dunif(5, 15)
    N4[1] ~ dunif(4, 14)
    N5[1] ~ dunif(3, 13)
    N6[1] ~ dunif(190, 220)
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T.count){
    
    nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
    N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
    N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
    N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
    N4[tt] ~ dbin(mean.phi.terrvis, round(N3[tt-1]))                                                       ### number of 4-year old survivors
    N5[tt] ~ dbin(mean.phi.terrvis, round(N4[tt-1]))                                                ### number of 5-year old survivors
    N6[tt] ~ dbin(mean.phi.terrvis, round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds
    
    } # tt
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (tlc in 1:T.count){
    Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
      for (s in 1:countries){			### start loop over every country
        Nterr.country[tlc,s] ~ dbin(country.prop[s],round(Nterr[tlc]))
        y.count[tlc,s] ~ dnorm(Nterr[tlc], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods (or binom for single broods after egg harvest)
    # -------------------------------------------------
    
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for juvenile survival from telemetry
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind.telemetry){
    
    for (t in f.telemetry[i]:(n.occasions.telemetry-1)){
    
    # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
    
    ps[1,i,t,1]<-1    ## dead birds stay dead
    ps[1,i,t,2]<-0
    ps[1,i,t,3]<-0
    
    ps[2,i,t,1]<-(1-phi.telemetry[i,t])
    ps[2,i,t,2]<-phi.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
    ps[2,i,t,3]<-phi.telemetry[i,t] * tag.fail.telemetry[i,t]
    
    ps[3,i,t,1]<-(1-phi.telemetry[i,t])
    ps[3,i,t,2]<-0
    ps[3,i,t,3]<-phi.telemetry[i,t]
    
    # Define probabilities of O(t) [last dim] given S(t)  [first dim]
    
    po[1,i,t,1]<-0
    po[1,i,t,2]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t]) * (1-p.found.dead.telemetry[i,t])
    po[1,i,t,3]<-0
    po[1,i,t,4]<-p.found.dead.telemetry[i,t]
    po[1,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t] * (1-p.found.dead.telemetry[i,t])
    
    po[2,i,t,1]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
    po[2,i,t,2]<-0
    po[2,i,t,3]<-0
    po[2,i,t,4]<-0
    po[2,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t]
    
    po[3,i,t,1]<-0
    po[3,i,t,2]<-0
    po[3,i,t,3]<-0
    po[3,i,t,4]<-0
    po[3,i,t,5]<-1
    
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind.telemetry){
    # Define latent state at first capture
    z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
    for (t in (f.telemetry[i]+1):n.occasions.telemetry){
    # State process: draw S(t) given S(t-1)
    z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
    # Observation process: draw O(t) given S(t)
    y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
    } #t
    } #i
    
    
    # -------------------------------------------------        
    # 2.5. Likelihood for adult survival from territory monitoring
    # -------------------------------------------------
    
    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
    #z.terrvis[ilterr,f.obsvis[ilterr]] <- firstobs[ilterr]
    
    for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
    z.terrvis[ilterr,klterr] ~ dbin(mean.phi.terrvis,z.terrvis[ilterr,klterr-1])
    } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    ## observation prob depends on monthly observation effort
    
    for (iobs in 1:nsite.terrvis){
    for (kobs in f.obsvis[iobs]:nprim.terrvis){
    for (jobs in 1:nrep.terrvis){
    
    linpred.terrvis[iobs,jobs,kobs] <- lmu.p.terrvis[iobs]+beta.obs.eff.terrvis*eff.terrvis[iobs,jobs,kobs] + rand.obs.terrvis[kobs]
    p.terrvis[iobs,jobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,jobs,kobs]))
    y.terrvis[iobs,jobs,kobs] ~ dbin(p.terrvis[iobs,jobs,kobs], z.terrvis[iobs,kobs])
    
    
    } 						# close kobs loop over monthly surveys
    } 						 # close jobs loop over year
    } 							# close iobs loop over sites
    
    
    
    # -------------------------------------------------        
    # 3. DERIVED PARAMETERS
    # -------------------------------------------------
    
    ### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES
    
    ## for WILD BIRDS
    for (ageprog in 1:36){
    logit(phi.wild.telemetry[ageprog]) <- lp.mean.telemetry +      ### intercept for mean survival 
    b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
    b.phi.age*ageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.juv.telemetry<-phi.wild.telemetry[1]*phi.wild.telemetry[2]*phi.wild.telemetry[3]*phi.wild.telemetry[4]*phi.wild.telemetry[5]*phi.wild.telemetry[6]*phi.wild.telemetry[7]*phi.wild.telemetry[8]*phi.wild.telemetry[9]*phi.wild.telemetry[10]*phi.wild.telemetry[11]*phi.wild.telemetry[12]					### multiply monthly survival from age 1-12
    ann.phi.sec.telemetry<-phi.wild.telemetry[13]*phi.wild.telemetry[14]*phi.wild.telemetry[15]*phi.wild.telemetry[16]*phi.wild.telemetry[17]*phi.wild.telemetry[18]*phi.wild.telemetry[19]*phi.wild.telemetry[20]*phi.wild.telemetry[21]*phi.wild.telemetry[22]*phi.wild.telemetry[23]*phi.wild.telemetry[24]					### multiply monthly survival from age 13-24
    ann.phi.third.telemetry<-phi.wild.telemetry[25]*phi.wild.telemetry[26]*phi.wild.telemetry[27]*phi.wild.telemetry[28]*phi.wild.telemetry[29]*phi.wild.telemetry[30]*phi.wild.telemetry[31]*phi.wild.telemetry[32]*phi.wild.telemetry[33]*phi.wild.telemetry[34]*phi.wild.telemetry[35]*phi.wild.telemetry[36]					### multiply monthly survival from age 25-36
    
    ## for CAPTIVE-REARED DELAYED RELEASE BIRDS
    for (captageprog in 1:36){
    logit(phi.capt.telemetry[captageprog]) <- lp.mean.telemetry +      ### intercept for mean survival
    b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
    b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
    b.phi.age*captageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.capt.rel.first.year<-phi.capt.telemetry[8]*phi.capt.telemetry[9]*phi.capt.telemetry[10]*phi.capt.telemetry[11]*phi.capt.telemetry[12]*phi.capt.telemetry[13]*phi.capt.telemetry[14]*phi.capt.telemetry[15]*phi.capt.telemetry[16]*phi.capt.telemetry[17]*phi.capt.telemetry[18]*phi.capt.telemetry[19]*phi.capt.telemetry[20]*phi.capt.telemetry[21]*phi.capt.telemetry[22]*phi.capt.telemetry[23]*phi.capt.telemetry[24]					### first year for delayed-release bird is longer, from month 8 to 24 
    
    ### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    
    
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
    
    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){
    
    ### FUTURE FECUNDITY IS BASED ON mu.fec - removed the chick removal scenario
    
    for (fut in 1:PROJECTION){
    fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
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
    N2wild.f[ncr,is,1] <- 0
    N2released.f[ncr,is,1] <- 0
    
    for (fut in 2:PROJECTION){
    
    fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean.phi.terrvis,1) ### invalid parent error if survival>1
    
    ### probabilistic formulation
    nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
    N1.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### +rescued[ncr,is,fut] number of 1-year old survivors
    N2wild.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))                ### number of 2-year old wild survivors
    N2released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.capt.rel.first.year),1),round(capt.release[fut,ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
    N2.f[ncr,is,fut] <-  N2wild.f[ncr,is,fut] + N2released.f[ncr,is,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
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
# SPECIFY INTEGRATED POPULATION MODEL WITH COUNTRY-SPECIFIC COUNT INPUT AND RANDOMLY VARIABLE ADULT SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2020_v3.jags")
cat("
  model {
    #-------------------------------------------------
    # integrated population model for the balkan population
    # - age structured model with 6 age classes: 
    # - age-specific probability to recruit at ages 4 or 5
    # - age-specific survival derived from tracking data
    # - adult survival based on territory occupancy; CHANGES in 2016
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
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
      sigma.obs.count[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    # Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
    #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(phi.telemetry[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
          b.phi.capt*(capt.telemetry[i]) +     ### survival dependent on captive-release (captive-raised or other)
          b.phi.mig*(mig.telemetry[i,t]) +     ### survival dependent on first migration across the sea
          b.phi.age*(age.telemetry[i,t])     ### survival dependent on age (juvenile or other)
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi.telemetry ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi.telemetry/(1 - mean.phi.telemetry))    # logit transformed survival intercept
    
    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for captive release on survival probability on logit scale
    
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(p.obs.telemetry[i,t]) <- base.obs.telemetry
        logit(tag.fail.telemetry[i,t]) <- base.fail.telemetry
        logit(p.found.dead.telemetry[i,t]) <- base.recover.telemetry
      } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs.telemetry ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    
    
    # Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    
    
    ## Priors for detection probability
    # Det prob varies for each site, det prob on logit scale, around 0.05 - 0.95 on normal scale
    for (ipterr in 1:nsite.terrvis){
      lmu.p.terrvis[ipterr] ~ dunif(-3, 2)				
    }							                        ## close k loop for det prob prior
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival
    ## annual survival must be greater than 80% on logit scale
    
    #for (nypterr in 1:2){   ## only 2 survival periods
    lm.phi.terrvis <- log(mean.phi.terrvis/(1 - mean.phi.terrvis))    # logit transformed survival intercept		      
    mean.phi.terrvis ~ dunif(0.75, 1)                                 # informative prior for annual survival probabilities
    #}    
    
    
    ### RANDOM OBSERVATION AND SURVIVAL EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
      rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
      rand.phi.terrvis[nyRpterr] ~ dnorm(rand.phi.offset[nyRpterr], tau.phi.terrvis)  ## include the decline in count data as prior mean for random effect
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)

    tau.phi.terrvis <- 1 / (sd.phi.terrvis * sd.phi.terrvis)
    sd.phi.terrvis ~ dunif(0, 0.05)
    
    
    # Priors for population process and immigration
    
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(70, 120)   ##changed from JUV
    N1[1] ~ dunif(15, 45)
    N2[1] ~ dunif(5, 20)
    N3[1] ~ dunif(5, 15)
    N4[1] ~ dunif(4, 14)
    N5[1] ~ dunif(3, 13)
    N6[1] ~ dunif(200, 250)
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T.count){
      
      nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
      N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
      N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
      N4[tt] ~ dbin(mean.phi.terrvis, round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N5[tt] ~ dbin(mean.phi.terrvis, round(N4[tt-1]))                                                ### number of 5-year old survivors
      N6[tt] ~ dbin(mean.phi.terrvis, round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds
      
    } # tt
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (tlc in 1:T.count){
      Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
      for (s in 1:countries){			### start loop over every country
        Nterr.country[tlc,s] ~ dbin(country.prop[s],round(Nterr[tlc]))
        y.count[tlc,s] ~ dnorm(Nterr[tlc], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }
      
      
      # -------------------------------------------------        
      # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods (or binom for single broods after egg harvest)
      # -------------------------------------------------
      
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for juvenile survival from telemetry
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind.telemetry){
      
      for (t in f.telemetry[i]:(n.occasions.telemetry-1)){
        
        # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
        
        ps[1,i,t,1]<-1    ## dead birds stay dead
        ps[1,i,t,2]<-0
        ps[1,i,t,3]<-0
        
        ps[2,i,t,1]<-(1-phi.telemetry[i,t])
        ps[2,i,t,2]<-phi.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
        ps[2,i,t,3]<-phi.telemetry[i,t] * tag.fail.telemetry[i,t]
        
        ps[3,i,t,1]<-(1-phi.telemetry[i,t])
        ps[3,i,t,2]<-0
        ps[3,i,t,3]<-phi.telemetry[i,t]
        
        # Define probabilities of O(t) [last dim] given S(t)  [first dim]
        
        po[1,i,t,1]<-0
        po[1,i,t,2]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t]) * (1-p.found.dead.telemetry[i,t])
        po[1,i,t,3]<-0
        po[1,i,t,4]<-p.found.dead.telemetry[i,t]
        po[1,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t] * (1-p.found.dead.telemetry[i,t])
        
        po[2,i,t,1]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
        po[2,i,t,2]<-0
        po[2,i,t,3]<-0
        po[2,i,t,4]<-0
        po[2,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t]
        
        po[3,i,t,1]<-0
        po[3,i,t,2]<-0
        po[3,i,t,3]<-0
        po[3,i,t,4]<-0
        po[3,i,t,5]<-1
        
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind.telemetry){
      # Define latent state at first capture
      z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
      for (t in (f.telemetry[i]+1):n.occasions.telemetry){
        # State process: draw S(t) given S(t-1)
        z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
      } #t
    } #i
    
    
    # -------------------------------------------------        
    # 2.5. Likelihood for adult survival from territory monitoring
    # -------------------------------------------------
    
    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
      #z.terrvis[ilterr,f.obsvis[ilterr]] <- firstobs[ilterr]
      
      for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
        z.terrvis[ilterr,klterr] ~ dbin(ann.phi.terrvis[klterr],z.terrvis[ilterr,klterr-1])
      } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites

    ### ANNUAL SURVIVAL PROBABILITY WITH RANDOM YEAR EFFECT
    for (klterr in 1:nprim.terrvis){
        linpred.phi.terrvis[klterr] <- lm.phi.terrvis + rand.phi.terrvis[klterr]
        ann.phi.terrvis[klterr] <- 1 / (1 + exp(-linpred.phi.terrvis[klterr]))
    } 						# close klterr loop over primary period  - years
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    ## observation prob depends on monthly observation effort
    
    for (iobs in 1:nsite.terrvis){
      for (kobs in f.obsvis[iobs]:nprim.terrvis){
        for (jobs in 1:nrep.terrvis){
          
          linpred.terrvis[iobs,jobs,kobs] <- lmu.p.terrvis[iobs]+beta.obs.eff.terrvis*eff.terrvis[iobs,jobs,kobs] + rand.obs.terrvis[kobs]
          p.terrvis[iobs,jobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,jobs,kobs]))
          y.terrvis[iobs,jobs,kobs] ~ dbin(p.terrvis[iobs,jobs,kobs], z.terrvis[iobs,kobs])
          
          
        } 						# close kobs loop over monthly surveys
      } 						 # close jobs loop over year
    } 							# close iobs loop over sites
    
    
    
    # -------------------------------------------------        
    # 3. DERIVED PARAMETERS
    # -------------------------------------------------
    
    ### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES
    
    ## for WILD BIRDS
    for (ageprog in 1:36){
      logit(phi.wild.telemetry[ageprog]) <- lp.mean.telemetry +      ### intercept for mean survival 
        b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*ageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.juv.telemetry<-phi.wild.telemetry[1]*phi.wild.telemetry[2]*phi.wild.telemetry[3]*phi.wild.telemetry[4]*phi.wild.telemetry[5]*phi.wild.telemetry[6]*phi.wild.telemetry[7]*phi.wild.telemetry[8]*phi.wild.telemetry[9]*phi.wild.telemetry[10]*phi.wild.telemetry[11]*phi.wild.telemetry[12]					### multiply monthly survival from age 1-12
    ann.phi.sec.telemetry<-phi.wild.telemetry[13]*phi.wild.telemetry[14]*phi.wild.telemetry[15]*phi.wild.telemetry[16]*phi.wild.telemetry[17]*phi.wild.telemetry[18]*phi.wild.telemetry[19]*phi.wild.telemetry[20]*phi.wild.telemetry[21]*phi.wild.telemetry[22]*phi.wild.telemetry[23]*phi.wild.telemetry[24]					### multiply monthly survival from age 13-24
    ann.phi.third.telemetry<-phi.wild.telemetry[25]*phi.wild.telemetry[26]*phi.wild.telemetry[27]*phi.wild.telemetry[28]*phi.wild.telemetry[29]*phi.wild.telemetry[30]*phi.wild.telemetry[31]*phi.wild.telemetry[32]*phi.wild.telemetry[33]*phi.wild.telemetry[34]*phi.wild.telemetry[35]*phi.wild.telemetry[36]					### multiply monthly survival from age 25-36
    
    ## for CAPTIVE-REARED DELAYED RELEASE BIRDS
    for (captageprog in 1:36){
      logit(phi.capt.telemetry[captageprog]) <- lp.mean.telemetry +      ### intercept for mean survival
        b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
        b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*captageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.capt.rel.first.year<-phi.capt.telemetry[8]*phi.capt.telemetry[9]*phi.capt.telemetry[10]*phi.capt.telemetry[11]*phi.capt.telemetry[12]*phi.capt.telemetry[13]*phi.capt.telemetry[14]*phi.capt.telemetry[15]*phi.capt.telemetry[16]*phi.capt.telemetry[17]*phi.capt.telemetry[18]*phi.capt.telemetry[19]*phi.capt.telemetry[20]*phi.capt.telemetry[21]*phi.capt.telemetry[22]*phi.capt.telemetry[23]*phi.capt.telemetry[24]					### first year for delayed-release bird is longer, from month 8 to 24 
    
    ### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    
    
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
    
    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){
      
      ### FUTURE FECUNDITY IS BASED ON mu.fec - removed the chick removal scenario
      
      for (fut in 1:PROJECTION){
        fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
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
        N2wild.f[ncr,is,1] <- 0
        N2released.f[ncr,is,1] <- 0
        
        for (fut in 2:PROJECTION){
          
          fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean.phi.terrvis,1) ### invalid parent error if survival>1
          
          ### probabilistic formulation
          nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
          N1.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### +rescued[ncr,is,fut] number of 1-year old survivors
          N2wild.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))                ### number of 2-year old wild survivors
          N2released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.capt.rel.first.year),1),round(capt.release[fut,ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
          N2.f[ncr,is,fut] <-  N2wild.f[ncr,is,fut] + N2released.f[ncr,is,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
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
# SPECIFY INTEGRATED POPULATION MODEL WITH COUNTRY-SPECIFIC COUNT INPUT AND TWO PHASE ADULT SURVIVAL CJS MODEL WITH RECRUITMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)

sink("EGVU_IPM_2020_v4.jags")
cat("
  model {
    #-------------------------------------------------
    # integrated population model for the balkan population
    # - age structured model with 6 age classes: 
    # - age-specific probability to recruit at ages 4 or 5
    # - age-specific survival derived from tracking data
    # - adult survival based on territory occupancy; CHANGES in 2016
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
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
      sigma.obs.count[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    # Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
    #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(phi.telemetry[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
          b.phi.capt*(capt.telemetry[i]) +     ### survival dependent on captive-release (captive-raised or other)
          b.phi.mig*(mig.telemetry[i,t]) +     ### survival dependent on first migration across the sea
          b.phi.age*(age.telemetry[i,t])     ### survival dependent on age (juvenile or other)
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi.telemetry ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi.telemetry/(1 - mean.phi.telemetry))    # logit transformed survival intercept
    
    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for captive release on survival probability on logit scale
    
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind.telemetry){
      for (t in f.telemetry[i]:(n.occasions.telemetry)){
        logit(p.obs.telemetry[i,t]) <- base.obs.telemetry
        logit(tag.fail.telemetry[i,t]) <- base.fail.telemetry
        logit(p.found.dead.telemetry[i,t]) <- base.recover.telemetry
      } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs.telemetry ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover.telemetry ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    
    
    # Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    
    
    ## Priors for detection probability
    lmu.p.terrvis <- log(mean.p.terrvis/(1 - mean.p.terrvis))    # logit transformed survival intercept
    mean.p.terrvis ~ dunif(0.5, 1)			
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival
    for (nypterr in 1:2){   ## only 2 survival periods
      mean.phi.terrvis[nypterr] ~ dunif(0.5, 1)   # informative prior for annual survival probabilities BEFORE/AFTER 2016
      mean.rec.terrvis[nypterr] ~ dunif(0, 1)   # informative prior for annual survival probabilities BEFORE/AFTER 2016
    } 
    
    
    ### RANDOM OBSERVATION AND SURVIVAL EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
      rand.obs.terrvis[nyRpterr] ~ dnorm(0, tau.obs.terrvis)
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    
    
    # Priors for population process and immigration
    
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(70, 120)   ##changed from JUV
    N1[1] ~ dunif(15, 45)
    N2[1] ~ dunif(5, 20)
    N3[1] ~ dunif(5, 15)
    N4[1] ~ dunif(4, 14)
    N5[1] ~ dunif(3, 13)
    N6[1] ~ dunif(190, 250)
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T.count){
      
      nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
      N1[tt]  ~ dbin(ann.phi.juv.telemetry, round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
      N2[tt] ~ dbin(ann.phi.sec.telemetry, round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.phi.third.telemetry, round(N2[tt-1]))                                                    ### number of 3-year old survivors
      N4[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N5[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round(N4[tt-1]))                                                ### number of 5-year old survivors
      N6[tt] ~ dbin(mean.phi.terrvis[phase[tt]], round((N5[tt-1]+N6[tt-1])))                                   ### number of 6-year or older (adult) birds
      
    } # tt
    
    
    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (tlc in 1:T.count){
      Nterr[tlc] <- N4[tlc] * 0.024 + N5[tlc] * 0.124 + N6[tlc]                                    ### number of observable territorial birds
      for (s in 1:countries){			### start loop over every country
        Nterr.country[tlc,s] ~ dbin(country.prop[s],round(Nterr[tlc]))
        y.count[tlc,s] ~ dnorm(Nterr[tlc], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }
      
      
      # -------------------------------------------------        
      # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods (or binom for single broods after egg harvest)
      # -------------------------------------------------
      
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for juvenile survival from telemetry
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind.telemetry){
      
      for (t in f.telemetry[i]:(n.occasions.telemetry-1)){
        
        # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
        
        ps[1,i,t,1]<-1    ## dead birds stay dead
        ps[1,i,t,2]<-0
        ps[1,i,t,3]<-0
        
        ps[2,i,t,1]<-(1-phi.telemetry[i,t])
        ps[2,i,t,2]<-phi.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
        ps[2,i,t,3]<-phi.telemetry[i,t] * tag.fail.telemetry[i,t]
        
        ps[3,i,t,1]<-(1-phi.telemetry[i,t])
        ps[3,i,t,2]<-0
        ps[3,i,t,3]<-phi.telemetry[i,t]
        
        # Define probabilities of O(t) [last dim] given S(t)  [first dim]
        
        po[1,i,t,1]<-0
        po[1,i,t,2]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t]) * (1-p.found.dead.telemetry[i,t])
        po[1,i,t,3]<-0
        po[1,i,t,4]<-p.found.dead.telemetry[i,t]
        po[1,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t] * (1-p.found.dead.telemetry[i,t])
        
        po[2,i,t,1]<-p.obs.telemetry[i,t] * (1-tag.fail.telemetry[i,t])
        po[2,i,t,2]<-0
        po[2,i,t,3]<-0
        po[2,i,t,4]<-0
        po[2,i,t,5]<-(1-p.obs.telemetry[i,t]) * tag.fail.telemetry[i,t]
        
        po[3,i,t,1]<-0
        po[3,i,t,2]<-0
        po[3,i,t,3]<-0
        po[3,i,t,4]<-0
        po[3,i,t,5]<-1
        
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind.telemetry){
      # Define latent state at first capture
      z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
      for (t in (f.telemetry[i]+1):n.occasions.telemetry){
        # State process: draw S(t) given S(t-1)
        z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
      } #t
    } #i
    
    
    # -------------------------------------------------        
    # 2.5. Likelihood for adult survival from territory monitoring
    # -------------------------------------------------
    
    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
      x.terrvis[ilterr,f.obsvis[ilterr]] <- 2 #firstobs[ilterr]
      
      for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
        z.terrvis[ilterr,klterr] ~ dbin(mean.phi.terrvis[phase[klterr]],x.terrvis[ilterr,klterr-1])
        mu.x[ilterr,klterr] ~ dbin(mean.rec.terrvis[phase[klterr]],(2-z.terrvis[ilterr,klterr]))
        x.terrvis[ilterr,klterr] <- mu.x[ilterr,klterr] + z.terrvis[ilterr,klterr]
        
      } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR
    
    for (iobs in 1:nsite.terrvis){
      for (kobs in f.obsvis[iobs]:nprim.terrvis){
        linpred.terrvis[iobs,kobs] <- lmu.p.terrvis+beta.obs.eff.terrvis*eff.terrvis[iobs,kobs] + rand.obs.terrvis[kobs]
        p.terrvis[iobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,kobs]))
        y.terrvis[iobs,kobs] ~ dbin(p.terrvis[iobs,kobs], x.terrvis[iobs,kobs])
        
      } 						 # close kobs loop over year
    } 							# close iobs loop over sites
    
    
    
    
    # -------------------------------------------------        
    # 3. DERIVED PARAMETERS
    # -------------------------------------------------
    
    ### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES
    
    ## for WILD BIRDS
    for (ageprog in 1:36){
      logit(phi.wild.telemetry[ageprog]) <- lp.mean.telemetry +      ### intercept for mean survival 
        b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*ageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.juv.telemetry<-phi.wild.telemetry[1]*phi.wild.telemetry[2]*phi.wild.telemetry[3]*phi.wild.telemetry[4]*phi.wild.telemetry[5]*phi.wild.telemetry[6]*phi.wild.telemetry[7]*phi.wild.telemetry[8]*phi.wild.telemetry[9]*phi.wild.telemetry[10]*phi.wild.telemetry[11]*phi.wild.telemetry[12]					### multiply monthly survival from age 1-12
    ann.phi.sec.telemetry<-phi.wild.telemetry[13]*phi.wild.telemetry[14]*phi.wild.telemetry[15]*phi.wild.telemetry[16]*phi.wild.telemetry[17]*phi.wild.telemetry[18]*phi.wild.telemetry[19]*phi.wild.telemetry[20]*phi.wild.telemetry[21]*phi.wild.telemetry[22]*phi.wild.telemetry[23]*phi.wild.telemetry[24]					### multiply monthly survival from age 13-24
    ann.phi.third.telemetry<-phi.wild.telemetry[25]*phi.wild.telemetry[26]*phi.wild.telemetry[27]*phi.wild.telemetry[28]*phi.wild.telemetry[29]*phi.wild.telemetry[30]*phi.wild.telemetry[31]*phi.wild.telemetry[32]*phi.wild.telemetry[33]*phi.wild.telemetry[34]*phi.wild.telemetry[35]*phi.wild.telemetry[36]					### multiply monthly survival from age 25-36
    
    ## for CAPTIVE-REARED DELAYED RELEASE BIRDS
    for (captageprog in 1:36){
      logit(phi.capt.telemetry[captageprog]) <- lp.mean.telemetry +      ### intercept for mean survival
        b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
        b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*captageprog     ### survival dependent on age (juvenile or other)
    }
    
    ann.phi.capt.rel.first.year<-phi.capt.telemetry[8]*phi.capt.telemetry[9]*phi.capt.telemetry[10]*phi.capt.telemetry[11]*phi.capt.telemetry[12]*phi.capt.telemetry[13]*phi.capt.telemetry[14]*phi.capt.telemetry[15]*phi.capt.telemetry[16]*phi.capt.telemetry[17]*phi.capt.telemetry[18]*phi.capt.telemetry[19]*phi.capt.telemetry[20]*phi.capt.telemetry[21]*phi.capt.telemetry[22]*phi.capt.telemetry[23]*phi.capt.telemetry[24]					### first year for delayed-release bird is longer, from month 8 to 24 
    
    ### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    
    
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
    
    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){
      
      ### FUTURE FECUNDITY IS BASED ON mu.fec - removed the chick removal scenario
      
      for (fut in 1:PROJECTION){
        fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
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
        N2wild.f[ncr,is,1] <- 0
        N2released.f[ncr,is,1] <- 0
        
        for (fut in 2:PROJECTION){
          
          fut.survival[ncr,is,fut] <-min(imp.surv[fut,is]*mean.phi.terrvis[2],1) ### invalid parent error if survival>1
          
          ### probabilistic formulation
          nestlings.f[ncr,is,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,is,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
          N1.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.juv.telemetry),1),round(nestlings.f[ncr,is,fut-1]))             ### +rescued[ncr,is,fut] number of 1-year old survivors
          N2wild.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.sec.telemetry),1),round(N1.f[ncr,is,fut-1]))                ### number of 2-year old wild survivors
          N2released.f[ncr,is,fut] ~ dbin(min((imp.surv[fut,is]*ann.phi.capt.rel.first.year),1),round(capt.release[fut,ncr]))             ### +rescued[ncr,is,fut] number of 1-year old survivors 
          N2.f[ncr,is,fut] <-  N2wild.f[ncr,is,fut] + N2released.f[ncr,is,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
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


  