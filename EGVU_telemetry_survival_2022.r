##########################################################################
#
# EGYPTIAN VULTURE MONTHLY SURVIVAL ANALYSIS FROM TELEMETRY - BALKANS ONLY
#
##########################################################################
# written by Steffen Oppel, July 2019
# data preparation by Evan Buechley, Ron Efrat, Louis Phipps, and Evan Buechley
# branched from "EGVU_telemetry_survival_analysis.r' on 9 Jan 2020

## MAJOR CHANGES ON 7 APRIL 2020
## removed all residents
## removed age effect and simply retained 2 age groups (<18 months and >18 months)
## test 3 models with breed,mig, and winter stage and west, Balkan and Caucasus populations, either additive or interactive

## REDUCED TO BALKANS BIRDS ONLY ON 14 APRIL 2020
## REQUIRES MAJOR MODEL ADJUSTMENT TO REMOVE NON-EXISTING OBSERVATION STATES (no bird with failed transmitter was ever resighted in Balkans)

## UPDATED 28 APRIL 2020 to include data import from database
## UPDATED 24 June 2020 to exclude Ethiopian birds and check difference in results

## REVISED IN JANUARY 2023 FOR FINAL REPORT AND MANUSCRIPT
## REVISED 6 MARCH 2023 to include longitudinal offset for survival

## REVISED 22 JUNE 2023 to include GoF test
## check utility of this: https://cran.r-project.org/web/packages/ggmcmc/vignettes/using_ggmcmc.html

## MAJOR REVISIONS ON 11 JULY TO MODIFY ANNUAL SURVIVAL ESTIMATES
## these depend on the value for longitude that is specified

library(jagsUI)
library(runjags)
library(tidyverse)
library(data.table)
library(lubridate)
library(tidyverse)
library(geosphere)
filter<-dplyr::filter
select<-dplyr::select
library(MCMCvis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD DATA FROM DATABASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### No new data 
# try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
# try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
# #system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE)
# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_telemetry_input.r")), wait = TRUE, invisible = FALSE)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Data"), silent=T)
# load("RODBC_EGVU_telemetry_input.RData")
# #load("EGVU_telemetry_survival_formatted_data.RData")
# head(locs)
# head(birds)
# birds$Fledge_date[is.na(birds$Fledge_date)]<-birds$Tag_date[is.na(birds$Fledge_date)] ### fill in 'fledge' which equals tag date for wild non-juveniles
# unique(birds$Age)
# birds<-birds %>% filter(Age %in% c("juv","2cal_year","3cal_year")) %>%
#   filter(Tag_date<ymd("2022-10-01")) %>%
#   filter(!Name=="Zighmund") %>%            ### remove single bird that was never free flying
#   #filter(!release_method %in% c("hacking","fostering")) %>%            ### remove hacked and fostered birds as we will not use that technique
#   #filter(origin=="wild")     %>%               ### USE ONLY WILD JUVENILES
#   filter(!(Region %in% c("Oromia","Afar","Syria")))  %>%               ### USE ONLY BALKAN BIRDS
#   dplyr::select(Name,Age,Tag_year,origin,release_method,Status,Fledge_date,Stop_date,Reason_death,Tag_long) %>%
#   arrange(Tag_year) %>%
#   
#   
#   # True States (S) - these are often unknown and cannot be observed, we just need them to initialise the model (best guess)
#   # 1 dead
#   # 2 alive with functioning tag
#   # 3 alive with defunct tag OR without tag (when tag was lost)
#   
#   mutate(TS= ifelse(Status=="Alive",2,1)) %>%
# 
# # Observed States (O) - these are based on the actual transmission history
# # 1 Tag ok, bird moving
# # 2 Tag ok, bird not moving (dead)
# # 3 Tag failed, bird observed alive
# # 4 Dead bird recovered
# # 5 No signal (=not seen)
#   mutate(OS= ifelse(Status=="Alive",1,
#                     ifelse(Status=="Unknown",5,
#                       ifelse(Reason_death %in% c("unknown","Unknown"),2,4))))
# 
# 
# ### CALCULATE DISTANCE BETWEEN SUCCESSIVE LOCATIONS ###
# ## this approach takes way too long to compute and was removed on 30 April 2020 ##
# ## reinstated in 2023 to have a more objective definition of when migration occurs
# #poss_dist <- possibly(geosphere::distm, otherwise = NA)
# orig_locs<-locs
# locs<- orig_locs %>% filter(Bird_ID %in% birds$Name) %>%
#   filter(!LocID %in% c(68850, 35863, 38421, 40832, 3238)) %>%						### manually enter non-sensical GPS locations to be excluded
#   mutate(Time=format(Time, format="%H:%M:%S")) %>%
#   mutate(Date=ymd_hms(paste(as.Date(Date),Time))) %>%
#   mutate(Month=month(Date)) %>%
#   mutate(PRIMOCC=paste(year(Date),Month, sep="_")) %>%     ## CREATE UNIQUE PRIMARY OCCASION (MONTH WITHIN YEARS)
#   mutate(date=format(Date, format="%m-%Y")) %>%
#   dplyr::select(-Time,-loc_type,-loc_quality,-temperature,-battery_voltage,-TxCount,-SatelliteCount,-duplicate) %>%
#   #nest(long, lat, .key = "coords") %>%
#   arrange(Bird_ID,Date) %>%
#   group_by(Bird_ID) %>%
#   #mutate(prev_coords = lag(coords)) %>%
#   mutate(distance = c(NA, geosphere::distVincentyEllipsoid(cbind(long, lat)))) %>%
#   #mutate(distance = map2_dbl(coords, prev_coords, poss_dist)) %>%
#   #unnest(coords) %>%
#   ungroup() %>%
#   select(LocID,Bird_ID,Date,date,long,lat, PRIMOCC, distance)
# head(locs)
# str(locs)
# 
# summary(locs$distance)
# locs %>% filter(distance>1000000)
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ### CREATE A TIME SERIES DATA FRAME ###
# mindate<-min(locs$Date)
# maxdate<-ymd_hms("2023-03-01 12:00:00")   ### curtail data collection to February 2023 ### max(locs$Date)+days(15)
# timeseries<-data.frame(date=seq(mindate, maxdate, "1 month")) %>%
#   mutate(month=month(date),year=year(date)) %>%
#   mutate(LIFE=ifelse(year(date)>2017,1,0)) %>%
#   mutate(date=format(date, format="%m-%Y")) %>%
#   mutate(season=ifelse(month %in% c(2,3,4,9,10), 'migration',ifelse(month %in% c(11,12,1),"winter","summer"))) %>%
#   mutate(col=seq_along(date)+1)
# dim(timeseries)
# 
# 
# ### CREATE BLANK MATRICES TO HOLD INFORMATION ABOUT TRUE AND OBSERVED STATES ###
# 
# EV.obs.matrix<-birds %>% select(Name) %>%
#   arrange(Name)
# EV.obs.matrix[,2:max(timeseries$col)]<-NA									
# 
# EV.state.matrix<-birds %>% select(Name) %>%
#   arrange(Name)
# EV.state.matrix[,2:max(timeseries$col)]<-NA
# 
# EV.age.matrix<-birds %>% select(Name) %>%
#   arrange(Name)
# EV.age.matrix[,2:max(timeseries$col)]<-NA
# 
# EV.mig.matrix<-birds %>% select(Name) %>%
#   arrange(Name)
# EV.mig.matrix[,2:max(timeseries$col)]<-0
# 
# EV.LIFE.matrix<-birds %>% select(Name) %>%
#   arrange(Name)
# EV.LIFE.matrix[,2:max(timeseries$col)]<-0
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # SUMMARISE MIGRATION STATE FOR BIRDS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ## define migration as >4000 km per month or 5 degree lat or long displacement
# 
# EV.mig.matrix.raw <- locs %>% filter(distance<1000000) %>%
#   filter(Date<as.Date(ymd("2023-03-01"))) %>%    ### curtail data collection to February 2023
#   group_by(Bird_ID,PRIMOCC, date) %>%
#   summarise(tot_dist=sum(distance, na.rm=T), lat_range=max(lat)-min(lat), long_range=max(long)-min(long)) %>%
#   mutate(mig=ifelse(tot_dist>4000000,1,ifelse(lat_range>5 | long_range>5,1,0))) %>%
#   left_join(timeseries, by='date')
# 
# hist(EV.mig.matrix.raw$tot_dist/1000, breaks=100)
# hist(EV.mig.matrix.raw$lat_range, breaks=30)
# hist(EV.mig.matrix.raw$long_range, breaks=30)
# 
# 
# ### THIS MATRIX IS INCOMPLETE BECAUSE MANY MONTHS DO NOT HAVE ANY LOCATIONS
# #EV.mig.matrix2 <- EV.mig.matrix.raw %>% dplyr::select(-tot_dist,-lat_range,-long_range,-LIFE,-season,-PRIMOCC,-month,-year,-date) %>%
# #  spread(key=col, value=mig, fill=0)
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # SUMMARISE SAMPLE SIZES FOR MANUSCRIPT
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# head(birds)
# table(birds$origin)
# table(birds$Tag_year)
# table(birds$release_method)
# table(birds$Status)
# table(birds$Reason_death[birds$Status=="Dead"])
# 16/34 ## natural cause of death
# 11/34 ## unknown cause
# 7/34 ## human cause of death
# 
# fwrite(birds,"EGVU_juv_survival_estimation.csv")
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # CREATE CAPTURE HISTORY FOR SURVIVAL ESTIMATIONS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ### CREATE A TIME SERIES DATA FRAME ###
# # mindate<-min(EV$start)
# # maxdate<-max(EV$end)
# # timeseries<-data.frame(date=seq(mindate, maxdate, "1 month")) %>%
# #   mutate(month=month(date),year=year(date)) %>%
# #   mutate(date=format(date, format="%m-%Y")) %>%
# #   mutate(season=ifelse(month %in% c(2,3,4,9,10), 'migration',ifelse(month %in% c(11,12,1),"winter","summer"))) %>%
# #   mutate(col=seq_along(date)+1)
# # dim(timeseries)
# 
# 
# ### FILL MATRICES WITH STATE INFORMATION ###
# for(n in EV.obs.matrix$Name){
#   
#   ### extract locations and start dates for each bird
#   xl<-locs[locs$Bird_ID==n,] %>% filter(Date<as.Date(ymd("2023-03-01")))   ### curtail data collection to February 2023
#   mindate<-as.Date(birds$Fledge_date[birds$Name==n])
#   mindate<-dplyr::if_else(is.na(mindate),as.Date(min(xl$Date)),mindate)
#   if(n %in% c("Elodie", "Odiseas")){mindate<-mindate-15}
#   mindate<-format(mindate, format="%m-%Y")
#   
#   ### extract end dates for each bird
#   if(is.na(birds$Stop_date[birds$Name==n])){maxdate<-as.Date(max(xl$Date))}else{			### for birds that are still alive
#     maxdate<-as.Date(birds$Stop_date[birds$Name==n])}
#   maxdate<-format(maxdate, format="%m-%Y")
# 
#   ### specify columns in matrix to populate
#   startcol<-timeseries$col[timeseries$date==mindate]
#   stopcol<-timeseries$col[timeseries$date==maxdate]
#   stopcol<-ifelse(stopcol<=startcol,min(startcol+1,max(timeseries$col)),stopcol)
# 
#   ## ASSIGN OBSERVED STATE
#   EV.obs.matrix[EV.obs.matrix$Name==n,startcol:(stopcol-1)]<-1
#   if(startcol==stopcol){EV.obs.matrix[EV.obs.matrix$Name==n,2:(stopcol-1)]<-NA} ## for the few cases where stopcol-1 is actually before startcol
#   EV.obs.matrix[EV.obs.matrix$Name==n,stopcol:max(timeseries$col)]<-birds$OS[birds$Name==n]   ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!  
#   
#   ## ASSIGN LIFE PROJECT STATE
#   EV.LIFE.matrix[EV.LIFE.matrix$Name==n,startcol:(stopcol-1)]<-timeseries$LIFE[(startcol-1):(stopcol-2)]
#   
#   ## ASSIGN INITIAL TRUE STATE (to initialise z-matrix of model)
#   EV.state.matrix[EV.state.matrix$Name==n,(startcol+1):(stopcol-1)]<-2      ## state at first capture is known, hence must be NA in z-matrix
#   EV.state.matrix[EV.state.matrix$Name==n,stopcol:max(timeseries$col)]<-birds$TS[birds$Name==n]    ## this assumes that state never changes - some birds may have gone off air and then been found dead - this needs manual adjustment!
#   EV.state.matrix[EV.state.matrix$Name==n,2:startcol]<-NA ## error occurs if z at first occ is not NA, so we need to specify that for birds alive for <1 month because stopcol-1 = startcol
#   
#   ## ASSIGN AGES FOR EACH MONTH
#   ## find age of bird at start of tracking ##
#   xage<-birds$Age[birds$Name==n]
#   xage<-ifelse(xage=="juv",1,as.numeric(str_extract_all(xage,"\\(?[0-9]+\\)?", simplify=TRUE)))
#   xage<-ifelse(xage==1,1,ifelse(xage==2,9,18))
#   agelength<-length(seq(startcol:max(timeseries$col)))
#   EV.age.matrix[EV.age.matrix$Name==n,startcol:max(timeseries$col)]<-seq(xage,(xage+agelength-1))      ## insert age progression
#   
#   # ## CALCULATE MIGRATION ACROSS THE SEA FROM LAT DISPLACEMENT FOR FIRST TIME MIGRANTS
#   # lat.mat<-xl %>% mutate(date=format(Date, format="%m-%Y")) %>%
#   #   left_join(timeseries, by="date") %>% group_by(col) %>% summarise(latdiff=max(lat)-min(lat), latmax=max(lat)) %>%
#   #   mutate(age=seq_along(col)+xage-1) %>%
#   #   mutate(seacross=ifelse(age<18 & latdiff>1.5 & latmax>33,1,0)) %>%
#   #   select(col,seacross) %>%
#   #   spread(key=col,value=seacross, fill=0)
#   # EV.mig.matrix[EV.age.matrix$Name==n,as.numeric(names(lat.mat))]<-lat.mat      ## insert migration
#   # if(birds$Reason_death[birds$Name==n] %in% c("drowned","Natural barrier")){
#   #   EV.mig.matrix[EV.age.matrix$Name==n,max(as.numeric(names(lat.mat))):max(timeseries$col)]<-1
#   # }
#   
#   ## USE AN ALTERNATIVE DEFINITION OF MIGRATION BASED ON DISTANCE AND DISPLACEMENT
#   mig.mat<-EV.mig.matrix.raw %>% filter(Bird_ID==n) %>% ungroup() %>%
#     dplyr::select(col,mig) %>%
#     spread(key=col,value=mig, fill=0)
#   EV.mig.matrix[EV.age.matrix$Name==n,as.numeric(names(mig.mat))]<-mig.mat      ## insert migration
#   if(birds$Reason_death[birds$Name==n] %in% c("drowned","Natural barrier")){
#     EV.mig.matrix[EV.age.matrix$Name==n,max(as.numeric(names(mig.mat))):max(timeseries$col)]<-1
#   }
# 
# }
# 
# 
# save.image("EGVU_telemetry_survival_formatted_data.RData")
load("EGVU_telemetry_survival_formatted_data.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE INPUT DATA FOR JAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### ENSURE SORT ORDER FOR ALL MATRICES IS IDENTICAL
EV.obs.matrix <- EV.obs.matrix %>% arrange(Name)
EV.state.matrix <- EV.state.matrix %>% arrange(Name)
EV.age.matrix <- EV.age.matrix %>% arrange(Name)
EV.mig.matrix <- EV.mig.matrix %>% arrange(Name)
EV.LIFE.matrix <- EV.LIFE.matrix %>% arrange(Name)
birds <- birds %>% arrange(Name)


#### Convert to numeric matrices that JAGS can loop over
y.telemetry<-as.matrix(EV.obs.matrix[,2:max(timeseries$col)])
z.telemetry<-as.matrix(EV.state.matrix[,2:max(timeseries$col)])
age.mat<-as.matrix(EV.age.matrix[,2:max(timeseries$col)])
mig.mat<-as.matrix(EV.mig.matrix[,2:max(timeseries$col)])
LIFE.mat<-as.matrix(EV.LIFE.matrix[,2:max(timeseries$col)])

#### SCALE THE AGE SO THAT NUMERICAL OVERFLOW DOES NOT OCCUR
max(age.mat, na.rm=T)
age.mat<-ifelse(age.mat>48,48,age.mat)
agescale<-scale(1:48)


#### create vector of first marking and of last alive record
get.first.telemetry<-function(x)min(which(!is.na(x)))
get.last.telemetry<-function(x)max(which(!is.na(x) & x==1))
f.telemetry<-apply(y.telemetry,1,get.first.telemetry)
l.telemetry<-apply(y.telemetry,1,get.last.telemetry)

## create "juv" matrix for first autumn migration, which differs for delayed release birds
juv.mat<-ifelse(age.mat<12,1,0)
birds$release_method<-ifelse(is.na(birds$release_method),"none",birds$release_method)
juv.mat[birds$release_method=="delay" & birds$Age=="2cal_year",]<-ifelse(age.mat[birds$release_method=="delay" & birds$Age=="2cal_year",]<18,1,0)
juv.mat[birds$release_method=="delay" & birds$Age=="3cal_year",]<-ifelse(age.mat[birds$release_method=="delay" & birds$Age=="3cal_year",]<26,1,0)


#### BUNDLE DATA INTO A LIST
INPUT.telemetry <- list(y = y.telemetry,
                        f = f.telemetry,
                        l = l.telemetry,
                        age = matrix(agescale[age.mat], ncol=ncol(age.mat), nrow=nrow(age.mat)),
                        capt = ifelse(birds$origin=="wild",0,1),
                        LIFE=LIFE.mat,
                        mig = mig.mat,
                        long = as.numeric(scale(birds$Tag_long)),
                        juv =juv.mat,   ## first migration for captive birds is later
                        nind = dim(y.telemetry)[1],
                        n.occasions = dim(y.telemetry)[2])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY AND SET UP MODEL RUNS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters monitored
parameters.telemetry <- c("base.obs","base.fail","base.recover","mean.phi","lp.mean","b.phi.age","b.phi.capt","b.phi.mig","b.phi.LIFE","b.phi.long")

# Initial values for some parameters
inits.telemetry <- function(){list(z = z.telemetry,
                                   mean.phi = runif(1, 0.9, 0.999), ### two intercepts for juvenile and adults
					                         base.obs = rnorm(1,0, 0.001),                # Prior for intercept of observation probability on logit scale
						                      base.fail = rnorm(1,0, 0.001),               # Prior for intercept of tag failure probability on logit scale
						                      base.recover = rnorm(1,0, 0.001))} 

# MCMC settings
ni <- 500
nt <- 4
nb <- 200
nc <- 3
nad<-100
ns<-1000  	## number of iterations (draws per chain)



### MIGRATION ONLY MODELS WITH NO GEOGRAPHIC STRUCTURE ############################
# Call JAGS from R (took 4 hrs on PC, DIC = 1817.493)
EGVU_juv_surv <- autojags(INPUT.telemetry, inits.telemetry, parameters.telemetry,
                                         #"C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\EGVU_telemetry_juv_survival_LIFE.jags",
                                         "C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_telemetry_juv_survival_LIFE_long.jags",
                                         n.chains = nc, n.thin = nt, n.burnin = nb, n.cores=nc, parallel=T) #, n.iter = ni)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
EGVU_juv_surv <- run.jags(data=INPUT.telemetry,
                       inits=inits.telemetry,
                       monitor=parameters.telemetry,
                       model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_telemetry_juv_survival_LIFE_long.jags",
                       #model="EGVU_telemetry_juv_survival_LIFE_long.jags",    ### C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\
                       n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                       method = "rjparallel")  



save.image("EGVU_LIFE_telemetry_survival_output_2023.RData")
load("EGVU_LIFE_telemetry_survival_output_2023.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ MADE UP GOF TEST  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### suggested by Pius Korner in July 2023
### simulate data from model
### compare frequency of states from simulated prediction to frequency in observed data

OBS<-table(as.factor(y.telemetry))
REPraw<-MCMCpstr(EGVU_juv_surv$mcmc, params=c("rep.states"), type="chains")
REP<-table(as.factor(as.numeric(REPraw$rep.states)))
chisq.test(OBS,REP)  ### should have a p-value >> 0.05 otherwise there would be disconcerting lack of fit



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ PRODUCE PARMETER ESTIMATE PLOT  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out<-as.data.frame(EGVU_juv_surv$summary) #%>% bind_cols(as.data.frame(EGVU_juv_surv$summary[[2]]))
out$parameter<-row.names(EGVU_juv_surv$summary)
#out<-as.data.frame(summary(EGVU_juv_surv))
#out$parameter<-row.names(summary(EGVU_juv_surv))

#names(out)[c(3,5,7)]<-c('lcl','median', 'ucl')
names(out)[c(1,2,3)]<-c('lcl','median', 'ucl')
out<-out %>%  select(parameter,Mean, median, lcl, ucl)

fwrite(out,"EGVU_LIFE_telemetry_surv_parm_est.csv")
#out<-fread("EGVU_LIFE_telemetry_surv_parm_est.csv")

### PLOT PARAMETERS ON LOGIT SCALE
plotdat<-out %>% filter(grepl("b.phi",parameter)) %>%
  mutate(parameter=ifelse(parameter=="b.phi.mig","migration",parameter)) %>%
  mutate(parameter=ifelse(parameter=="b.phi.age","age",parameter)) %>%
  mutate(parameter=ifelse(parameter=="b.phi.capt","origin",parameter)) %>%
  mutate(parameter=ifelse(parameter=="b.phi.LIFE","LIFE",parameter)) %>%
  mutate(parameter=ifelse(parameter=="b.phi.long","longitude",parameter))

ggplot(plotdat)+
  geom_point(aes(x=parameter, y=median))+
  geom_point(data=plotdat[plotdat$parameter=="LIFE",],aes(x=parameter, y=median), size=2, colour="red")+
  geom_errorbar(aes(x=parameter, ymin=lcl, ymax=ucl), width=.1) +
  geom_point(data=plotdat[plotdat$parameter=="LIFE",],aes(x=parameter, y=median), size=2, colour="red")+
  geom_hline(aes(yintercept=0), colour="darkgrey") +
  
  ## format axis ticks
  xlab("Monthly survival parameter") +
  ylab("estimate (logit scale)") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("EGVU_parameter_estimates_LIFE.jpg", height=7, width=10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ PRODUCE MONTHY SURVIVAL OUTPUT GRAPH SHOWING LIFE IMPACT  #############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### PREPARE RAW MCMC OUTPUT
parmcols<-dimnames(EGVU_juv_surv$mcmc[[1]])[[2]]

# back-scale longitude
longscale<-scale(birds$Tag_long)

# ### COMBINE SAMPLES ACROSS CHAINS
MCMCout<-rbind(EGVU_juv_surv$mcmc[[1]],EGVU_juv_surv$mcmc[[2]],EGVU_juv_surv$mcmc[[3]])
# str(MCMCout)

### SET UP ANNUAL TABLE FOR PLOTTING THE MONTHLY SURVIVAL GRAPH
## does not include migration

AnnTab<-expand.grid(capt=0,age=seq(1,36), LIFE=c(0,1), long = 25.5)
AnnTab<-expand.grid(capt=1,age=seq(9,23), LIFE=c(0,1), long = 25.5) %>%     ## for delayed release juveniles
  bind_rows(AnnTab) %>%
  mutate(mig=0) %>% ### to ensure there is not a weird notch in the graph
  mutate(juv=ifelse(capt==0,ifelse(age <13,1,0),
                    ifelse(age <18,1,0))) %>% 
  mutate(scaleage=(age-attr(agescale, 'scaled:center'))/attr(agescale, 'scaled:scale')) %>%
  mutate(scalelong=(long-attr(longscale, 'scaled:center'))/attr(longscale, 'scaled:scale'))%>%
  mutate(AgeGroup=ifelse(age>24,"third",ifelse(age<13 & capt==0,"juv","second")))

Xin<-AnnTab

### CALCULATE PREDICTED VALUE FOR EACH SAMPLE

MCMCpred<-data.frame()
for(s in 1:nrow(MCMCout)) {

  X<-  Xin %>%
    
    ##CALCULATE MONTHLY SURVIVAL
    mutate(logit.surv=as.numeric(MCMCout[s,match("lp.mean",parmcols)])+
             as.numeric(MCMCout[s,match("b.phi.age",parmcols)])*scaleage +
             as.numeric(MCMCout[s,match("b.phi.mig",parmcols)])*mig +
             as.numeric(MCMCout[s,match("b.phi.LIFE",parmcols)])*LIFE +
             as.numeric(MCMCout[s,match("b.phi.capt",parmcols)])*capt +
             as.numeric(MCMCout[s,match("b.phi.long",parmcols)])*scalelong*juv) %>%
    
    ##BACKTRANSFORM TO NORMAL SCALE
    mutate(surv=plogis(logit.surv)) %>%
    
    ##CALCULATE ANNUAL SURVIVAL
    mutate(Origin=ifelse(capt==0,"wild","captive-reared")) %>%
    mutate(LIFE=ifelse(LIFE==0,"before LIFE","with LIFE")) %>%
    mutate(simul=s)              
  
  
  MCMCpred<-rbind(MCMCpred,as.data.frame(X)) 
  
}

### CREATE PLOT

MCMCpred %>%   rename(raw.surv=surv) %>%
  group_by(LIFE,Origin, age) %>%
  summarise(surv=quantile(raw.surv,0.5),surv.lcl=quantile(raw.surv,0.025),surv.ucl=quantile(raw.surv,0.975)) %>%
  mutate(LIFE=if_else(LIFE=="before LIFE", "before flyway project","since flyway project started")) %>%
	filter(Origin=="wild") %>% ##filter(LIFE=="before LIFE")
  ggplot()+
  geom_ribbon(aes(x=age, ymin=surv.lcl, ymax=surv.ucl, fill=LIFE), alpha=0.2) +   ##, type=Origin
  geom_line(aes(x=age, y=surv, color=LIFE),linewidth=2)+     ## , linetype=Origin

  ## format axis ticks
  scale_x_continuous(name="Age in months", limits=c(1,12), breaks=seq(1,12,1), labels=seq(1,12,1)) +
  #scale_y_continuous(name="Monthly survival probability", limits=c(0.8,1), breaks=seq(0.,1,0.05)) +
  labs(y="Monthly survival probability",fill="Period", colour="Period") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +

  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=16, color="black"),
        legend.position=c(0.82,0.18), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))


ggsave("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\MANUSCRIPTS\\in_prep\\EGVU_LIFE\\FIG_3.jpg", height=7, width=11)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE ANNUAL SURVIVAL FOR PVA PROJECTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



### SET UP ANNUAL TABLE FOR ESTIMATING ANNUAL SURVIVAL
### includes migration cost

AnnTab<-expand.grid(capt=0,age=seq(1,36), LIFE=c(0,1), long = 25.5)
AnnTab<-expand.grid(capt=1,age=seq(9,23), LIFE=c(0,1), long = 25.5) %>%     ## for delayed release juveniles
  bind_rows(AnnTab) %>%
  # mutate(mig=ifelse(capt==0,ifelse(age %in% c(4,21,27,32),1,0),
  #                   ifelse(age %in% c(12,21),1,0))) %>%
  mutate(mig=ifelse(capt==0,ifelse(age %in% c(3,4,
                                              20,21,22,26,27,28,32,33,34),1,0),
                    ifelse(age %in% c(14,15,16,20,21,22),1,0))) %>% 
  mutate(juv=ifelse(capt==0,ifelse(age <13,1,0),
                    ifelse(age <18,1,0))) %>% 
  mutate(scaleage=(age-attr(agescale, 'scaled:center'))/attr(agescale, 'scaled:scale')) %>%
  mutate(scalelong=(long-attr(longscale, 'scaled:center'))/attr(longscale, 'scaled:scale'))%>%
  mutate(AgeGroup=ifelse(age>24,"third",ifelse(age<13 & capt==0,"juv","second")))

Xin<-AnnTab

### CALCULATE PREDICTED VALUE FOR EACH SAMPLE
# takes too long, hence reduced to use summary output instead
# but using summary values of parameters creates weird confidence intervals for before LIFE period

MCMCpred<-data.frame()
for(s in 1:nrow(MCMCout)) {
  
  X<-  Xin %>%
    
    ##CALCULATE MONTHLY SURVIVAL
    mutate(logit.surv=as.numeric(MCMCout[s,match("lp.mean",parmcols)])+
             as.numeric(MCMCout[s,match("b.phi.age",parmcols)])*scaleage +
             as.numeric(MCMCout[s,match("b.phi.mig",parmcols)])*mig +
             as.numeric(MCMCout[s,match("b.phi.LIFE",parmcols)])*LIFE +
             as.numeric(MCMCout[s,match("b.phi.capt",parmcols)])*capt +
             as.numeric(MCMCout[s,match("b.phi.long",parmcols)])*scalelong*juv) %>%
    
    ##BACKTRANSFORM TO NORMAL SCALE
    mutate(surv=plogis(logit.surv)) %>%
    
    ##CALCULATE ANNUAL SURVIVAL
    mutate(Origin=ifelse(capt==0,"wild","captive-reared")) %>%
    mutate(LIFE=ifelse(LIFE==0,"before LIFE","with LIFE")) %>%
    mutate(simul=s)              
  
  
  MCMCpred<-rbind(MCMCpred,as.data.frame(X)) 
  
}


### CALCULATE PREDICTED SURVIVAL BASED ON FINAL MODEL

surv.parm<-  MCMCpred %>% 
  
  ### CALCULATE ANNUAL SURVIVAL
  group_by(Origin,LIFE,simul,AgeGroup) %>%
  summarise(ann.surv=prod(surv)) %>%
  ungroup() %>%
  
  ### CALCULATE CREDIBLE INTERVALS
  group_by(LIFE,Origin,AgeGroup) %>%
  summarise(mean.surv=mean(ann.surv),sd=sd(ann.surv),med.surv=quantile(ann.surv,0.5),lcl.surv=quantile(ann.surv,0.025),ucl.surv=quantile(ann.surv,0.975)) %>%
  arrange(LIFE,Origin)

surv.parm
fwrite(surv.parm,"LIFE_EGVU_juvenile_surv_estimates.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE PROPORTIONAL INCREASE IN SUVIVAL FOR FINAL REPORT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCMCpred %>% group_by(LIFE,AgeGroup,Origin) %>%
  summarise(med.surv=quantile(surv,0.5)) %>%
  spread(key=LIFE, value=med.surv) %>%
  mutate(prop.increase=((`with LIFE`-`before LIFE`)/`before LIFE`) *100)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT / SAVE INTERMEDIATE OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival")
try(setwd("C:\\Users\\steffenoppel\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
save.image("EGVU_LIFE_telemetry_survival_output_2023.RData")
#load("EGVU_LIFE_telemetry_survival_output_2023.RData")









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SETTING UP JAGS CODE FOR LIFE PROJECT MODEL (needs to be run first)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Specify model in BUGS language
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival")
setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival")
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
sink("EGVU_telemetry_juv_survival_LIFE_long.jags")
cat("
  model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: monthly survival probability intercept
    # tag.fail: probability that tag will fail
    
    # p.obs: probability to be tracked with functioning tag (=1)
    # p.found.dead: probability for carcass to be recovered
    # p.seen.alive: probability to be observed alive despite the tag being defunct
    
    # -------------------------------------------------
    # States (S):
    # 1 dead
    # 2 alive with functioning tag
    # 3 alive with defunct tag or tag lost
    
    # Observations (O):
    # 1 Tag ok, bird moving
    # 2 Tag ok, bird not moving (dead, or tag lost and no longer on bird)
    # 3 Tag failed, bird observed alive
    # 4 Dead bird recovered
    # 5 No signal (=not seen)
    
    # -------------------------------------------------
    
    # Priors and constraints
    
    
    #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind){
      for (t in f[i]:(n.occasions)){
        logit(phi[i,t]) <- lp.mean +      ### intercept for mean survival 
            b.phi.capt*(capt[i]) +     ### survival dependent on captive-release (captive-raised or other)
            b.phi.mig*(mig[i,t]) +     ### survival dependent on migration or not (in this month)
            b.phi.age*(age[i,t]) +    ### survival dependent on age (juvenile or other)
            b.phi.LIFE*(LIFE[i,t]) +     ### survival dependent on LIFE project
            b.phi.long*(long[i])*(juv[i,t]) ### survival varies for juvenile first-time migrants when they originate from farther west
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean <- log(mean.phi/(1 - mean.phi))    # logit transformed survival intercept

    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dnorm(0, 0.001)           # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for captive release on survival probability on logit scale
    b.phi.LIFE ~ dnorm(0, 0.001)         # Prior for LIFE effect on survival probability on logit scale    
    b.phi.long ~ dnorm(0, 0.001)         # Prior for longitudinal effect on survival probability on logit scale  
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind){
      for (t in f[i]:(n.occasions)){
        logit(p.obs[i,t]) <- base.obs
        logit(tag.fail[i,t]) <- base.fail #+ beta3*tfail[i] #### probability of TAG FAILURE is influenced by tag type and tag age
        logit(p.found.dead[i,t]) <- base.recover
      } #t
    } #i
    
    
    ##### SLOPE PARAMETERS FOR OBSERVATION PROBABILITY
    base.obs ~ dnorm(0, 0.001)                # Prior for intercept of observation probability on logit scale
    base.fail ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    base.recover ~ dnorm(0, 0.001)               # Prior for intercept of tag failure probability on logit scale
    #beta3 ~ dnorm(0, 0.001)T(-10, 10)         # Prior for slope parameter for fail probability with tag movement during last 10 GPS fixes
    
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind){
      
      for (t in f[i]:(n.occasions-1)){
        
        # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
        
        ps[1,i,t,1]<-1    ## dead birds stay dead
        ps[1,i,t,2]<-0
        ps[1,i,t,3]<-0
        
        ps[2,i,t,1]<-(1-phi[i,t])
        ps[2,i,t,2]<-phi[i,t] * (1-tag.fail[i,t])
        ps[2,i,t,3]<-phi[i,t] * tag.fail[i,t]
        
        ps[3,i,t,1]<-(1-phi[i,t])
        ps[3,i,t,2]<-0
        ps[3,i,t,3]<-phi[i,t]
        
        # Define probabilities of O(t) [last dim] given S(t)  [first dim]
        
        po[1,i,t,1]<-0
        po[1,i,t,2]<-p.obs[i,t] * (1-tag.fail[i,t]) * (1-p.found.dead[i,t])
        po[1,i,t,3]<-0
        po[1,i,t,4]<-p.found.dead[i,t]
        po[1,i,t,5]<-(1-p.obs[i,t]) * tag.fail[i,t] * (1-p.found.dead[i,t])
        
        po[2,i,t,1]<-p.obs[i,t] * (1-tag.fail[i,t])
        po[2,i,t,2]<-0
        po[2,i,t,3]<-0
        po[2,i,t,4]<-0
        po[2,i,t,5]<-(1-p.obs[i,t]) * tag.fail[i,t]
        
        po[3,i,t,1]<-0
        po[3,i,t,2]<-0
        po[3,i,t,3]<-0
        po[3,i,t,4]<-0
        po[3,i,t,5]<-1
        
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 2 ## alive when first marked
      for (t in (f[i]+1):n.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1,])
        rep.states[i,t] ~ dcat(po[z[i,t], i, t-1,])    ## include this as goodness of fit measure
      } #t
    } #i
  }
    ",fill = TRUE)
sink()



