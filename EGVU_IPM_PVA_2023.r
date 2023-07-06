##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS
#
##########################################################################
# based on Lieury et al. 2015
# published by Steffen Oppel in 2021 (Journal of Applied Ecology)
# goal is to assess population viability with and without captive released juveniles on Balkans

# update in June 2023 to address reviewer comment

### REVISED IN JUNE 2023 to include GOF TEST FOR SSM and implement the changes from Animal Conservation reviewer
### removed future scenarios to just have single future projection

library(readxl)
library(jagsUI)
library(runjags)
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
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Data"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
head(breed)
countrytrendinput<- fread("EVGU_countrytrendinput2019.csv") %>%
  filter(year<2018)

countrytrendinput<- occu %>% filter(year>2017) %>%
  filter(!(Country=="Niger")) %>%    # introduced in 2019 because database now has data from Niger
  group_by(year,Country) %>%
  summarise(N=sum(occupancy)) %>%
  spread(key=Country,value=N) %>%
  bind_rows(countrytrendinput) %>%
  arrange(year)
countrytrendinput


## calculate mean proportion of each country to total count
country.props<-colMeans(prop.table(as.matrix(countrytrendinput[,-1]), margin = 1), na.rm=T)


## PRODUCTIVITY INPUT
trendinput<- occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  group_by(year) %>%
  summarise(N=sum(occupancy), R=sum(breeding, na.rm=T), J=sum(fledglings, na.rm=T))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE JUVENILE TRACKING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
rm(locs, orig.locs)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE TERRITORY MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Data"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
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
  select(territory_name,date,visit_duration,n_adults) %>%
  mutate(Year=year(date),JDAY=yday(date),MONTH=month(date)) %>%
  mutate(PRIMOCC=paste(Year,MONTH, sep="_")) %>%     ## CREATE UNIQUE PRIMARY OCCASION (MONTH WITHIN YEARS)
  #filter(Year < 2020) %>%
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
  group_by(territory_name,Year,MONTH,PRIMOCC) %>%
  summarise(N=max(n_adults, na.rm=T),effort=sum(visit_duration, na.rm=T)) %>%
  ungroup() %>%
  mutate(yearnum=as.numeric(Year)-2005) %>%
  mutate(terrnum=as.numeric(as.factor(territory_name))) %>%
  arrange(territory_name)

head(EGVUsum)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP SIMPLE ENCOUNTER HISTORY AND EFFORT MATRIX
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

enc.terrvis <- EGVUsum %>%
  group_by(terrnum, Year) %>%
  summarise(occ=max(N,na.rm=T)) %>%
  spread(key=Year, value=occ, fill = 0) %>%
  arrange(terrnum)

obs.terrvis <- EGVUsum %>%
  group_by(terrnum, Year) %>%
  summarise(effort=sum(effort,na.rm=T)) %>%
  spread(key=Year, value=effort, fill = 0) %>%
  arrange(terrnum)

enchist.terrvis<-as.matrix(enc.terrvis[,2:dim(enc.terrvis)[2]])
effort.terrvis<-as.matrix(obs.terrvis[,2:dim(obs.terrvis)[2]])

# Create vector with first occasion of observation
get.first.terrvis <- function(x) min(which(x>0))
f.obsvis <- apply(effort.terrvis, 1, get.first.terrvis)
first.obs <- vector()
for (l in 1:dim(enchist.terrvis)[1]){
  first.obs[l]<-enchist.terrvis[l,f.obsvis[l]]
}


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


## scale observation effort

range(effort.terrvis)
eff.scale<-scale(0:1993)


# PREPARING THE REQUIRED INPUT DATA:
R.terrvis<-length(unique(EGVUsum$terrnum))		# Number of individuals (= Territory-year combinations)
K.terrvis<-length(unique(EGVUsum$Year))	# Number of primary occasions



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE MATRIX OF RELEASE AND SURVIVAL SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
capt.release=5
project.time=20
imp.surv=c(1)
PROJECTION.years<-seq(1,50,1)

capt.rel.mat<- expand.grid(PROJECTION.years,capt.release,project.time) %>%
  rename(Year=Var1,REL=Var2,DUR=Var3) %>%
  mutate(scenario=paste(REL,DUR, sep="_")) %>%
  mutate(REL=ifelse(Year>DUR,0,REL)) %>%
  select(-DUR) %>%
  spread(key=scenario, value=REL)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DATA AND INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Bundle data across all of the data sources

INPUT <- list(### ADULT SURVIVAL MODEL
              y.terrvis = enchist.terrvis,
              nsite.terrvis = nrow(enchist.terrvis),
              nprim.terrvis = ncol(enchist.terrvis),
              phase=c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),  ## test of before and during project
              eff.terrvis=scale(effort.terrvis),
              f.obsvis=f.obsvis,
              
              ### STATE SPACE MODEL FOR TREND
              y.count=as.matrix(countrytrendinput[,2:5]),
              T.count=dim(countrytrendinput)[1],
              country.prop=as.numeric(country.props),  ## proportion of count data that come from countries ALB, BG, GR, NMD
              countries=4,
              period=ifelse(trendinput$year<2018,1,2),
              
              ### PRODUCTIVITY MODEL
              R.fec=trendinput$R,
              J.fec=trendinput$J,

              ### MULTISTATE MODEL FOR JUVENILE SURVIVAL
              y.telemetry = y.telemetry,
              f.telemetry = f.telemetry,
              age = matrix(agescale[age.mat], ncol=ncol(age.mat), nrow=nrow(age.mat)),
              capt = ifelse(birds$origin=="wild",0,1),
              LIFE=LIFE.mat,
              mig = mig.mat,
              long = as.numeric(scale(birds$Tag_long)),
              juv =juv.mat,   ## first migration for captive birds is later
              nind = dim(y.telemetry)[1],
              n.occasions = dim(y.telemetry)[2],
              
              ### Future Projection and SCENARIOS FOR TRAJECTORY
              PROJECTION=30,                ## used 10 and 50 years previously, now trying 30 years
              scen.capt.release=ncol(capt.rel.mat)-1,
              captprog.age=c(rep(0,10),1,rep(0,25)),
              migprog.age=c(0,1,rep(0,34)),
              capt.release=capt.rel.mat[,2:ncol(capt.rel.mat)])



## Parameters to be estimated ('monitored') by JAGS
paraIPM<-c("mu.fec","lambda.t","b.phi.age","b.phi.capt","b.phi.mig","ann.phi.capt.rel.first.year",
           "ann.phi.juv.telemetry","ann.phi.sec.telemetry","ann.phi.third.telemetry", "ann.phi.terrvis", "base.fail.telemetry",    #"breed.prop4","breed.prop5",
           "mean.phi.terrvis","baseline.lambda","LIFE.lambda","fut.lambda","Nterr", "Nterr.f",
           "Dmape.rep","Dmape.obs","fit.rep.terr","fit.terr")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


initIPM <- function(){list(z.telemetry = z.telemetry,
                           mean.phi = runif(1, 0.9, 0.999), ### two intercepts for juvenile and adults
                           base.obs = rnorm(1,0, 0.001),                # Prior for intercept of observation probability on logit scale
                           base.fail = rnorm(1,0, 0.001),               # Prior for intercept of tag failure probability on logit scale
                           base.recover = rnorm(1,0, 0.001),
                           
                           #sigma.proc.count=runif(2,0.1,5),
                           #mean.lambda.count=runif(2,0.1,2),
                           sigma.obs.count=runif(4,0.1,100),
                          
                           #z.terrvis = cjs.init.z(enchist.terrvis, f.obsvis),
                           mean.p.terrvis=runif(1,0.75, 1),
                           mean.phi.terrvis=runif(2,0.75, 1))}   



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MCMC settings
ns<-1000  	## number of iterations (draws per chain)
nt<-25		  ## thinning rate
nb<-500	  ## length of burn-in (number of iterations discarded at the start of each chain until convergence is reached)
nc<-3		    ## number of chains
nad<-100



try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
EGVU_PVA<- run.jags(data=INPUT, inits=initIPM, monitor=paraIPM,
                      model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_IPM_2023.jags",
                      n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                      method = "rjparallel") 

failed.jags(c('model'))

out<-as.data.frame(summary(EGVU_PVA))
out$parameter<-rownames(summary(EGVU_PVA))
fwrite(out,"EGVU_IPM_estimates.csv")

MCMCsummary(EGVU_PVA$mcmc)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TEST  OF SSM MODEL 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EV.gof<-as.data.frame(summary(EGVU_PVA, vars=c("Dmape.rep","Dmape.obs")))
EV.gof$parameter<-as.character(row.names(EV.gof))

OBS <- EGVU_PVA$mcmc[[1]][,141]
REP <- EGVU_PVA$mcmc[[1]][,140]
tibble(Rep=REP,obs=OBS) %>%
  ggplot(aes(x=REP,y=OBS)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

PVAL=round(mean(REP > OBS),2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TEST OF TERRITORY SURVIVAL WITH MCMCvis package
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## USING THE CALCULATED FIT VALUES FROM THE JAGS MODEL
OBS <- MCMCpstr(EGVU_PVA, params=c("fit.terr"), type="chains")
REP <- MCMCpstr(EGVU_PVA, params=c("fit.rep.terr"), type="chains")
GOF<-tibble(Rep=as.numeric(REP[[1]]),Obs=as.numeric(OBS[[1]])) %>%
  mutate(P=ifelse(Obs>Rep,1,0))

ggplot(GOF,aes(x=Rep,y=Obs)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  annotate("text",label=as.character(mean(GOF$P)),x=50,y=60)

mean(GOF$P)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORT THE MODEL OUTPUT [INCLUDES MODEL INPUT DATA]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#load("EGVU_IPM_input.RData")
try(setwd("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\EGVU_PVA_CaptiveRelease"), silent=T)
load("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\EGVU_PVA_CaptiveRelease\\EGVU_IPM2020_output_v4_FINAL.RData")
out<-as.data.frame(NeoIPM.SIMPLE$summary)  ## changed from NeoIPMbasic
out$parameter<-row.names(NeoIPM.SIMPLE$summary) ## changed from NeoIPMbasic



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOOK UP SCENARIOS FROM MATRIX OF RELEASE AND SURVIVAL SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## give the projections proper scenario labels
## fut.lam[capt.release,surv.imp] - need to extract scenarios from those indices
#capt.release=seq(0,15,1)
#imp.surv=c(1,1.02,1.04,1.06,1.08,1.10)
dim(capt.rel.mat)
head(surv.inc.mat)

## MELT MATRICES TO RECREATE SCENARIOS FROM INDEX NUMBERS
ncr.lu<-capt.rel.mat %>% gather(key=Scenario, value=n.released,-Year) %>%
  filter(Year==1) %>%
  separate(Scenario,sep="_",into=c("n.rel","n.years")) %>%
  mutate(capt.index=seq_along(Year)) %>%
  select(capt.index,n.rel,n.years)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TABLE 1 FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out$parameter
head(out)
TABLE1<-out %>% filter(parameter %in% c('mean.lambda','mean.phi.terrvis[1]','mean.phi.terrvis[2]',
                                        'ann.phi.juv.telemetry[1]',"ann.phi.capt.rel.first.year[1]",'ann.phi.sec.telemetry[1]','ann.phi.third.telemetry[1]',
                                        'ann.phi.juv.telemetry[17]',"ann.phi.capt.rel.first.year[17]",'ann.phi.sec.telemetry[17]','ann.phi.third.telemetry[17]',
                                        'mu.fec')) %>%
  select(parameter,Median,Lower95,Upper95)
TABLE1$Parameter<-c("fecundity","delayed release first year survival","wild first year survival","wild second year survival", "wild third year survival","adult survival (good year)","adult survival (poor year)","population growth rate")
TABLE1
#fwrite(TABLE1,"Table1.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
EV.past<-out[(grep("Nterr\\[",out$parameter)),c(12,2,1,3)] %>%
  mutate(Year=seq(2006,(2005+INPUT$T.count),1)) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  dplyr::select(parm,Year,median,lcl,ucl)

## retrieve the population projections and insert proper scenario labels
EV.fut<-out[(grep("Nterr.f",out$parameter)),c(12,2,1,3)] %>%
  mutate(Year=rep(seq(max(EV.past$Year)+1,(max(EV.past$Year)+INPUT$PROJECTION),1),each=dim(ncr.lu)[1])) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  dplyr::select(parm,Year,median,lcl,ucl)

### SUMMARISE FOR BASELINE TRAJECTORY
EV.base <- EV.fut %>% filter(n.rel==0 & surv.inc<1.0001 & n.years==10 & lag.time==10) %>%
  select(parm, median, lcl, ucl, Year) %>%
  bind_rows(EV.past) %>%
  arrange(Year)

### FILTER FOR FUTURE PLOT TRAJECTORY
EV.plot <- EV.fut %>% 
  bind_rows(EV.past) %>%
  arrange(Year)

### ADD SUMMARISED COUNT DATA
EV.count<-countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  mutate(prop=ifelse(is.na(N),0,country.props[match(Country,names(country.props))])) %>%
  group_by(year) %>%
  summarise(N=sum(N,na.rm=T),prop=sum(prop)) %>%
  mutate(TOT=N/prop) %>%
  mutate(surv=dplyr::lead(TOT)/TOT) 


## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("cornflowerblue", "firebrick"))
#colfunc <- colorRampPalette(c("goldenrod", "darkolivegreen"))
#colfunc <- colorRampPalette(c("black", "lightgray"))


ggplot()+
  geom_ribbon(data=EV.plot,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_line(data=EV.plot, aes(x=Year, y=median),linewidth=1)+
  geom_point(data=EV.count,aes(x=year, y=TOT), size=2,col='darkblue')+

  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,220),breaks=seq(0,210,30), labels=as.character(seq(0,210,30)))+
  scale_x_continuous(name="Year", breaks=seq(2005,2050,5), labels=as.character(seq(2005,2050,5)))+

  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=16),
        legend.text=element_text(size=12, color="black"),
        legend.background=element_blank(),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"),
        legend.position = c(0.90, 0.82))

try(setwd("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\EGVU_PVA_CaptiveRelease\\JApplEcol_attempt"), silent=T)

ggsave("Figure_1.jpg", width=180,height=131, quality=100, units="mm")
ggsave("Figure_1.eps", width=180,height=131, device=cairo_ps, fallback_resolution = 1000, units="mm")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 2: PLOT FUT POP GROWTH AGAINST NUMBER OF RELEASED CHICKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the future lambdas
selcol<-grep("fut.lambda",dimnames(EGVU_PVA$mcmc[[1]])[[2]])
fut.lambda<-EGVU_PVA$mcmc[[1]][,selcol]
for (c in 2:nc){
  fut.lambda<-rbind(fut.lambda,EGVU_PVA$mcmc[[c]][,selcol])
}


FUTLAM<-as.data.frame(fut.lambda) %>% gather(key="parm",value="f.lam") %>%
  group_by(parm) %>%
  summarise(median=quantile(f.lam,0.5),lcl=quantile(f.lam,0.025),ucl=quantile(f.lam,0.975)) %>%
  arrange(median)
tail(FUTLAM)
#fwrite(FUTLAM,"TableX1_FutureGrowthRates.csv")

### CREATE PLOT ###
FUTLAM %>% #mutate(n.rel=as.numeric(n.chicks.released.per.year)) %>%
  filter(lag.time=="10") %>%
  mutate(surv.inc.ord=factor(surv.inc, levels = c("none","+ 2 %","+ 4 %","+ 6 %","+ 8 %"))) %>%
  mutate(n.years=sprintf("for %s years",n.years)) %>%
  arrange(n.rel) %>%
  
  
  ggplot(futlam)+
  geom_hline(aes(yintercept=1), color='red', size=1)+
  geom_line(aes(x=n.rel,y=median),size=1)+
  geom_ribbon(aes(x=n.rel, ymin=lcl,ymax=ucl),alpha=0.2)+
  facet_grid(n.years~surv.inc.ord) +
  
  ylab("Future population growth rate") +
  xlab("Number of captive-bred birds released per year") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("Figure_3.jpg", width=180,height=131, quality=100, units="mm")
ggsave("Figure_3.eps", width=180,height=131, units="mm", device=cairo_ps, fallback_resolution = 1000)








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TABLE 2: EXTINCTION PROBABILITY WITHOUT SURVIVAL IMPROVEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## THIS PROCESS WILL TAKE 5 MIN
extprop <- data.frame()

for(scen in 1:nrow(ncr.lu)){
  
  ### FIND COLUMS WE NEED
  colname<-sprintf("Nterr.f\\[%s,",scen)
  selcol<-grep(colname,dimnames(FOCMOD$samples[[1]])[[2]])
  
  allchainsamples <- data.frame()
  for(chain in 1:4) {
    
    ### EXTRACT AND SUMMARISE DATA
    samplesout<-as.data.frame(FOCMOD$samples[[1]][,selcol]) %>% gather(key="parm", value="value")
    allchainsamples <- rbind(allchainsamples,as.data.frame(samplesout))
  }
  
  ### CALCULATE EXTINCTION PROBABILITY
  allchainsamples<- allchainsamples %>%
    mutate(capt.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])) %>%
    mutate(surv.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
    mutate(Year=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,3])+(max(trendinput$year))) %>%
    
    mutate(n=1, inc=ifelse(value<25,1,0)) %>%
    group_by(capt.index,surv.index,Year) %>%
    summarise(ext.prob=sum(inc)/sum(n))
  
  extprop <- rbind(extprop,as.data.frame(allchainsamples))
  print(scen)
}

head(samplesout)
head(extprop)
dim(extprop)


## create factors for plot labels and order them appropriately

extprop<- extprop %>%
  left_join(ncr.lu, by="capt.index") %>%
  left_join(surv.lu, by="surv.index") %>%
  filter(lag.time==10) %>%
  mutate(surv.inc=ifelse(as.numeric(surv.inc)>1,paste("+",as.integer((as.numeric(surv.inc)-1)*100),"%"),"none")) %>%
  mutate(surv.inc.ord=factor(surv.inc, levels = c("none","+ 2 %","+ 4 %","+ 6 %","+ 8 %","+ 10 %"))) %>%
  mutate(n.rel.ord=factor(n.rel, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))) %>%
  mutate(n.years=sprintf("for %s years",n.years)) %>%
  arrange(n.rel,Year) %>%
  mutate(release=paste(n.rel,n.years," "))

dim(extprop)



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########## 
### CREATE FIGURE 3 INSTEAD OF TABLE 2 (REQUESTED BY REVIEWER) ###
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########## 

TABLE2 %>% mutate(ext.prob=ifelse(n.rel==0,47.6,ext.prob)) %>%
  mutate(n.rel=as.numeric(n.rel)) %>%
  
  ggplot()+
  geom_line(aes(x=n.rel,y=ext.prob, color=n.years),size=1)+
  ylab("Extinction probability in 2049") +
  xlab("Number of captive-bred birds released per year") +
  
  ## format axis ticks
  scale_y_continuous(name="Extinction probability in 2049", limits=c(0,50),breaks=seq(0,50,10), labels=as.character(paste(seq(0,50,10),"%",sep="")))+
  scale_x_continuous(name="Number of captive-bred birds released per year", breaks=seq(0,15,3), labels=as.character(seq(0,15,3)))+
  guides(color=guide_legend(title="Duration \n of releases"))+
  scale_color_viridis_d(begin = .4, direction = -1) +
  #scale_colour_manual(palette=colfunc)+
  #scale_fill_manual(palette=colfunc)+
  
  ## other formatting
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"),
        legend.position = c(0.88, 0.84))


ggsave("Figure_2.jpg", width=180,height=144, units="mm",quality=100)
ggsave("Figure_2.eps", width=180,height=144, units="mm",device=cairo_ps, fallback_resolution = 1000)












#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INTEGRATED POPULATION MODEL WITH COUNTRY-SPECIFIC COUNT INPUT AND RANDOMLY VARIABLE ADULT SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)

sink("EGVU_IPM_2023.jags")
cat("
  model {
    #-------------------------------------------------
    # integrated population model for the balkan population
    # - age structured model with 6 age classes: 
    # - age-specific probability to recruit at ages 4 or 5
    # - age-specific survival derived from tracking data
    # - adult survival based on territory occupancy; CHANGED in 2023 to not include monthly occasions and effort
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - FUTURE PROJECTION WITH ONE SCENARIO
    # - productivity supplemented by captive-bred juveniles (0-10 per year)
    # - included GOF test for SSM
    #-------------------------------------------------
    
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    # Priors and constraints FOR FECUNDITY
    mu.fec ~ dunif(0,2)           # Priors on fecundity can range from 0- 2 chicks per pair (uninformative)
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
      sigma.obs.count[s] ~ dunif(0.1,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    # Priors and constraints FOR JUVENILE SURVIVAL FROM TELEMETRY
       #### MONTHLY SURVIVAL PROBABILITY
    for (i in 1:nind){
      for (t in f.telemetry[i]:(n.occasions)){
        logit(phi[i,t]) <- lp.mean.telemetry +      ### intercept for mean survival 
            b.phi.capt*(capt[i]) +     ### survival dependent on captive-release (captive-raised or other)
            b.phi.mig*(mig[i,t]) +     ### survival dependent on migration or not (in this month)
            b.phi.age*(age[i,t]) +    ### survival dependent on age (juvenile or other)
            b.phi.LIFE*(LIFE[i,t]) +     ### survival dependent on LIFE project
            b.phi.long*(long[i])*(juv[i,t]) ### survival varies for juvenile first-time migrants when they originate from farther west
      } #t
    } #i
    
    #### BASELINE FOR SURVIVAL PROBABILITY (wild adult stationary from east)
    mean.phi ~ dunif(0.9, 1)   # uninformative prior for all MONTHLY survival probabilities
    lp.mean.telemetry <- log(mean.phi/(1 - mean.phi))    # logit transformed survival intercept

    #### SLOPE PARAMETERS FOR SURVIVAL PROBABILITY
    b.phi.age ~ dunif(0, 5)           # Prior for age effect on survival probability on logit scale - should be positive
    b.phi.capt ~ dnorm(0, 0.001)         # Prior for captive release on survival probability on logit scale
    b.phi.mig ~ dunif(-2,0)         # Prior for COST OF MIGRATION migration on survival probability on logit scale
    b.phi.LIFE ~ dnorm(0, 0.001)         # Prior for LIFE effect on survival probability on logit scale    
    b.phi.long ~ dnorm(0, 0.001)         # Prior for longitudinal effect on survival probability on logit scale  
    
    #### TAG FAILURE AND LOSS PROBABILITY
    for (i in 1:nind){
      for (t in f.telemetry[i]:(n.occasions)){
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
    
    
        # Priors and constraints FOR ADULT SURVIVAL FROM TERRITORY MONITORING
    
    ## Priors for detection probability
    lmu.p.terrvis <- log(mean.p.terrvis/(1 - mean.p.terrvis))    # logit transformed detection intercept
    mean.p.terrvis ~ dunif(0, 1)			
    
    # Det prob relationship with observation effort
    beta.obs.eff.terrvis ~ dnorm(0, 0.0001)
    
    
    # Priors for survival
    for (nypterr in 1:2){   ## only 2 survival periods
      mean.phi.terrvis[nypterr] ~ dunif(0.5, 1)   # informative prior for annual survival probabilities BEFORE/AFTER 2016
      #mean.rec.terrvis[nypterr] ~ dunif(0, 1)   # informative prior for annual recruitment probabilities BEFORE/AFTER 2016
    }    
    # Prior for recruitment
   mean.rec.terrvis ~ dunif(0, 1)   # informative prior for annual recruitment probabilities BEFORE/AFTER 2016 
    
    
    ### RANDOM OBSERVATION EFFECT FOR EACH TERRITORY AND YEAR
    
    for (nyRpterr in 1:nprim.terrvis){
      for (nsit in 1:nsite.terrvis){
        rand.obs.terrvis[nsit,nyRpterr] ~ dnorm(0, tau.obs.terrvis)
        }
    }
    
    tau.obs.terrvis <- 1 / (sd.obs.terrvis * sd.obs.terrvis)
    sd.obs.terrvis ~ dunif(0, 3)
    
    
    # Priors for population process and immigration
    
    # Initial population sizes for first year of monitoring
    nestlings[1] ~ dunif(50, 100)   ##changed from JUV
    N1[1] ~ dunif(10, 35)
    N2[1] ~ dunif(5, 20)
    N3[1] ~ dunif(3, 15)
    N4[1] ~ dunif(0, 10)
    N5[1] ~ dunif(0, 5)
    N6[1] ~ dunif(190, 250)
    
    
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T.count){
      
      nestlings[tt] <- mu.fec * 0.5 * Nterr[tt]                                                              ### number of local recruits
      N1[tt]  ~ dbin(ann.phi.juv.telemetry[period[tt]], round(nestlings[tt-1]))                                                    ### number of 1-year old survivors - add CAPT.ADD in here
      N2[tt] ~ dbin(ann.phi.sec.telemetry[period[tt]], round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.phi.third.telemetry[period[tt]], round(N2[tt-1]))                                                    ### number of 3-year old survivors
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
        y.count[tlc,s] ~ dnorm(Nterr.country[tlc,s], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }

      

      # -------------------------------------------------        
      # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods (or binom for single broods after egg harvest)
      # -------------------------------------------------
      
      J.fec[tlc] ~ dpois(rho.fec[tlc])
      rho.fec[tlc] <- R.fec[tlc]*mu.fec
    } #	close loop over every year in which we have count and fecundity data
    
    
    
      # -------------------------------------------------        
      # 2.4. Goodness of Fit test for the SSM part (and the whole IPM)
      # -------------------------------------------------
          
    
        # Assessing the fit of the state-space model
    # 1. Compute fit statistic for observed data
    # Discrepancy measure: mean absolute error
    for (tgof in 1:T.count){
      C.exp[tgof] <- Nterr[tgof]                            # Expected counts
      C.act[tgof] <- sum(y.count[tgof,])
      Dssm.obs[tgof] <- abs((C.act[tgof] - C.exp[tgof]) / C.act[tgof])           # Discrepancy measure
    }
    Dmape.obs <- sum(Dssm.obs)

    # 2. Compute fit statistic for replicate data
    # Discrepancy measure: mean absolute error
    for (tgof in 1:T.count){
        for (sgof in 1:countries){			### start loop over every country
         C.rep.c[tgof,sgof] ~ dnorm(Nterr.country[tgof,sgof], tau.obs.count[sgof])								# Distribution for random error in observed numbers (counts)
        }
      C.rep[tgof] <- sum(C.rep.c[tgof,])                 # Generate replicate data
      Dssm.rep[tgof] <- abs((C.rep[tgof] - C.exp[tgof]) / C.rep[tgof])   # Discrepancy measure
    }
    Dmape.rep <- sum(Dssm.rep)
    
    
    
    # -------------------------------------------------        
    # 2.5. Likelihood for juvenile survival from telemetry
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind){
      
      for (t in f.telemetry[i]:(n.occasions-1)){
        
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
      z.telemetry[i,f.telemetry[i]] <- 2 ## alive when first marked
      for (t in (f.telemetry[i]+1):n.occasions){
        # State process: draw S(t) given S(t-1)
        z.telemetry[i,t] ~ dcat(ps[z.telemetry[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y.telemetry[i,t] ~ dcat(po[z.telemetry[i,t], i, t-1,])
      } #t
    } #i

    
    # -------------------------------------------------        
    # 2.6. Likelihood for adult survival from territory monitoring
    # -------------------------------------------------
    
    ### ECOLOGICAL STATE MODEL WITH ESTIMATE OF SURVIVAL
    
    for (ilterr in 1:nsite.terrvis){
      x.terrvis[ilterr,f.obsvis[ilterr]] <- 2 #firstobs[ilterr]
    
      for (klterr in (f.obsvis[ilterr]+1):nprim.terrvis){
        z.terrvis[ilterr,klterr] ~ dbin(mean.phi.terrvis[phase[klterr]],x.terrvis[ilterr,klterr-1])
        mu.x[ilterr,klterr] ~ dbin(mean.rec.terrvis,(2-z.terrvis[ilterr,klterr]))
        x.terrvis[ilterr,klterr] <- mu.x[ilterr,klterr] + z.terrvis[ilterr,klterr]

      } 						# close klterr loop over primary period  - years
    } 							# close ilterr loop over sites
    
    ### OBSERVATION MODEL WITH RANDOM EFFECT OF YEAR

    for (iobs in 1:nsite.terrvis){
      for (kobs in f.obsvis[iobs]:nprim.terrvis){
        linpred.terrvis[iobs,kobs] <- lmu.p.terrvis+beta.obs.eff.terrvis*eff.terrvis[iobs,kobs] + rand.obs.terrvis[iobs,kobs]
        p.terrvis[iobs,kobs] <- 1 / (1 + exp(-linpred.terrvis[iobs,kobs]))
        y.terrvis[iobs,kobs] ~ dbin(p.terrvis[iobs,kobs], x.terrvis[iobs,kobs])
        
        
        ## GOODNESS OF FIT TEST SECTION
        ## Expected occupancy at site iobs, survey kobs
        exp.terr[iobs,kobs] <- p.terrvis[iobs,kobs] * x.terrvis[iobs,kobs] 
    
        ## Discrepancy 
        ## (note small value added to denominator to avoid potential divide by zero)
        E.obs.terr[iobs,kobs] <- pow((y.terrvis[iobs,kobs] - exp.terr[iobs,kobs]), 2) / (exp.terr[iobs,kobs] + 0.5)

        ## Simulate new count from model
        y.rep.terr[iobs,kobs] ~ dbin(p.terrvis[iobs,kobs], x.terrvis[iobs,kobs])

        ## Discrepancy 
        E.rep.terr[iobs,kobs] <- pow((y.rep.terr[iobs,kobs] - exp.terr[iobs,kobs]), 2) / (exp.terr[iobs,kobs] + 0.5)

      } 						 # close kobs loop over year
      fit.terr.temp[iobs] <- sum(E.obs.terr[iobs,f.obsvis[iobs]:nprim.terrvis])
      fit.rep.terr.temp[iobs] <- sum(E.rep.terr[iobs,f.obsvis[iobs]:nprim.terrvis])
    } 							# close iobs loop over sites
    
    # chi-squared test statistics
    fit.terr <- sum(fit.terr.temp[])
    fit.rep.terr <- sum(fit.rep.terr.temp[])
    
    
    
    # -------------------------------------------------        
    # 3. DERIVED PARAMETERS
    # -------------------------------------------------
    
    ### 3.1 TELEMETRY DERIVED SURVIVAL ESTIMATES

  for(LIF in c(1:2)){    
    ## for WILD BIRDS
    for (ageprog in 1:36){
      logit(phi.wild.telemetry[ageprog,LIF]) <- lp.mean.telemetry +      ### intercept for mean survival 
        b.phi.mig*(migprog.age[ageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*ageprog +    ### survival dependent on age (juvenile or other)
        b.phi.LIFE*(LIF-1) +     ### survival dependent on LIFE project
        b.phi.long*(0.39219364)
    }

    ## for CAPTIVE-REARED DELAYED RELEASE BIRDS
    for (captageprog in 1:36){
      logit(phi.capt.telemetry[captageprog,LIF]) <- lp.mean.telemetry +      ### intercept for mean survival
        b.phi.capt +     ### survival dependent on captive-release (captive-raised or other)
        b.phi.mig*(captprog.age[captageprog]) +     ### survival dependent on migration (first-time crossing of sea)
        b.phi.age*captageprog +    ### survival dependent on age (juvenile or other)
        b.phi.LIFE*(LIF-1)  +     ### survival dependent on LIFE project
        b.phi.long*(0.39219364)
    }
    
    
    #### SUMMARISE ANNUAL SURVIVAL PROBABILITY
    ann.phi.juv.telemetry[tt]<-prod(phi.wild.telemetry[1:12,LIF])
    # phi.wild.telemetry[2,phase[tt]]*
    # phi.wild.telemetry[3,phase[tt]]*
    # phi.wild.telemetry[4,phase[tt]]*
    # phi.wild.telemetry[5,phase[tt]]*
    # phi.wild.telemetry[6,phase[tt]]*
    # phi.wild.telemetry[7,phase[tt]]*
    # phi.wild.telemetry[8,phase[tt]]*
    # phi.wild.telemetry[9,phase[tt]]*
    # phi.wild.telemetry[10,phase[tt]]*
    # phi.wild.telemetry[11,phase[tt]]*
    # phi.wild.telemetry[12,phase[tt]]					### multiply monthly survival from age 1-12
    
    ann.phi.sec.telemetry[tt]<-prod(phi.wild.telemetry[13:24,LIF])
    # phi.wild.telemetry[14,phase[tt]]*
    # phi.wild.telemetry[15,phase[tt]]*
    # phi.wild.telemetry[16,phase[tt]]*
    # phi.wild.telemetry[17,phase[tt]]*
    # phi.wild.telemetry[18,phase[tt]]*
    # phi.wild.telemetry[19,phase[tt]]*
    # phi.wild.telemetry[20,phase[tt]]*
    # phi.wild.telemetry[21,phase[tt]]*
    # phi.wild.telemetry[22,phase[tt]]*
    # phi.wild.telemetry[23,phase[tt]]*
    # phi.wild.telemetry[24,phase[tt]]					### multiply monthly survival from age 13-24
    ann.phi.third.telemetry[tt]<-prod(phi.wild.telemetry[25:36,LIF]) ##*phi.wild.telemetry[26,phase[tt]]*phi.wild.telemetry[27,phase[tt]]*phi.wild.telemetry[28,phase[tt]]*phi.wild.telemetry[29,phase[tt]]*phi.wild.telemetry[30,phase[tt]]*phi.wild.telemetry[31,phase[tt]]*phi.wild.telemetry[32,phase[tt]]*phi.wild.telemetry[33,phase[tt]]*phi.wild.telemetry[34,phase[tt]]*phi.wild.telemetry[35,phase[tt]]*phi.wild.telemetry[36,phase[tt]]					### multiply monthly survival from age 25-36
    
    ann.phi.capt.rel.first.year[tt]<-prod(phi.capt.telemetry[8:24,LIF]) ##*phi.capt.telemetry[9,phase[tt]]*phi.capt.telemetry[10,phase[tt]]*phi.capt.telemetry[11,phase[tt]]*phi.capt.telemetry[12,phase[tt]]*phi.capt.telemetry[13,phase[tt]]*phi.capt.telemetry[14,phase[tt]]*phi.capt.telemetry[15,phase[tt]]*phi.capt.telemetry[16,phase[tt]]*phi.capt.telemetry[17,phase[tt]]*phi.capt.telemetry[18,phase[tt]]*phi.capt.telemetry[19,phase[tt]]*phi.capt.telemetry[20,phase[tt]]*phi.capt.telemetry[21,phase[tt]]*phi.capt.telemetry[22,phase[tt]]*phi.capt.telemetry[23,phase[tt]]*phi.capt.telemetry[24,phase[tt]]					### first year for delayed-release bird is longer, from month 8 to 24 
    
  } ### end loop over LIFE-specific survival
    
    
    ### 3.2 POPULATION GROWTH DERIVED FROM COUNTS    
    
    # Annual population growth rate
    for (ipop in 1:(T.count-1)){
      lambda.t[ipop] <- Nterr[ipop+1] / max(Nterr[ipop],1)           ### inserted max of 1 to avoid this becoming an invalid parent error
      loglambda.t[ipop]<-log(lambda.t[ipop])## for calculating geometric mean of overall population growth rate
    }
    
    #### OVERALL POPULATION GROWTH RATE  #########
    
    baseline.lambda<-exp((1/(T.count-6))*sum(loglambda.t[1:(T.count-6)]))   # Geometric mean
    LIFE.lambda<-exp((1/5)*sum(loglambda.t[(T.count-5):(T.count-1)]))   # Geometric mean    
    
    
    
    # -------------------------------------------------        
    # 4. PREDICTION INTO THE FUTURE
    # -------------------------------------------------
    ### INCLUDE SCENARIOS FOR N CAPTIVE BIRDS AND SURVIVAL IMPROVEMENT
    
    # CAPTIVE RELEASE OF JUVENILE BIRDS
    for (ncr in 1:scen.capt.release){
      
      ### FUTURE FECUNDITY IS BASED ON mu.fec
      
      for (fut in 1:PROJECTION){
        fut.fec[fut,ncr] <-mu.fec  ## no chicks taken after 5 years anymore
        fut.survival[ncr,fut] <-mean.phi.terrvis[2] ### future survival remains as under LIFE project
      }
        
        ## POPULATION PROCESS
        ## need to copy previous array elements
        nestlings.f[ncr,1] <- round(nestlings[T.count])   ##JUV[T.count]
        N1.f[ncr,1] <- N1[T.count]
        N2.f[ncr,1] <- N2[T.count]
        N3.f[ncr,1] <- N3[T.count]
        N4.f[ncr,1] <- N4[T.count]
        N5.f[ncr,1] <- N5[T.count]
        N6.f[ncr,1] <- N6[T.count]
        Nterr.f[ncr,1] <- Nterr[T.count]
        N2wild.f[ncr,1] <- 0
        N2released.f[ncr,1] <- 0
        
        for (fut in 2:PROJECTION){
          
          ### probabilistic formulation
          nestlings.f[ncr,fut] <- (fut.fec[fut,ncr] * 0.5 * Nterr.f[ncr,fut])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
          N1.f[ncr,fut] ~ dbin(ann.phi.juv.telemetry[2],round(nestlings.f[ncr,fut-1]))             ### +rescued[ncr,fut] number of 1-year old survivors
          N2wild.f[ncr,fut] ~ dbin(ann.phi.sec.telemetry[2],round(N1.f[ncr,fut-1]))                ### number of 2-year old wild survivors
          N2released.f[ncr,fut] ~ dbin(ann.phi.capt.rel.first.year[2],round(capt.release[fut]))             ### +rescued[ncr,fut] number of 1-year old survivors 
          N2.f[ncr,fut] <-  N2wild.f[ncr,fut] + N2released.f[ncr,fut]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
          N3.f[ncr,fut] ~ dbin(ann.phi.third.telemetry[2],round(N2.f[ncr,fut-1]))                                                    ### number of 3-year old survivors
          N4.f[ncr,fut] ~ dbin(fut.survival[ncr,fut],round(N3.f[ncr,fut-1]))                                                       ### number of 4-year old survivors
          N5.f[ncr,fut] ~ dbin(fut.survival[ncr,fut],round(N4.f[ncr,fut-1]))                                                       ### number of 5-year old survivors
          N6.f[ncr,fut] ~ dbin(fut.survival[ncr,fut],round((N5.f[ncr,fut-1]+N6.f[ncr,fut-1])))                                   ### number of 6-year or older (adult) birds
          Nterr.f[ncr,fut] <- round((N4.f[ncr,fut] * 0.024) + (N5.f[ncr,fut] * 0.124) + (N6.f[ncr,fut]))
          
          
        } # fut
        
        for (fut2 in 1:(PROJECTION-1)){
          lambda.t.f[ncr,fut2] <- Nterr.f[ncr,fut2+1] / max(Nterr.f[ncr,fut2],1)
          loglambda.t.f[ncr,fut2]<-log(lambda.t.f[ncr,fut2])## for calculating geometric mean of overall population growth rate
        } # fut2
        
        
        #### FUTURE POPULATION GROWTH RATE  #########
        fut.lambda[ncr]<-exp((1/(PROJECTION-1))*sum(loglambda.t.f[ncr,1:(PROJECTION-1)]))   # Geometric mean
      
    } # end scenario of n capt released
    
  } # END of the JAGS model
    ",fill = TRUE)
sink()
