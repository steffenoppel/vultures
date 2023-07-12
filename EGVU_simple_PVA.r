##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS
#
##########################################################################
# based on Oppel et al 2021 (Journal of Applied Ecology)
# update in June 2023 to address reviewer comment

### after implementation of IPM failed, used simple projection matrix in addition to the SSM
### projected survival values are highly dependent on longitude

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
library(MCMCvis)
library(popbio)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN OUTPUT FROM DIFFERENT SUBMODELS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE"), silent=T)
load("Analysis/EGVU_LIFE_telemetry_survival_output_2023.RData")
load("Analysis/EGVU_LIFE_territory_survival_output_2022.RData")

### COMBINE ALL ESTIMATES INTO A SINGLE MATRIX
surv.parm<-surv.parm %>% filter(AgeGroup=="third") %>% ungroup() %>%
  mutate(AgeGroup="adult") %>%
  mutate(mean.surv=as.numeric(EV_ad_surv$summary$statistics[1:2,1]),
         sd=as.numeric(EV_ad_surv$summary$statistics[1:2,2])) %>%
  bind_rows(surv.parm) %>%
  select(LIFE,Origin,AgeGroup,mean.surv,sd) %>%
  ungroup() %>%
  arrange(AgeGroup,Origin,LIFE)
surv.parm<-as.data.frame(surv.parm)  ## necessary because the indexing doesn't work on tibbles
rm(list=setdiff(ls(), c("ad.phi.mean","ad.phi.sd","surv.parm","filter","select")))
fwrite(surv.parm,"output/EGVU_survival_estimates_LIFE.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN INPUT DATA FOR TREND MODEL 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Data"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)
head(breed)
countrytrendinput<- fread("EVGU_countrytrendinput2019.csv") %>%
  dplyr::filter(year<2018)

countrytrendinput<- occu %>% filter(year>2017) %>%
  filter(!(Country=="Niger")) %>%    # introduced in 2019 because database now has data from Niger
  group_by(year,Country) %>%
  summarise(N=sum(occupancy)) %>%
  spread(key=Country,value=N) %>%
  bind_rows(countrytrendinput) %>%
  arrange(year)


# ### ADD DATA FROM 2022
# #countrytrendinput<- fread("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\SMARTbirds\\EGVU_D1_Terr_Occupancy_2022.csv") %>%
# countrytrendinput<- fread("EGVU_D1_Terr_Occupancy_2022.csv") %>%
#   mutate(year=2022) %>%
#   group_by(year,Country) %>%
#   summarise(N=sum(N_adults)) %>%
#   spread(key=Country,value=N) %>%
#   rename(Albania=ALB,Bulgaria=BG,Greece=GR,`North Macedonia`=NMD) %>%
#   bind_rows(countrytrendinput) %>%
#   arrange(year)
# countrytrendinput


## calculate mean proportion of each country to total count
country.props<-colMeans(prop.table(as.matrix(countrytrendinput[,-1]), margin = 1), na.rm=T)


### PLOT RAW COUNT DATA 

countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  ggplot()+
  geom_line(aes(x=year, y=N, col=Country), linewidth=1)

### INSPECT TOTAL

totcount<-countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  mutate(prop=ifelse(is.na(N),0,country.props[match(Country,names(country.props))])) %>%
  group_by(year) %>%
  summarise(N=sum(N,na.rm=T),prop=sum(prop)) %>%
  mutate(TOT=N/prop) %>%
  mutate(surv=dplyr::lead(TOT)/TOT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET FUTURE PROJECTION SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

capt.release=5
project.time=20
imp.surv=c(1)
PROJECTION.years<-seq(1,30,1)

capt.rel.mat<- expand.grid(PROJECTION.years,capt.release,project.time) %>%
  rename(Year=Var1,REL=Var2,DUR=Var3) %>%
  mutate(scenario=paste(REL,DUR, sep="_")) %>%
  mutate(REL=ifelse(Year>DUR,0,REL)) %>%
  select(-DUR) %>%
  spread(key=scenario, value=REL)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DATA AND SAVE WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### create function to estimate alpha and beta from mean and sd
estBetaParams <- function(mu, sd) {
  alpha <- ((1 - mu) / (sd^2) - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

### test whether it works
# estBetaParams(mu=0.333, sd=0.162)
# hist(rbeta(100,2.485,4.978))


## Bundle data across all of the data sources

INPUT <- list(### STATE SPACE MODEL FOR TREND
              y.count=as.matrix(countrytrendinput[,2:5]),
              T.count=dim(countrytrendinput)[1],
              #country.prop=as.numeric(country.props),  ## proportion of count data that come from countries ALB, BG, GR, NMD
              countries=4,
              period=ifelse(countrytrendinput$year<2018,1,2),
              
              ### Future Projection and SCENARIOS FOR TRAJECTORY
              ## provide data for beta distribution where mean=a and sd=b
              ad.phi.mean=as.numeric(estBetaParams(surv.parm[1:2,4],surv.parm[1:2,5])$alpha),
              ad.phi.sd=as.numeric(estBetaParams(surv.parm[1:2,4],surv.parm[1:2,5])$beta),   ## precision=1/variance and variance = SD^2
              juv.phi.mean=as.numeric(estBetaParams(surv.parm[3:4,4],surv.parm[3:4,5])$alpha),
              juv.phi.sd=as.numeric(estBetaParams(surv.parm[3:4,4],surv.parm[3:4,5])$beta),   ## precision=1/variance and variance = SD^2
              sec.phi.mean=as.numeric(estBetaParams(surv.parm[7:8,4],surv.parm[7:8,5])$alpha),
              sec.phi.sd=as.numeric(estBetaParams(surv.parm[7:8,4],surv.parm[7:8,5])$beta),   ## precision=1/variance and variance = SD^2
              third.phi.mean=as.numeric(estBetaParams(surv.parm[9:10,4],surv.parm[9:10,5])$alpha),
              third.phi.sd=as.numeric(estBetaParams(surv.parm[9:10,4],surv.parm[9:10,5])$beta),   ## precision=1/variance and variance = SD^2
              sec.capt.phi.mean=as.numeric(estBetaParams(surv.parm[5:6,4],surv.parm[5:6,5])$alpha),
              sec.capt.phi.sd=as.numeric(estBetaParams(surv.parm[5:6,4],surv.parm[5:6,5])$beta),   ## precision=1/variance and variance = SD^2
              
              PROJECTION=30,                ## used 10 and 50 years previously, now trying 30 years
              capt.release=capt.rel.mat[,2:ncol(capt.rel.mat)])


### test whether the beta priors make sense works
hist(rbeta(100,INPUT$ad.phi.mean,INPUT$ad.phi.sd))
hist(rbeta(100,INPUT$sec.phi.mean,INPUT$sec.phi.sd))

## Parameters to be estimated ('monitored') by JAGS
paraPVA<-c("mean.lambda.count","Ntot","fut.lambda","Nterr.f",
           "Dmape.rep","Dmape.obs") #,"phi.ad","phi.juv","phi.sec","phi.third","mean.fec",
           #"N1.f","N2.f","N3.f","N4.f","N5.f","N6.f")



try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE"), silent=T)
save.image("Data/EGVU_LIFE_PVA_input_2023.RData")
load("Data/EGVU_LIFE_PVA_input_2023.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inits.count<- function() {list(sigma.proc.count=runif(2,0.1,5),
                                mean.lambda.count=runif(2,0.1,2),
                                sigma.obs.count=runif(4,0.1,100))}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## MCMC settings
ns<-50000  	## number of iterations (draws per chain)
nt<-5		  ## thinning rate
nb<-15000	  ## length of burn-in (number of iterations discarded at the start of each chain until convergence is reached)
nc<-3		    ## number of chains
nad<-5000



try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
EGVU_PVA<- run.jags(data=INPUT, inits=inits.count, monitor=paraPVA,
                      model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_PVA_2023.jags",
                      n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                      method = "rjparallel")
# EGVU_PVA<- autojags(data=INPUT, inits=inits.count, monitor=paraPVA,
#                     model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_PVA_2023.jags",
#                     n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T) 

failed.jags(c('model'))

out<-as.data.frame(summary(EGVU_PVA))
out$parameter<-rownames(summary(EGVU_PVA))
fwrite(out,"C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output/EGVU_IPM_estimates.csv")

MCMCsummary(EGVU_PVA$mcmc)
MCMCdiag(EGVU_PVA$mcmc,
         round = 3,
         file_name = 'EGVU_PVA_10Jul2023',
         dir = 'C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output',
         mkdir = 'PVA_10Jul2023',
         add_field = '1.0',
         add_field_names = 'Data version',
         save_obj = TRUE,
         obj_name = 'pva-fit-10Jul2023',
         add_obj = list(INPUT, sessionInfo()),
         add_obj_names = c('ipm-data-7Jul2023', 'session-info-7Jul2023'))

try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output"), silent=T)
save.image("EGVU_LIFE_PVA_output_2023.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TEST  OF SSM MODEL 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OBS <- MCMCpstr(EGVU_PVA$mcmc, params=c("Dmape.obs"), type="chains")
REP <- MCMCpstr(EGVU_PVA$mcmc, params=c("Dmape.rep"), type="chains")
GOF<-tibble(Rep=as.numeric(REP[[1]]),Obs=as.numeric(OBS[[1]])) %>%
  mutate(P=ifelse(Obs>Rep,1,0))

ggplot(GOF,aes(x=Rep,y=Obs)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  annotate("text",label=as.character(mean(GOF$P)),x=1.5,y=1.6)

mean(GOF$P)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TABLE 1 FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out$parameter
head(out)
TABLE1<-out %>% filter(parameter %in% c('baseline.lambda','LIFE.lambda','mean.phi.terrvis[1]','mean.phi.terrvis[2]',
                                        'ann.phi.juv.telemetry[1]',"ann.phi.capt.rel.first.year[1]",'ann.phi.sec.telemetry[1]','ann.phi.third.telemetry[1]',
                                        'ann.phi.juv.telemetry[2]',"ann.phi.capt.rel.first.year[2]",'ann.phi.sec.telemetry[2]','ann.phi.third.telemetry[2]',
                                        'mu.fec')) %>%
  select(parameter,Median,Lower95,Upper95)
TABLE1$Parameter<-c("fecundity","delayed release first year survival","wild first year survival","wild second year survival", "wild third year survival","adult survival (good year)","adult survival (poor year)","population growth rate")
TABLE1
#fwrite(TABLE1,"Table1.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
EV.past<-out[(grep("Ntot\\[",out$parameter)),c(12,2,1,3)] %>%
  mutate(Year=seq(2006,(2005+INPUT$T.count),1)) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  dplyr::select(parm,Year,median,lcl,ucl)

## retrieve the population projections and insert proper scenario labels
EV.fut<-out[(grep("Nterr.f",out$parameter)),c(12,2,1,3)] %>%
  mutate(Year=seq(max(EV.past$Year)+1,(max(EV.past$Year)+INPUT$PROJECTION))) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  dplyr::select(parm,Year,median,lcl,ucl)

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
  geom_vline(aes(xintercept=2040),linewidth=1, linetype="dotted")+
  geom_ribbon(data=EV.plot,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_line(data=EV.plot, aes(x=Year, y=median),linewidth=1)+
  geom_point(data=EV.count,aes(x=year, y=TOT), size=2,col='darkblue')+

  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,500),breaks=seq(0,500,50), labels=as.character(seq(0,500,50)))+
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

# try(setwd("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\EGVU_PVA_CaptiveRelease\\JApplEcol_attempt"), silent=T)
# 
# ggsave("Figure_1.jpg", width=180,height=131, quality=100, units="mm")
# ggsave("Figure_1.eps", width=180,height=131, device=cairo_ps, fallback_resolution = 1000, units="mm")
# 



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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TROUBLESHOOT POPULATION STRUCTURE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(popbio)

#########################################################################
# CALCULATE CRUDE STABLE AGE DISTRIBUTION TO ALLOCATE MORTALITY TO AGE CLASSES BY EXPECTED PROPORTION OF ABUNDANCE
#########################################################################

## take mean survival and fecundity values from paper
#Sx <- c(0.296,0.484,0.580,rep(0.932,3))  ## estimates from Oppel et al. 2021
#Sx <- c(0.296*1.06,0.484*1.06,0.580*1.06,rep(0.957,3))  ## estimates from Oppel et al. 2021 with 6% improvement
Sx <- c(as.numeric(surv.parm[4,4]),as.numeric(surv.parm[8,4]),as.numeric(surv.parm[10,4]),rep(as.numeric(surv.parm[2,4]),3))
Fx <- c(rep(0,3),0.024*1.005*0.5,0.124*1.005*0.5,1.005*0.5)
les.mat <- matrix(rep(0,36),nrow=6)
les.mat[1,] <- Fx
for(i in 1:5){
  les.mat[(i+1),i] <- Sx[i]
}
les.mat[6,6] <- 0.944
n<-as.numeric(c(50),rep(NA,5))   # produce simple initial vector of population sizes based on mean N of juveniles and survival
for (l in 2:5){
  n[l]<-round(n[l-1]*Sx[l-1])
}
n[6]<-100-n[4]*0.024-n[5]*0.124 # populate adults from the sum of estimated total pop size in 2021
n
pop.pro <- pop.projection(A=les.mat, n, 100)
pop.pro$stable.stage
pop.pro$lambda
sum(n)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY POPULATION MODEL WITH COUNTRY-SPECIFIC COUNT INPUT AND RANDOMLY VARIABLE ADULT SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# adapted from IPM for French Egyptian vulture populations Lieury et al. 2016

try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)

sink("EGVU_PVA_2023_2scen.jags")
cat("
  model {

    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
        # Priors and constraints
    
    N.est.count[1,1] ~ dunif(20,40)   ## draw random value from a uniform distribution between 50 and 200 for initial population size 
    N.est.count[1,2] ~ dunif(60,100)   ## draw random value from a uniform distribution between 50 and 200 for initial population size 
    N.est.count[1,3] ~ dunif(30,60)   ## draw random value from a uniform distribution between 50 and 200 for initial population size 
    N.est.count[1,4] ~ dunif(40,70)   ## draw random value from a uniform distribution between 50 and 200 for initial population size 
    
  for (t in 1:2) {   ###allow two separate time trends
    mean.lambda.count[t] ~ dunif(0.1,2)	#Prior for mean growth rate
    sigma.proc.count[t]  ~ dunif(0.1,10)	#Prior for SD of state process (annual variation in pop size)
    sigma2.proc.count[t] <-pow(sigma.proc.count[t] ,2)
    tau.proc.count[t] <-pow(sigma.proc.count[t] ,-2)
  }
    
    # Priors and constraints FOR POPULATION COUNTS OBSERVATION
    for (s in 1:countries){			### start loop over every country
      sigma.obs.count[s] ~ dunif(0.1,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs.count[s]<-pow(sigma.obs.count[s],-2)
    }
    
    
    
    ##### Likelihood function
    
    ## State process for entire time series
    
    for (t in 1:(T.count-1)){
      for (s in 1:countries){			### start loop over every country
        lambda.count[t,s] ~ dnorm(mean.lambda.count[period[t]], tau.proc.count[period[t]])								# Distribution for random error of growth rate
        N.est.count[t+1,s]<-N.est.count[t,s]*lambda.count[t,s]										# Linear predictor (population size based on past pop size and change rate)
      }
    }
    ## Observation process
    for (t in 1:T.count){
      for (s in 1:countries){			### start loop over every country
        y.count[t,s] ~ dnorm(N.est.count[t,s], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }
    }


    ##### Derived quantity
    for (t in 1:T.count){
      Ntot[t]<-sum(N.est.count[t,])
    }
    
    
  # Assessing the fit of the state-space model
  # 1. Compute fit statistic for observed data
  # Discrepancy measure: mean absolute error
  for (t in 1:T.count){
    C.exp[t] <- Ntot[t]                            # Expected counts
    C.act[t] <- sum(y.count[t,])
    Dssm.obs[t] <- abs((C.act[t] - C.exp[t]) / C.act[t])           # Discrepancy measure
  }
  Dmape.obs <- sum(Dssm.obs)

  # 2. Compute fit statistic for replicate data
  # Discrepancy measure: mean absolute error
  for (t in 1:T.count){
        for (s in 1:countries){			### start loop over every country
         C.rep.c[t,s] ~ dnorm(N.est.count[t,s], tau.obs.count[s])								# Distribution for random error in observed numbers (counts)
      }
    C.rep[t] <- sum(C.rep.c[t,])                 # Generate replicate data
    Dssm.rep[t] <- abs((C.rep[t] - C.exp[t]) / C.rep[t])   # Discrepancy measure
  }
  Dmape.rep <- sum(Dssm.rep)
    
    # -------------------------------------------------        
    # 2. specify parameters for projection
    # -------------------------------------------------
    
    mean.fec ~ dunif(0.914,1.109)
  for (scen in 1:2) {    ### explore two scenarios with pre-LIFE estimates and post-LIFE estimates  
    phi.ad[scen] ~ dbeta(ad.phi.mean[scen],ad.phi.sd[scen])
    phi.juv[scen] ~ dbeta(juv.phi.mean[scen],juv.phi.sd[scen])
    phi.capt.juv[scen] ~ dbeta(sec.capt.phi.mean[scen],sec.capt.phi.sd[scen])
    phi.sec[scen] ~ dbeta(sec.phi.mean[scen],sec.phi.sd[scen])
    phi.third[scen] ~ dbeta(third.phi.mean[scen],third.phi.sd[scen])

    
              
    # -------------------------------------------------        
    # 3. PREDICTION INTO THE FUTURE
    # -------------------------------------------------


        ## POPULATION PROCESS
        ## need to set up population structure
        nestlings.f[1,scen] <- round(Nterr.f[1,scen]*mean.fec * 0.5)   ##JUV[T.count]
        N1.f[1,scen] ~ dbin((mean.fec * 0.5*phi.juv[scen]),round(Ntot[T.count-1]))
        N2.f[1,scen] ~ dbin((mean.fec * 0.5*phi.juv[scen]*phi.sec[scen]),round(Ntot[T.count-2]))
        N3.f[1,scen] ~ dbin((mean.fec * 0.5*phi.juv[scen]*phi.sec[scen]*phi.third[scen]),round(Ntot[T.count-3]))
        N4.f[1,scen] ~ dbin((mean.fec * 0.5*phi.juv[scen]*phi.sec[scen]*phi.third[scen]*phi.ad[scen]),round(Ntot[T.count-4]))
        N5.f[1,scen] ~ dbin((mean.fec * 0.5*phi.juv[scen]*phi.sec[scen]*phi.third[scen]*phi.ad[scen]*phi.ad[scen]),max(1,round(Ntot[T.count-5])))
        # N6.f[1] ~ dbin(phi.ad[scen],round(Ntot[T.count]*0.44+(mean.fec * 0.5*phi.juv[scen]*phi.sec[scen]*phi.third[scen]*phi.ad[scen]*phi.ad[scen]*phi.ad[scen]*Ntot[T.count-6])))                                  ### number of 6-year or older (adult) birds
        # Nterr.f[1] <- round((N4.f[1,scen] * 0.024) + (N5.f[1] * 0.124) + (N6.f[1,scen]))
        N6.f[1,scen] <- round(Nterr.f[1,scen]-(N4.f[1,scen] * 0.024) - (N5.f[1,scen] * 0.124))
        Nterr.f[1,scen] <- Ntot[T.count]
        N2wild.f[1,scen] <- 0
        N2released.f[1,scen] <- 0
        
        for (fut in 2:PROJECTION){
          
          ### probabilistic formulation
          nestlings.f[fut,scen] <- (mean.fec * 0.5 * Nterr.f[fut,scen])           ### number of local recruits calculated as REDUCED fecundity times number of territorial pairs
          N1.f[fut,scen] ~ dbin(phi.juv[scen],round(nestlings.f[fut-1,scen]))             ### +rescued[fut] number of 1-year old survivors
          N2wild.f[fut,scen] ~ dbin(phi.sec[scen],max(1,round(N1.f[fut-1,scen])))                ### number of 2-year old wild survivors
          N2released.f[fut,scen] ~ dbin(phi.capt.juv[scen],round(capt.release[fut]))             ### +rescued[fut] number of 1-year old survivors 
          N2.f[fut,scen] <-  N2wild.f[fut,scen] + N2released.f[fut,scen]       ### sum of the N1 cohort derived from wild and captive-bred juveniles 
          N3.f[fut,scen] ~ dbin(phi.third[scen],max(1,round(N2.f[fut-1,scen])))                                                    ### number of 3-year old survivors
          N4.f[fut,scen] ~ dbin(phi.ad[scen],max(1,round(N3.f[fut-1,scen])))                                                       ### number of 4-year old survivors
          N5.f[fut,scen] ~ dbin(phi.ad[scen],max(1,round(N4.f[fut-1,scen])))                                                       ### number of 5-year old survivors
          N6.f[fut,scen] ~ dbin(phi.ad[scen],max(1,round((N5.f[fut-1,scen]+N6.f[fut-1,scen]))))                                   ### number of 6-year or older (adult) birds
          Nterr.f[fut,scen] <- round((N4.f[fut,scen] * 0.024) + (N5.f[fut,scen] * 0.124) + (N6.f[fut,scen]))
          
        } # fut
        
        for (fut2 in 1:(PROJECTION-1)){
          lambda.t.f[fut2,scen] <- Nterr.f[fut2+1,scen] / max(Nterr.f[fut2,scen],1)
          loglambda.t.f[fut2,scen]<-log(lambda.t.f[fut2,scen])## for calculating geometric mean of overall population growth rate
        } # fut2
        
        
        #### FUTURE POPULATION GROWTH RATE  #########
        fut.lambda[scen]<-exp((1/(PROJECTION-1))*sum(loglambda.t.f[1:(PROJECTION-1),scen]))   # Geometric mean
    } # end of future projection scenarios
      
  } # END of the JAGS model
    ",fill = TRUE)
sink()


EGVU_PVA<- run.jags(data=INPUT, inits=inits.count, monitor=paraPVA,
                    model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_PVA_2023_2scen.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel")

out<-as.data.frame(summary(EGVU_PVA))
out$parameter<-rownames(summary(EGVU_PVA))
## retrieve the past population estimates (2006-2019)
EV.past<-out[(grep("Ntot\\[",out$parameter)),c(12,2,1,3)] %>%
  mutate(Year=seq(2006,(2005+INPUT$T.count),1)) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  #mutate(Scenario=ifelse(Year<2017,"before LIFE","with LIFE")) %>%
  mutate(Scenario="projection without management") %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)

## retrieve the population projections and insert proper scenario labels
EV.fut<-out[(grep("Nterr.f",out$parameter)),c(12,2,1,3)] %>%
  mutate(Scenario=base::rep(x=c("projection without management","projection with management"), each=INPUT$PROJECTION)) %>%
  mutate(Year=base::rep(seq(max(EV.past$Year)+1,(max(EV.past$Year)+INPUT$PROJECTION)), 2)) %>%
  rename(parm=parameter,median=Median,lcl=Lower95,ucl=Upper95) %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)

### FILTER FOR FUTURE PLOT TRAJECTORY
EV.plot <- EV.fut %>% 
  bind_rows(EV.past) %>%
  arrange(Year) %>%
  mutate(ucl=ifelse(ucl>400,400,ucl))

### ADD SUMMARISED COUNT DATA
EV.count<-countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  mutate(prop=ifelse(is.na(N),0,country.props[match(Country,names(country.props))])) %>%
  group_by(year) %>%
  summarise(N=sum(N,na.rm=T),prop=sum(prop)) %>%
  mutate(TOT=N/prop) %>%
  mutate(surv=dplyr::lead(TOT)/TOT) 

ggplot()+
  geom_vline(aes(xintercept=2040),linewidth=1, linetype="dotted")+
  geom_ribbon(data=EV.plot,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  geom_line(data=EV.plot, aes(x=Year, y=median,colour=Scenario),linewidth=1)+
  geom_point(data=EV.count,aes(x=year, y=TOT), size=2,col='darkblue')+
  
  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,400),breaks=seq(0,300,50), labels=as.character(seq(0,300,50)))+
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
        legend.position = c(0.15, 0.92))


