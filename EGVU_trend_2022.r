##########################################################################
#
# EGYPTIAN VULTURE POPULATION TREND AND BREEDING SUCCESS
#
##########################################################################
## based on Velevski et al. 2015
## SSM model based on Kery and Schaub 2012, Ch 5
## simple Poisson model for productivity, test for ingredients of IPM

## UPDATED on 6 April 2018 to streamline data prep and prepare for IPM
## loads data via a 32-bit R script "ReadEGVUpopdata.r"

## updated 16 May 2018 to prepare for IPM - re-labelled all input features and nodes in the model to avoid duplication with count and reproduction models
## updated 8 June to find priors for IPM - included data from ALB and NMD

## modified on Sept 7 2021 to prepare abstract for EBCC conference
## updated 23 Feb 2022 to include results for conference and build a broken-stick model

### REVISED IN JUNE 2023 to include GOF TEST FOR SSM and implement the changes from Animal Conservation reviewer
### removed future scenarios to just have single future projection

library(data.table)
library(jagsUI)
library(runjags)
library(tidyverse)
library(ggplot2)
filter<-dplyr::filter
select<-dplyr::select


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE POPULATION MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
# 
# 
# try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
# try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
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


### ADD DATA FROM 2022
#countrytrendinput<- fread("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\SMARTbirds\\EGVU_D1_Terr_Occupancy_2022.csv") %>%
countrytrendinput<- fread("EGVU_D1_Terr_Occupancy_2022.csv") %>%
  mutate(year=2022) %>%
  group_by(year,Country) %>%
  summarise(N=sum(N_adults)) %>%
  spread(key=Country,value=N) %>%
  rename(Albania=ALB,Bulgaria=BG,Greece=GR,`North Macedonia`=NMD) %>%
  bind_rows(countrytrendinput) %>%
  arrange(year)
countrytrendinput


## calculate mean proportion of each country to total count
country.props<-colMeans(prop.table(as.matrix(countrytrendinput[,-1]), margin = 1), na.rm=T)


### PLOT RAW COUNT DATA 

countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  ggplot()+
  geom_line(aes(x=year, y=N, col=Country), size=1)

### INSPECT TOTAL

totcount<-countrytrendinput %>% gather(key="Country", value="N",-year) %>%
  mutate(prop=ifelse(is.na(N),0,country.props[match(Country,names(country.props))])) %>%
  group_by(year) %>%
  summarise(N=sum(N,na.rm=T),prop=sum(prop)) %>%
  mutate(TOT=N/prop) %>%
  mutate(surv=dplyr::lead(TOT)/TOT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARIES FOR MANUSCRIPT 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max(totcount$N)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODIFY THE DATA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


trendinput<- occu %>% filter(year>2005) %>%
  filter(Country %in% c("Bulgaria","Greece")) %>%    # introduced in 2019 because database now has data from albania and Macedonia
  group_by(year) %>%
  summarise(N=sum(occupancy), R=sum(breeding, na.rm=T), J=sum(fledglings, na.rm=T))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY TREND MODEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### MODEL WITH TWO SLOPES AND COUNTRY-SPECIFIC ABUNDANCE
sink("EGVU_Balkan_trend2_country.jags")
cat("
    model{
    
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
    
    
}														# close the model loop
    ", fill=TRUE)
sink()
###############  End of model in BUGS language ############







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DATA AND INITS AND RUN MODEL  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Bundle data
# INPUT.count<-list(y.count=trendinput$N, T.count=length(trendinput$N), period=ifelse(trendinput$year<2016,1,2))		## year is standardized so that covariate values are not too far away from zero
# 
# ## Initial values
# inits.count<- function() {list(sigma.proc.count=runif(1,0,5),
#                          mean.lambda.count=runif(1,0.1,2),
#                          sigma.obs.count=runif(1,0,10),
#                          N.est.count=c(runif(1,90,150),rep(NA,(length(trendinput$N)-1))))}


## Parameters to be estimated ('monitored') by WinBUGS
params.count<-c("mean.lambda.count","Ntot","Dmape.rep","Dmape.obs")


## MCMC settings
ns<-10000  	## number of iterations (draws per chain)
nt<-25		  ## thinning rate
nb<-5000	  ## length of burn-in (number of iterations discarded at the start of each chain until convergence is reached)
nc<-3		    ## number of chains
nad<-1000

## call JAGS from R   ###
## label model with name of country:
# EGVUtrend<- jagsUI(data=INPUT.count, inits=inits.count, parameters.to.save=params.count,
#                    model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\EGVU_trend.jags",
#                    n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)



#### RUN MODEL FOR ALL 4 COUNTRIES #####

INPUT.count<-list(y.count=as.matrix(countrytrendinput[,2:5]), T.count=dim(countrytrendinput)[1],
                  country.prop=as.numeric(country.props),  ## proportion of count data that come from countries ALB, BG, GR, NMD
                  countries=4, period=ifelse(trendinput$year<2018,1,2))
# inits.count<- function() {list(sigma.proc.count=runif(1,0.1,5),
#                                mean.lambda.count=runif(1,0.1,2),
#                                sigma.obs.count=runif(4,0.1,100),
#                                N.est.count=c(runif(1,190,250),rep(NA,(length(trendinput$N)-1))))}
# 
# EGVUtrend<- autojags(data=INPUT.count, inits=inits.count, parameters.to.save=params.count,
#                    model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\EGVU_Balkan_trend.jags",
#                    n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)
# 
# inits.count2<- function() {list(sigma.proc.count=runif(2,0.1,5),
#                                mean.lambda.count=runif(2,0.1,2),
#                                sigma.obs.count=runif(4,0.1,100),
#                                N.est.count=c(runif(1,190,250),rep(NA,(length(trendinput$N)-1))))}
# 
# EGVUtrend2<- autojags(data=INPUT.count, inits=inits.count2, parameters.to.save=params.count,
#                      model.file="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\EGVU_Balkan_trend2.jags",
#                      n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T)

inits.count3<- function() {list(sigma.proc.count=runif(2,0.1,5),
                                mean.lambda.count=runif(2,0.1,2),
                                sigma.obs.count=runif(4,0.1,100))}

# EGVUtrend3<- run.jags(data=INPUT.count, inits=inits.count3, monitor=params.count,
#                       model="C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\EGVU_Balkan_trend2_country.jags",
#                       n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
#                       method = "rjparallel") 

try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis"), silent=T)
EGVUtrend3<- run.jags(data=INPUT.count, inits=inits.count3, monitor=params.count,
                      model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_Balkan_trend2_country.jags",
                      n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                      method = "rjparallel") 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TEST  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EV.gof<-as.data.frame(summary(EGVUtrend3, vars=c("Dmape.rep","Dmape.obs")))
EV.gof$parameter<-as.character(row.names(EV.gof))

OBS <- EGVUtrend3$mcmc[[1]][,22]
REP <- EGVUtrend3$mcmc[[1]][,21]
tibble(Rep=REP,obs=OBS) %>%
  ggplot(aes(x=REP,y=OBS)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

PVAL=round(mean(REP > OBS),2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT OUTPUT  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EV.ssm<-as.data.frame(summary(EGVUtrend3, vars='Ntot'))
EV.ssm$parameter<-as.character(row.names(EV.ssm))


EV.ssm<-EV.ssm[(grep("Ntot",EV.ssm$parameter)),c(1,2,3)] %>%
  mutate(Year=seq(2006,2022,1))
names(EV.ssm)[c(2,1,3)]<-c('mean','lcl','ucl')


## DRAW FIGURE

ggplot()+
  geom_line(data=EV.ssm, aes(x=Year, y=mean), size=3,col='darkred')+
  geom_ribbon(data=EV.ssm,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.1, fill='darkred')+
  geom_point(data=totcount, aes(x=year,y=TOT),col="cornflowerblue",size=4) +

  ## format axis ticks
  scale_y_continuous(name="Number of adult territorial Egyptian Vultures", limits=c(75,225),breaks=seq(80,220,20), labels=as.character(seq(80,220,20)))+
  scale_x_continuous(name="Year", breaks=seq(2006,2022,2), labels=as.character(seq(2006,2022,2)))+
  geom_segment(aes(x = 2017.5, y = 80, xend = 2017.5, yend = 220), colour="gray15", linetype = "dashed", size=1)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=20, color="black"), 
        axis.title=element_text(size=24), 
        strip.text.x=element_text(size=20, color="black"), 
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"),
        strip.background=element_rect(fill="white", colour="black"))

ggsave("EGVU_pop_trend_Balkans_2022.jpg", width=10, height=8)



######################################################################################
# FIGURE 3 CALCULATE PROBABILITY OF LAMBDA <1
######################################################################################

## get all mcmc samples for fut.growth.rate[1-3]
which(dimnames(EGVUtrend3$mcmc[[1]])[[2]]=="mean.lambda.count[1]")
which(dimnames(EGVUtrend3$mcmc[[1]])[[2]]=="mean.lambda.count[2]")


## collate all samples
fut.lam.samp<-data.frame()
for(ch in 1:nc){
  fut.lam.samp<-bind_rows(fut.lam.samp,as.data.frame((EGVUtrend3$mcmc[[ch]])[,1:2]))
}
head(fut.lam.samp)
dim(fut.lam.samp)

## PROBABILITY OF FUTURE GROWTH
fut.lam.samp %>% gather(key="Period",value="lam") %>% mutate(pos=ifelse(lam>1,1,0)) %>%
  group_by(Period) %>%
  summarise(prob=sum(pos)/dim(fut.lam.samp)[1])

## plot histograms for future pop growth rate
fut.lam.samp %>% rename(before=`mean.lambda.count[1]`,during=`mean.lambda.count[2]`) %>%
  gather(key="Period",value="N") %>% 
  mutate(Period=if_else(Period=="before", "before flyway project","since flyway project started")) %>%
  
  ggplot(aes(x = N, fill = Period)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = ..density..), color="black") +
  geom_density(alpha=0.5) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 30), colour="gray15", linetype = "dashed", size=1)+
  
  labs(x="Population growth rate", y="Probability density",
       fill="Period") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=20, color="black"), 
        axis.title=element_text(size=24),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("EGVU_pop_growth_rate_comparison.jpg", width=14, height=8)



######################################################################################
# CALCULATE INCREASE IN LAMBDA FOR REPORT
######################################################################################

((summary(EGVUtrend3)[2,2]-summary(EGVUtrend3)[1,2])/summary(EGVUtrend3)[1,2]) *100

EV.ssm<-summary(EGVUtrend3)
fwrite(EV.ssm,"Balkans_EGVU_poptrend_output2022.csv")