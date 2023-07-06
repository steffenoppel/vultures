##########################################################################
#
# EGYPTIAN VULTURE ADULT SURVIVAL FROM TERRITORY MONITORING
#
##########################################################################
# based on Oppel et al 2016
# updated to use only 1 occasion per year (rather than monthly survival)
# removed occupancy model formulation (discarded all monthly visits) and switched to simple CJS model formulation
# updated 30 April 2020: switched from using visits data to occupancy table in a binary model (rather than 2/1/0 occupancy)

# updated 5 December 2022: revisited to analyse data from project
## deleted all survival models except the one used for PVA model paper (Oppel et al. 2021)
## trialled various phase cut-offs, but retained pre-LIFE cutoff (including previous LIFE project)

## attempted to include GoF test in June 2023
## GoF test failed - trying to fit in jagsUI to use MCMCvis package
## may need to change model to ensure better fit

### ON 5 JULY I TRIED TO REVISE THE MODEL, AND IF INITS ARE STARTING AT 0.5 NOTHING WORKS

library(readxl)
library(jagsUI)
library(runjags)
library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select
library(jagshelper)
library(MCMCvis)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD AND MANIPULATE TERRITORY MONITORING DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE)
# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend\\ReadEGVUpopdata.r")), wait = TRUE, invisible = FALSE, intern = T)
# try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
# try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationTrend"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Data"), silent=T)
load("EGVU_poptrend_input.RData")
head(occu)

# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_surveys.r")), wait = TRUE, invisible = FALSE)
# system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\RODBC_surveys.r")), wait = TRUE, invisible = FALSE, intern = T)
# 
# try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
# try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
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





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET UP DATA AND INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## scale observation effort
#explored various options to find reasonable scale that avoids invalid parent error
range(effort.terrvis)
eff.scale<-scale(0:1993)

## OPTION 1 - standard scaling function
scaled.effort.matrix<-scale(effort.terrvis)
# dimnames(scaled.effort.matrix) <-NULL
# attr(scaled.effort.matrix,"scaled:scale") <-NULL
# attr(scaled.effort.matrix,"scaled:center") <-NULL
# str(scaled.effort.matrix)

## OPTIONS 2-4 - other scalars
scaled.effort.matrix<-sqrt(effort.terrvis) ## still quite large range
scaled.effort.matrix<-log(effort.terrvis) ## cannot deal with 0 values
scaled.effort.matrix<-effort.terrvis/max(effort.terrvis,na.rm=T)


## Bundle data across all of the data sources

INPUT <- list(y.terrvis = enchist.terrvis,
              nsite.terrvis = nrow(enchist.terrvis),
              nprim.terrvis = ncol(enchist.terrvis),
              #phase=c(1,2,2,1,2,1,1,1,1,2,1,2,1,1),  ## good and bad mixture up to 2020
              phase=c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),  ## test of before and during project
              eff.terrvis=scale(effort.terrvis),
              #firstobs=first.obs,
              f.obsvis=f.obsvis)



## Parameters to be estimated ('monitored') by JAGS
params<-c("mean.phi.terrvis","mean.p.terrvis","mean.rec.terrvis","fit.rep.terr","fit.terr")
#params<-c("y.rep.terr")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SPECIFY INITIAL VALUES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initial values
inits <- function(){list(#z.terrvis = cjs.init.z(enchist.terrvis, f.obsvis),
  mean.p.terrvis=runif(1,0.75, 1),
  mean.phi.terrvis=runif(2,0.75, 1))}  

# MCMC settings
nc <- 3
nt <- 4
ni <- 5000
nb <- 1000
nad<-1000
ns<-10000 	## number of iterations (draws per chain)

### THIS MODEL RUNS IN jagsUI
EV_ad_surv <- jagsUI(data=INPUT,
                       inits=inits,
                       parameters.to.save=params,
                       model.file="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_adult_surv_terrmon_REVISED.jags",
                       n.chains=nc, n.thin=nt, n.burnin=nb, parallel=T,n.iter=ni)
 
### THIS MODEL RUNS IN runjags
# EV_ad_surv <- run.jags(data=INPUT,
#                        inits=inits,
#                        monitor=params,
#                        #model="C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival\\EGVU_adult_surv_terrmon_CJS_recruit.jags",
#                        model="C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/Analysis/EGVU_adult_surv_terrmon_CJS_recruit.jags",
#                        n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
#                        method = "rjparallel")

EV_ad_surv
#extract(EV_ad_surv,what='dic')
#out<-summary(EV_ad_surv[[1]]) %>%   ### for runjags output
out<-as.data.frame(EV_ad_surv$summary)  ## for jagsUI output
out$parameter<-rownames(EV_ad_surv$summary)
fwrite(out,"C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output/EGVU_adult_surv_LIFE.csv")

MCMCdiag(EV_ad_surv$samples,
         round = 3,
         file_name = 'EGVU_adult_surv_LIFE_5Jul2023',
         dir = 'C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output',
         mkdir = 'model-5Jul2023',
         add_field = '1.0',
         add_field_names = 'Data version',
         save_obj = TRUE,
         obj_name = 'model-fit-5Jul2023',
         add_obj = list(INPUT, sessionInfo()),
         add_obj_names = c('model-data-5Jul2023', 'session-info-5Jul2023'))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TEST WITH MCMCvis package
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MCMCplot(EV_ad_surv, params=c("mean.phi.terrvis","mean.rec.terrvis"))

## USING THE CALCULATED FIT VALUES FROM THE JAGS MODEL
OBS <- MCMCpstr(EV_ad_surv, params=c("fit.terr"), type="chains")
REP <- MCMCpstr(EV_ad_surv, params=c("fit.rep.terr"), type="chains")
GOF<-tibble(Rep=as.numeric(REP[[1]]),Obs=as.numeric(OBS[[1]])) %>%
  mutate(P=ifelse(Obs>Rep,1,0))

ggplot(GOF,aes(x=Rep,y=Obs)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  annotate("text",label=as.character(mean(GOF$P)),x=50,y=60)

mean(GOF$P)

  
  
#   
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # GOF TEST WITH mcmc OUTPUT FROM RUNJAGS PACKAGE
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ## TO DO: NEED TO FIGURE OUT WHY y.rep.ter values are >>2 - should not be possible!!
# 
# ## USING THE CALCULATED FIT VALUES FROM THE JAGS MODEL
# # which(dimnames(EV_ad_surv$mcmc[[1]])[[2]]=="E.rep.terr[3,1]")
# # which(dimnames(EV_ad_surv$mcmc[[1]])[[2]]=="E.rep.terr[89,17]")
# which(dimnames(EV_ad_surv$mcmc[[1]])[[2]]=="fit.terr")
# which(dimnames(EV_ad_surv$mcmc[[1]])[[2]]=="fit.rep.terr")
# # OBS <- apply(EV_ad_surv$mcmc[[1]][,6:1153],1,sum,na.rm=T)
# # REP <- apply(EV_ad_surv$mcmc[[1]][,1154:2301],1,sum,na.rm=T)
# OBS <- c(EV_ad_surv$mcmc[[1]][,6],EV_ad_surv$mcmc[[2]][,6],EV_ad_surv$mcmc[[3]][,6])
# REP <- c(EV_ad_surv$mcmc[[1]][,7],EV_ad_surv$mcmc[[2]][,7],EV_ad_surv$mcmc[[3]][,7])
# tibble(Rep=REP,obs=OBS) %>%
#   ggplot(aes(x=REP,y=OBS)) + geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# 
# round(mean(REP > OBS),2)
# 
# 
# ## CALCULATING THE FIT STATISTIC MANUALLY IN R
# ### EXTRACT THE GOF replicated y matrix to compare to input matrix
# dim(enchist.terrvis)
# 
# jags_plist<-function (x, p = NULL) 
# {
#   x_dflist <- lapply(x$mcmc, as.data.frame)
#   x2 <- lapply(1:length(x_dflist[[1]]), function(x) sapply(x_dflist, 
#                                                            "[[", x))
#   names(x2) <- names(x_dflist[[1]])
#   these <- rep(F, length(x2))
#   if (!is.null(p)) {
#     for (i in 1:length(p)) {
#       these[substr(names(x2), 1, nchar(p[i])) == p[i]] <- T
#     }
#     x2 <- x2[these]
#   }
#   if (length(x2) == 0) 
#     warning("No parameters with matching names, returning empty list")
#   return(x2)
# }
# 
# out_plist <- jags_plist(x=EV_ad_surv, p="y.rep.terr")
# exp_plist <- jags_plist(x=EV_ad_surv, p="exp.terr")
# 
# ## trying to recreate matrix
# arraylocs<-tibble(parm=names(out_plist)) %>%
#   mutate(row=sapply(strsplit(names(out_plist), ","), "[", 1)) %>%
#   mutate(col=sapply(strsplit(names(out_plist), ","), "[", 2)) %>%
#   mutate(col=as.numeric(gsub("\\D", "", col))) %>%
#   mutate(row=as.numeric(gsub("\\D", "", row))) 
# 
# 
# ### loop over chains and samples
# OBS<-enchist.terrvis
# EXP_array<-list()
# REP_array<-list()
# fit.obs<-list()
# fit.rep<-list()
# for(ch in c(1,2,3)) {
#   for(it in 1:1000) {
#     s=ifelse(ch==1,it,ifelse(ch==2,it+1000,it+2000))
#     EXP_array[[s]]<-matrix(0, nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis))
#     REP_array[[s]]<-matrix(0, nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis))
#     fit.obs[[s]]<-matrix(0, nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis))
#     fit.rep[[s]]<-matrix(0, nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis))
#     for(n in 1:length(arraylocs$parm)){
#       EXP_array[[s]][arraylocs$row[n],arraylocs$col[n]]<-ifelse(exp_plist[[n]][it,ch]>2,2,exp_plist[[n]][it,ch])
#       REP_array[[s]][arraylocs$row[n],arraylocs$col[n]]<-ifelse(out_plist[[n]][it,ch]>2,2,out_plist[[n]][it,ch])
#       
#       ### CALCULATE DISCREPANCY MEASURE
#       fit.obs[[s]][arraylocs$row[n],arraylocs$col[n]]<-((OBS[arraylocs$row[n],arraylocs$col[n]] - EXP_array[[s]][arraylocs$row[n],arraylocs$col[n]])^2) / (EXP_array[[s]][arraylocs$row[n],arraylocs$col[n]] + 0.5)
#       fit.rep[[s]][arraylocs$row[n],arraylocs$col[n]]<-((REP_array[[s]][arraylocs$row[n],arraylocs$col[n]] - EXP_array[[s]][arraylocs$row[n],arraylocs$col[n]])^2) / (EXP_array[[s]][arraylocs$row[n],arraylocs$col[n]] + 0.5)
#       
#     }
#     
#   }
#   
# }
# 
# 
# 
# ### summarise FIT STATISTICS AND CREATE PLOT
# OBS <- sapply(fit.obs,sum, na.rm=T)
# REP <- sapply(fit.rep,sum, na.rm=T)
# tibble(Rep=REP,obs=OBS) %>%
#   ggplot(aes(x=REP,y=OBS)) + geom_point() +
#   geom_abline(intercept = 0, slope = 1)
# 
# PVAL=round(mean(REP > OBS),2)
# 
# 
# ### trying to understand REP array and why the values are so high
# # s=3000
# # hist(fit.rep[[s]])
# # hist(fit.obs[[s]])
# # REP_array[[s]][which(fit.rep[[s]]>10)]
# # EXP_array[[s]][which(fit.rep[[s]]>10)]
# # OBS[which(fit.rep[[s]]>10)]
# # out_plist[which(sapply(out_plist, FUN=function(X) X>100))]
# # out_plist[[1100]]
# 
# ### trying to populate matrix from maximum output value (across chains and sims)
# matrix(unlist(lapply(out_plist,max)), nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis), byrow=F)
# REP_mat<-EXP_mat<-matrix(0, nrow=nrow(enchist.terrvis),ncol=ncol(enchist.terrvis), byrow=F)
# for(n in 1:length(arraylocs$parm)){
#   #EXP_mat[arraylocs$row[n],arraylocs$col[n]]<-min(exp_plist[[n]])
#   REP_mat[arraylocs$row[n],arraylocs$col[n]]<-min(out_plist[[n]])
# }
# 
# REP_mat[which(REP_mat>2,arr.ind = TRUE)[,1],]
# EXP_mat[which(REP_mat>2,arr.ind = TRUE)[,1],]
# enchist.terrvis[which(REP_mat>2,arr.ind = TRUE)[,1],]
# 
# 
# 
# ####### ALTERNATIVE GOF TEST ####################################################
# ### CODE FROM PHIL PATTON ##
# # https://github.com/philpatton/cowbird/tree/main/R
# # based on paper https://link.springer.com/article/10.1007/s10530-022-02825-3
# 
# ### essentially just compare raw observed (enchist.terrvis) and replicated observations (y.rep.terr)
# GOF_test<-tibble(REP=as.numeric(REP_mat),OBS=as.numeric(enchist.terrvis)) %>%
# 		filter(!is.na(REP)) %>%
# 		filter(!is.na(OBS)) %>%
# 		mutate(P=ifelse(REP>OBS,1,0))
# mean(GOF_test$P)
# GOF_test %>% gather(key="Data",value="value",-P) %>%
#   ggplot(aes(x = value, fill = Data)) +                       # Draw overlaying histogram
#   geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = ..density..), color="black") +
#   geom_density(alpha=0.5)
# 
# 





######################################################################################
# FIGURE TO SHOW POSTERIOR PROBABILITY DISTRIBUTION
######################################################################################

## get all mcmc samples for fut.growth.rate[1-3]
which(dimnames(EV_ad_surv$samples[[1]])[[2]]=="mean.phi.terrvis[1]")
which(dimnames(EV_ad_surv$samples[[1]])[[2]]=="mean.phi.terrvis[2]")


## collate all samples
phi.samp<-data.frame()
for(ch in 1:nc){
  phi.samp<-bind_rows(phi.samp,as.data.frame((EV_ad_surv$samples[[ch]])[,1:2]))
}
head(phi.samp)
dim(phi.samp)

## PROBABILITY OF FUTURE GROWTH
phi.samp %>% gather(key="Period",value="lam") %>% mutate(pos=ifelse(lam>1,1,0)) %>%
  group_by(Period) %>%
  summarise(prob=sum(pos)/dim(phi.samp)[1])

## plot histograms for future pop growth rate
phi.samp %>% rename(before=`mean.phi.terrvis[1]`,during=`mean.phi.terrvis[2]`) %>%
  gather(key="Period",value="N") %>% 
  mutate(Period=if_else(Period=="before", "before LIFE projects","since LIFE projects started")) %>%
  
  ggplot(aes(x = N, fill = Period)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = after_stat(density)), color="black") +
  geom_density(alpha=0.5) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  #geom_segment(aes(x = 1, y = 0, xend = 1, yend = 30), colour="gray15", linetype = "dashed", size=1)+
  
  labs(x="Annual adult survival probability", y="Probability density",
       fill="Period") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=20, color="black"), 
        axis.title=element_text(size=24),
        legend.text=element_text(size=20),
        legend.title = element_text(size=22),
        legend.position=c(0.2,0.88),
        panel.grid.major = element_line(linewidth=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output/EGVU_survival_comparison_LIFE.jpg", width=11, height=8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE PROPORTIONAL INCREASE IN SUVIVAL FOR FINAL REPORT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## CALCULATING DIFFERENCES FOR EACH SAMPLE
((phi.samp %>% rename(before=`mean.phi.terrvis[1]`,during=`mean.phi.terrvis[2]`) %>%
  mutate(diff=during-before) %>%
  summarise(mean=mean(diff), lcl=quantile(diff, 0.025), ucl=quantile(diff, 0.975))) /mean(phi.samp[,1]))*100

## CALCULATING DIFFERENCES FOR EACH SAMPLE
(phi.samp %>% rename(before=`mean.phi.terrvis[1]`,during=`mean.phi.terrvis[2]`) %>%
  mutate(increase=ifelse(during>before,1,0)) %>%
  summarise(sum(increase))) /dim(phi.samp)[1]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT / SAVE FINAL OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival")
try(setwd("C:\\Users\\steffenoppel\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Survival"), silent=T)
try(setwd("C:/Users/sop/Documents/Steffen/RSPB/EGVU/EGVU_LIFE/output"), silent=T)
save.image("EGVU_LIFE_territory_survival_output_2022.RData")
load("EGVU_LIFE_territory_survival_output_2022.RData")




