##########################################################################
#
# EGYPTIAN VULTURE POPULATION VIABILITY ANALYSIS - GENERATING OUTPUT
#
##########################################################################
# model is run in EGVU_IPM_PVA_2019.r
# model output saved as workspace and imported here
# script to deliver desired output by HOS

## email by Vic Saravia (HOS) on 29 July 2019 specifies the following graph requirements: 

# 1-A population trajectory of the EV in the Balkans showing the expected year when the population will go extinct 
# 2-A graph like the previous one showing how population numbers would be affected if we take all second wild chicks from the wild population every year for 5 years in order to add them to the captive-bred pool. This graph should not include the effect of releases of birds or increased survival. We want to determine if harvesting juveniles from the wild has a significant negative effect on the population in order to assess if we should use this method or just rely on birds coming from zoos to create our captive bred pool.
# 3-The extinction probability graphs as you have sent them in your previous mail, but without including birds released that have been harvested from the wild.
# 4-I find the violin graph you sent somehow difficult to understand especially for those who are not used to reading that kind of graphs. 
# 5-In general, all graphs (particularly those presenting 4 scenarios to compare) would be easier to read if a grid was added.

## ADDED REQUEST ON 5 Aug 2019 to combine graph 1 and 2
## inconsistency in model output: resolved definition of extinction probability
## changed baseline graph to output from multi-scenario model

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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORT THE MODEL OUTPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
load("EGVU_IPM_output2019_v2.RData")
out<-as.data.frame(NeoIPM.ALL$summary)  ## changed from NeoIPMbasic
out$parameter<-row.names(NeoIPM.ALL$summary) ## changed from NeoIPMbasic
#write.table(out,"EGVU_IPM_estimates_v3.csv", sep=",", row.names=F)

out<-as.data.frame(NeoIPMeggredNoRescue$summary)
out$parameter<-row.names(NeoIPMeggredNoRescue$summary)
#write.table(out,"EGVU_IPM_future_v2.csv", sep=",", row.names=F)

out<-as.data.frame(NeoIPMeggredRescue$summary)
out$parameter<-row.names(NeoIPMeggredRescue$summary)


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
# GRAPH 1: POPULATION TRAJECTORY UNDER THE NO MANAGEMENT SCENARIO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out<-as.data.frame(NeoIPM.ALL$summary)  ## changed from NeoIPMbasic
out$parameter<-row.names(NeoIPM.ALL$summary) ## changed from NeoIPMbasic

## retrieve the population projections
EV.fut<-out[(grep("Nterr",out$parameter)),c(12,5,3,7)] %>%
  mutate(Year=c(trendinput$year,rep(seq(2019,2068,1),each=4*7)))
names(EV.fut)[1:4]<-c('parm','median','lcl','ucl')

## give the projections proper scenario labels
capt.release=c(0,2,4,6)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12)
EV.fut <- EV.fut %>% mutate(capt.release=0, imp.surv=0)
EV.fut$capt.release[grep(",",EV.fut$parm)]<-capt.release[as.numeric(substr(EV.fut$parm[grep(",",EV.fut$parm)],9,9))]
EV.fut$imp.surv[grep(",",EV.fut$parm)]<-imp.surv[as.numeric(substr(EV.fut$parm[grep(",",EV.fut$parm)],11,11))]
EV.fut <- EV.fut %>% filter(capt.release==0 & imp.surv<1.0001)

### produce plot FOR BASELINE TRAJECTORY

#pdf("EV_population_projection_BASELINE.pdf", width=10, height=7)
#jpeg("EV_population_projection_BASELINE.jpg", width=9, height=6, units="in", res=600, quality=100)

ggplot()+
  geom_line(data=EV.fut, aes(x=Year, y=median), color="cornflowerblue",size=1)+
  geom_ribbon(data=EV.fut,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_point(data=trendinput, aes(x=year+0.1, y=N), size=1,col='darkblue')+

  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,130),breaks=seq(0,130,30), labels=as.character(seq(0,130,30)))+
  scale_x_continuous(name="Year", breaks=seq(2006,2068,5), labels=as.character(seq(2006,2068,5)))+
  
  ## ADD LINES FOR EXTINCTION
  geom_vline(xintercept=EV.fut$Year[min(which(EV.fut$lcl<5))],linetype='dashed', size=1,colour="firebrick")+
  geom_vline(xintercept=EV.fut$Year[min(which(EV.fut$median<5))],linetype='dashed', size=1,colour="firebrick")+
  
  ## ADD LABELS FOR EXTINCTION
  geom_text(aes(y=125,x=EV.fut$Year[min(which(EV.fut$lcl<5))],label=paste("5% probability \n in ",xintercept=EV.fut$Year[min(which(EV.fut$lcl<5))])), size=5, colour="firebrick", hjust=1)+
  geom_text(aes(y=125,x=EV.fut$Year[min(which(EV.fut$median<5))],label=paste("50% probability \n in ",xintercept=EV.fut$Year[min(which(EV.fut$median<5))])), size=5, colour="firebrick", hjust=0)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

dev.off()








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 2: POPULATION TRAJECTORY WHEN REMOVING ALL SECOND CHICKS FOR 5 YEARS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out<-as.data.frame(NeoIPMeggredNoRescue$summary)
out$parameter<-row.names(NeoIPMeggredNoRescue$summary)

## retrieve the population projections
EV.fut.chick<-out[(grep("Nterr",out$parameter)),c(12,5,3,7)] %>%
  mutate(Year=c(trendinput$year,rep(seq(2019,2068,1),each=4*7)))
names(EV.fut.chick)[1:4]<-c('parm','median','lcl','ucl')

## give the projections proper scenario labels
capt.release=c(0,2,4,6)
imp.surv=c(1,1.02,1.04,1.06,1.08,1.1,1.12)
EV.fut.chick <- EV.fut.chick %>% mutate(capt.release=0, imp.surv=0)
EV.fut.chick$capt.release[grep(",",EV.fut.chick$parm)]<-capt.release[as.numeric(substr(EV.fut.chick$parm[grep(",",EV.fut.chick$parm)],9,9))]
EV.fut.chick$imp.surv[grep(",",EV.fut.chick$parm)]<-imp.surv[as.numeric(substr(EV.fut.chick$parm[grep(",",EV.fut.chick$parm)],11,11))]
EV.fut.chick <- EV.fut.chick %>% filter(capt.release==0 & imp.surv<1.0001)



### COMBINE BOTH TRAJECTORIES FOR A PLOT OF BOTH SCENARIOS
EV.fut$Scenario<-"Baseline"
EV.fut.chick$Scenario<-"Chick removal"
plotdat<- rbind(EV.fut, EV.fut.chick) 



#pdf("EV_population_projection_TwoScenarios.pdf", width=10, height=7)
#jpeg("EV_population_projection_TwoScenarios.jpg", width=9, height=6, units="in", res=600, quality=100)

ggplot()+
  geom_line(data=plotdat, aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=plotdat,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  geom_point(data=trendinput, aes(x=year+0.1, y=N), size=1,col='darkblue')+
  
  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,130),breaks=seq(0,130,30), labels=as.character(seq(0,130,30)))+
  scale_x_continuous(name="Year", breaks=seq(2006,2068,5), labels=as.character(seq(2006,2068,5)))+
  
  ## ADD LINES FOR EXTINCTION
  geom_vline(xintercept=EV.fut$Year[min(which(EV.fut$lcl<5))],linetype='dashed', size=1,colour="#F8766D")+
  geom_vline(xintercept=EV.fut$Year[min(which(EV.fut$median<5))],linetype='dashed', size=1,colour="#F8766D")+
  geom_vline(xintercept=EV.fut.chick$Year[min(which(EV.fut.chick$lcl<5))],linetype='dotted', size=1,colour="#00BFC4")+
  geom_vline(xintercept=EV.fut.chick$Year[min(which(EV.fut.chick$median<5))],linetype='dotted', size=1,colour="#00BFC4")+
  
  ## ADD LABELS FOR EXTINCTION
  geom_text(aes(y=125,x=EV.fut.chick$Year[min(which(EV.fut.chick$lcl<5))]-1),label="5% probability", size=5, colour="darkblue", hjust=1)+
  geom_text(aes(y=125,x=EV.fut$Year[min(which(EV.fut$median<5))]+1),label="50% probability", size=5, colour="darkblue", hjust=0)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 3: EXTINCTION PROBABILITY OVER TIME WITHOUT RESCUED CHICKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## OUTPUT FROM COMBINED MODEL
## NO CHICK REMOVAL FOR BASELINE SCENARIO
## CHICK REMOVAL FOR ALL SCENARIOS WITH CAPTIVE RELEASES
# replaced NeoIPMeggredNoRescue with NeoIPM.ALL


samplesout<-as.data.frame(rbind(NeoIPM.ALL$samples[[1]],NeoIPM.ALL$samples[[2]],NeoIPM.ALL$samples[[3]],NeoIPM.ALL$samples[[4]]))
head(samplesout)

extprop <- samplesout %>% gather(key="parm", value="value") %>%
  filter(grepl("Nterr.f",parm)) %>%
  mutate(capt.release=as.numeric(substr(parm,9,9)), imp.surv=as.numeric(substr(parm,11,11)))  %>%
  filter(imp.surv %in% c(1,2,3,4)) %>% ### reduce plotting options
  mutate(Year=ifelse(nchar(parm)==14,substr(parm,13,13),substr(parm,13,14))) %>%
  mutate(n=1, inc=ifelse(value<5,1,0)) %>%
  group_by(imp.surv,capt.release,Year) %>%
  summarise(ext.prob=sum(inc)/sum(n)) %>%
  mutate(Year=as.numeric(Year)+2018) 


## create factors for plot labels and order them appropriately
extprop$capt.release <- factor(extprop$capt.release, labels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
extprop$imp.surv <- factor(extprop$imp.surv, labels = c("no improvement","surv +2%", "surv +4%", "surv +6%"))



### produce plot with 4 panels and multiple lines per year

#pdf("EV_extinction_probability_C3.pdf", width=10, height=7)
#jpeg("EV_extinction_probability_C3.jpg", width=9, height=6, units="in", res=600, quality=100)
ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=capt.release), size=1)+
  facet_wrap(~imp.surv,ncol=2) +
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,1),breaks=seq(0,1,0.2), labels=as.character(seq(0,100,20)))+
  scale_x_continuous(name="Year", breaks=seq(2019,2068,5), labels=as.character(seq(2019,2068,5)))+
  guides(color=guide_legend(title="N captive releases"),fill=guide_legend(title="N captive releases"))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_line(colour="darkgrey"),
        panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12))

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPARE WHY EXTINCTION PROBABILITY OF THE GRAPHS 1 and 3 DO NOT MATCH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EV.fut$Year[min(which(EV.fut$median<5))]
extprop %>% filter(imp.surv=="no improvement") %>%
  filter(capt.release=="no captive releases") %>%
  filter(Year>2040)


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










