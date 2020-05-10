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

## CHANGED MODEL IN SEPT 2019 to INCLUDE MORE CHICKS RELEASED AND FEWER SURVIVAL INCREASE SCENARIOS
## significant complication to extract correct output as indices now have either 1 or 2 characters

## MODIFIED ON 14 OCT 2019 after expanding projections

## MODIFIED ON 21 NOV TO INCLUDE 2019 DATA

## email by Vic Saravia (HOS) on 26 NOV 2019 specifies the following graph requirements: 

# -	Would it be possible to add the 50% extinction risk line too? 5% is good because it gives the idea of immediacy, but the 50% transmits better the idea of inevitability.
# -	We would also like to have the previous graph with the 5 year chick removal factor which you had done in the past.


## MODIFIED 20 FEB 2020 TO PRESENT DATA FOR MANUSCRIPT

## COPIED INTO NEW FILE ON 10 MAY 2020 after complete revision of IPM


library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)
filter<-dplyr::filter
select<-dplyr::select



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPORT THE MODEL OUTPUT [INCLUDES MODEL INPUT DATA]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel\\vultures"), silent=T)
#load("EGVU_IPM_output2019_v4.RData")
load("EGVU_IPM2020_output_v1.RData")
out<-as.data.frame(NeoIPM.ALL$summary)  ## changed from NeoIPMbasic
out$parameter<-row.names(NeoIPM.ALL$summary) ## changed from NeoIPMbasic
#write.table(out,"EGVU_IPM_estimates_v3.csv", sep=",", row.names=F)


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

surv.lu<-surv.inc.mat %>% gather(key=Scenario, value=n.released,-Year) %>%
  filter(Year==30) %>%
  separate(Scenario,sep="_",into=c("surv.inc","lag.time")) %>%
  mutate(surv.index=seq_along(Year)) %>%
  select(surv.index,surv.inc,lag.time)

head(ncr.lu)
head(surv.lu)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR REPORT /MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try(setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\PopulationModel"), silent=T)
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.lambda','mean.phi.terrvis[1]','mean.phi.terrvis[2]','ann.phi.juv.telemetry',"ann.phi.capt.rel.first.year",'ann.phi.sec.telemetry','ann.phi.third.telemetry','mu.fec')) %>%
  select(parameter,c(5,3,7))


# ## retrieve the adult survival estimates averaged across all years
# selcol<-grep("mean.phi.terrvis",dimnames(NeoIPM.ALL$samples[[1]])[[2]])
# ann.surv.terrvis<-numeric()
# for (c in 1:nc){
#   ann.surv.terrvis<-c(ann.surv.terrvis,as.numeric(NeoIPM.ALL$samples[[c]][,selcol]))
# }

#TABLE1[6,]<-c("adult survival",quantile(ann.surv.terrvis,0.5),quantile(ann.surv.terrvis,0.025),quantile(ann.surv.terrvis,0.975))
names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("fecundity","delayed release first year survival","wild first year survival","wild second year survival", "wild third year survival","adult survival (before 2015)","adult survival (after 2015)","population growth rate")
fwrite(TABLE1,"EGVU_IPM_demographic_parameter_estimates_IPM2020_v1.csv")





########## CREATE TABLE OF ALL FUTURE POPULATION GROWTH RATES ############################


## retrieve the future lambdas
selcol<-grep("fut.lambda",dimnames(NeoIPM.ALL$samples[[1]])[[2]])
fut.lambda<-NeoIPM.ALL$samples[[1]][,selcol]
for (c in 2:nc){
  fut.lambda<-rbind(fut.lambda,NeoIPM.ALL$samples[[c]][,selcol])
}



FUTLAM<-as.data.frame(fut.lambda) %>% gather(key="parm",value="f.lam") %>%
  group_by(parm) %>%
  summarise(median=quantile(f.lam,0.5),lcl=quantile(f.lam,0.025),ucl=quantile(f.lam,0.975)) %>%
  mutate(capt.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])) %>%
  mutate(surv.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%  
  left_join(ncr.lu, by="capt.index") %>%
  left_join(surv.lu, by="surv.index") %>%
  mutate(surv.inc=ifelse(as.numeric(surv.inc)>1,paste("+",as.integer((as.numeric(surv.inc)-1)*100),"%"),"none")) %>%
  #select(n.rel,n.years,surv.inc,lag.time,median,lcl,ucl) %>%
  arrange(median)

FUTLAM
fwrite(FUTLAM,"EGVU_fut_pop_growth_rate_all_scenarios_IPM2020_v1.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER THE NO MANAGEMENT SCENARIO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
EV.past<-out[(grep("Nterr\\[",out$parameter)),c(12,5,3,7)] %>%
  mutate(Year=trendinput$year)
names(EV.past)[1:4]<-c('parm','median','lcl','ucl')

## retrieve the population projections and insert proper scenario labels
EV.fut<-out[(grep("Nterr.f",out$parameter)),c(12,5,3,7)] %>%
  mutate(Year=rep(seq(max(trendinput$year)+1,(max(trendinput$year)+INPUT$PROJECTION),1),each=dim(ncr.lu)[1]*dim(surv.lu)[1])) %>%
  mutate(capt.index=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])) %>%
  mutate(surv.index=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%  
  left_join(ncr.lu, by="capt.index") %>%
  left_join(surv.lu, by="surv.index") %>%
  rename(parm=parameter,median=`50%`,lcl=`2.5%`,ucl=`97.5%`) %>%
  dplyr::select(parm,n.rel,n.years,surv.inc,lag.time,Year,median,lcl,ucl)


#fwrite(EV.fut,"EGVU_fut_pop_size_all_scenarios.csv")
#fwrite(EV.past,"EGVU_past_pop_size.csv")




### SUMMARISE FOR BASELINE TRAJECTORY
EV.base <- EV.fut %>% filter(n.rel==0 & surv.inc<1.0001 & n.years==10 & lag.time==10) %>%
  select(parm, median, lcl, ucl, Year) %>%
  bind_rows(EV.past) %>%
  arrange(Year)



### CREATE PLOT FOR BASELINE TRAJECTORY

ggplot()+
  geom_line(data=EV.base, aes(x=Year, y=median), color="cornflowerblue",size=1)+
  geom_ribbon(data=EV.base,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_point(data=trendinput, aes(x=year+0.1, y=N), size=1,col='darkblue')+

  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,125),breaks=seq(0,120,30), labels=as.character(seq(0,120,30)))+
  scale_x_continuous(name="Year", breaks=seq(2005,2050,5), labels=as.character(seq(2005,2050,5)))+
  
  ## ADD LINES FOR EXTINCTION
  geom_vline(xintercept=EV.base$Year[min(which(EV.base$lcl<10))],linetype='dashed', size=1,colour="firebrick")+
  geom_vline(xintercept=EV.base$Year[min(which(EV.base$median<10))],linetype='dashed', size=1,colour="firebrick")+
  
  ## ADD LABELS FOR EXTINCTION
  geom_text(aes(y=125,x=EV.base$Year[min(which(EV.base$lcl<10))]-1),label=paste("5% probability \n in ",EV.base$Year[min(which(EV.base$lcl<10))]), size=5, colour="firebrick", hjust=1)+
  geom_text(aes(y=125,x=EV.base$Year[min(which(EV.base$median<10))]-1),label=paste("50% probability \n in ",xintercept=EV.base$Year[min(which(EV.base$median<10))]), size=5, colour="firebrick", hjust=1)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("EGVU_population_projection_BASELINE_IPM2020_v1.jpg", width=9, height=6)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 2: PLOT POPULATION TREND FOR DIFFERENT SCENARIOS OF IMPROVEMENT AND CAPT RELEASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(EV.fut)

## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("cornflowerblue", "firebrick"))


## SELECT ONLY 10 YEAR SURVIVAL LAG
## modify factors for printing

EV.fut %>%
  filter(lag.time==10) %>%
  mutate(surv.inc=ifelse(as.numeric(surv.inc)>1,paste("+",as.integer((as.numeric(surv.inc)-1)*100),"%"),"none")) %>%
  mutate(surv.inc.ord=factor(surv.inc, levels = c("none","+ 2 %","+ 4 %","+ 6 %","+ 8 %","+ 10 %"))) %>%
  mutate(n.rel.ord=factor(n.rel, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))) %>%
  #mutate(lag.time=sprintf("after %s years",lag.time)) %>%
  mutate(n.years=sprintf("for %s years",n.years)) %>%
  arrange(n.rel,Year) %>%
  #mutate(release=paste(n.rel,n.years," ")) %>% 
  #mutate(n.years=as.factor(n.years)) %>%

  
    
### produce plot with 18 panels and multiple lines per year

ggplot()+
  geom_line(aes(x=Year, y=median, color=n.rel.ord),size=1)+
  #geom_ribbon(data=EV.fut,aes(x=Year, ymin=lcl,ymax=ucl, fill=capt.release),alpha=0.2)+
  #geom_line(data=trendinput, aes(x=year, y=N), size=1,col='cornflowerblue')+
  facet_grid(n.years~surv.inc.ord) +
  
  ## format axis ticks
  scale_y_continuous(name="N territorial Egyptian Vultures", limits=c(0,200),breaks=seq(0,200,50), labels=as.character(seq(0,200,50)))+
  scale_x_continuous(name="Year", breaks=seq(2020,2050,5), labels=as.character(seq(2020,2050,5)))+
  guides(color=guide_legend(title="N chicks \n released"))+
  scale_colour_manual(palette=colfunc)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("EGVU_population_projection_allScenarios_IPM2020_v1.jpg", width=16, height=13)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 3: PLOT FUT POP GROWTH AGAINST NUMBER OF RELEASED CHICKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(FUTLAM)

FUTLAM %>% mutate(n.rel=as.numeric(n.rel)) %>%
  filter(lag.time==10) %>%
  #mutate(surv.inc=ifelse(as.numeric(surv.inc)>1,paste("+",as.integer((as.numeric(surv.inc)-1)*100),"%"),"none")) %>%
  mutate(surv.inc.ord=factor(surv.inc, levels = c("none","+ 2 %","+ 4 %","+ 6 %","+ 8 %","+ 10 %"))) %>%
  mutate(n.rel.ord=factor(n.rel, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))) %>%
  #mutate(lag.time=sprintf("after %s years",lag.time)) %>%
  mutate(n.years=sprintf("for %s years",n.years)) %>%
  arrange(n.rel) %>%
  
  ### CREATE PLOT ###
  ggplot()+
  geom_hline(aes(yintercept=1), color='red', size=1)+
  geom_line(aes(x=n.rel,y=median),size=1)+
  geom_ribbon(aes(x=n.rel, ymin=lcl,ymax=ucl),alpha=0.2)+
  facet_grid(n.years~surv.inc.ord) +
  
  #guides(color=guide_legend(title="N years \n releases"),fill=guide_legend(title="N years \n releases"))+
  
  
  ylab("Future population growth rate") +
  xlab("Number of chicks released per year") +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("EGVU_future_growth_rate_allScenarios_IPM2020_v1.jpg", width=10, height=8)









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 4: EXTINCTION PROBABILITY OVER TIME WITHOUT RESCUED CHICKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## OUTPUT FROM COMBINED MODEL
## NO CHICK REMOVAL FOR BASELINE SCENARIO
## CHICK REMOVAL FOR ALL SCENARIOS WITH CAPTIVE RELEASES
# replaced NeoIPMeggredNoRescue with NeoIPM.ALL

### COMMENTED OUT ON 29 NOV 2019 - simply read in csv file

rm(list=setdiff(ls(), c("NeoIPM.ALL","ncr.lu","surv.lu","trendinput")))


### CANNOT PROCESS ALL DATA AT ONCE, BECAUSE MEMORY OVERFLOW. NEED TO LOOP OVER EACH SCENARIO
extprop <- data.frame()

for(scen in 1:nrow(ncr.lu)){

  ### FIND COLUMS WE NEED
  colname<-sprintf("Nterr.f\\[%s,",scen)
  selcol<-grep(colname,dimnames(NeoIPM.ALL$samples[[1]])[[2]])

  allchainsamples <- data.frame()
  for(chain in 1:4) {

      ### EXTRACT AND SUMMARISE DATA
      samplesout<-as.data.frame(NeoIPM.ALL$samples[[1]][,selcol]) %>% gather(key="parm", value="value")
      allchainsamples <- rbind(allchainsamples,as.data.frame(samplesout))
    }

  ### CALCULATE EXTINCTION PROBABILITY
    allchainsamples<- allchainsamples %>%
      mutate(capt.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])) %>%
      mutate(surv.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
      mutate(Year=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,3])+(max(trendinput$year))) %>%

      mutate(n=1, inc=ifelse(value<10,1,0)) %>%
      group_by(capt.index,surv.index,Year) %>%
      summarise(ext.prob=sum(inc)/sum(n))

    extprop <- rbind(extprop,as.data.frame(allchainsamples))
    print(scen)
}

head(samplesout)
head(extprop)
dim(extprop)

fwrite(extprop,"EGVU_ext_prob_all_scenarios_IPM2020_v1.csv")
extprop<-fread("EGVU_ext_prob_all_scenarios_IPM2020_v1.csv")



## create factors for plot labels and order them appropriately

extprop<- extprop %>%
  left_join(ncr.lu, by="capt.index") %>%
  left_join(surv.lu, by="surv.index") %>%
  filter(lag.time==10) %>%
  mutate(surv.inc=ifelse(as.numeric(surv.inc)>1,paste("+",as.integer((as.numeric(surv.inc)-1)*100),"%"),"none")) %>%
  mutate(surv.inc.ord=factor(surv.inc, levels = c("none","+ 2 %","+ 4 %","+ 6 %","+ 8 %","+ 10 %"))) %>%
  mutate(n.rel.ord=factor(n.rel, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))) %>%
  #mutate(lag.time=sprintf("after %s years",lag.time)) %>%
  mutate(n.years=sprintf("for %s years",n.years)) %>%
  arrange(n.rel,Year) %>%
  mutate(release=paste(n.rel,n.years," ")) 
#extprop$capt.release <- factor(extprop$capt.release, labels = c("no captive releases","+ 2 chicks/year","+ 4 chicks/year","+ 6 chicks/year"))
#extprop$imp.surv <- factor(extprop$imp.surv, labels = c("no improvement","surv +2%", "surv +4%", "surv +6%"))
  
dim(extprop)






### produce plot with 4 panels and multiple lines per year

## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("cornflowerblue", "firebrick"))


ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=n.rel.ord), size=1)+
  facet_grid(n.years~surv.inc.ord) +
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,1),breaks=seq(0,1,0.2), labels=as.character(seq(0,100,20)))+
  scale_x_continuous(name="Year", breaks=seq(2020,2050,5), labels=as.character(seq(2020,2050,5)))+
  guides(color=guide_legend(title="N captive releases"),fill=guide_legend(title="N captive releases"))+
  scale_colour_manual(palette=colfunc)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
  
ggsave("EGVU_extinction_probability_allScenarios_IPM2020_v1.jpg", width=9, height=6)



### CREATE TABLE 2 FOR MANUSCRIPT ###
## added 18 Feb 2020
head(extprop)
TABLE2<- extprop %>%
  filter(surv.inc.ord=="none") %>%
  filter(Year==2049) %>%
  select(ext.prob,n.rel.ord,n.years) %>%
  spread(key=n.years,value=ext.prob)
fwrite(TABLE2,"TABLE2.csv")