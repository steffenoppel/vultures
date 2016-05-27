##########################################################################
# MIGRATION DATA FOR EGYPTIAN VULTURES
# Oppel, S., Dobrev, V., Arkumarev, V., Saravia, V., Bounas, A., Kret, E., Velevski, M., Stoychev, S. & Nikolov, S.C. 2015. High juvenile mortality during migration in a declining population of a long-distance migratory raptor. Ibis, 157(3): 545-557.
# paper available at http://lifeneophron.eu/files/docs/1455988483_251.pdf
# original script written by Steffen Oppel
# August 2013
# last modification 5 July 2014: switched from Transmitter_ID to Bird_ID as identifier, because one transmitter was used on two separate birds in 2010 and 2011
# 7 July 2014: gave up on the dynBBMM approach, as no sensible output maps; reverted to kernel analysis at population level
# 23 September 2014: outsourced all data preparation to script EGVU_migration_data_prep.r
# 29 Sept 2014: replaced object 'tracks' with 'migration'
# LAST UPDATE 22 Jan 2015: removed body mass as not specifically requested by reviewers
##########################################################################

# Load necessary library
#library(BBMM)
library(maptools)
library(sp)
library(rgdal)
library(adehabitatHR)
require(maps)
require(mapdata)
require(adehabitatLT)
require(maptools)
require(geosphere)
require(rgeos)
library(raster)
library(move)
library(ggmap)
library(colorspace)
library(lattice)
library(plotKML)
library(ggplot2)
library(adehabitatHR)
library(randomForest)
library(party)
library(verification)
library(reshape)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PREPARED WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")

load("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_input2.RData")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_input2.RData")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE START AND END DATES OF FALL MIGRATION FOR INDIVIDUAL ANIMALS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anim<-unique(birds$Name)
mig_summary<-birds[!birds$Transmitter_ID %in% c(118149,130675),c(1:13,18)]
mig_summary$start_mig<-NA
mig_summary$end_mig<-NA
mig_summary$mig_days<-0
mig_summary$mig_dist<-0
mig_summary$mig_speed<-0
mig_summary$mig_completed<-0


### PLOT OF CUMULATIVE DISTANCE TRAVELLED #######

par (mfrow=c(6,3), mar=c(3,4,1,0), oma=c(0,0,0,0))
for (a in anim){
input<-migration[migration$Bird_ID == a,]
plottitle<-paste(a,birds$Transmitter_ID[match(a,birds$Name)], sep=" ")
plot(cumul_dist~Date, type='l', frame=F, data=input, main=plottitle)
}


### LOOP TO CALCULATE START AND END DATES #######

all_migdata<-data.frame()

for (a in mig_summary$Name){
input<-migration[migration$Bird_ID == a,]
dailydist<-aggregate(step_dist~Date, data=input, FUN=sum)
start<-as.Date(min(dailydist$Date[dailydist$step_dist>50]), format="%Y-%m-%d")+1		## migration starts when bird travels >50 km per day
#subSaharaArr<-min(input$Date[input$lat<22])
subSahara<-input[input$lat<22,]
if (dim(subSahara)[1]>0){									## calculate end of migration only for birds that crossed the Sahara
southward<-aggregate(home_dist~Date, data=subSahara, FUN=mean)			## calculate average distance to home for each day
southward$reversal<-0
mig_summary$mig_completed[mig_summary$Name==a]<-1
for (d in 2:dim(southward)[1]){
southward$reversal[d]<-southward$home_dist[d]-southward$home_dist[d-1]		## calculate difference in 
}
end<-as.Date(min(southward$Date[southward$reversal<0], na.rm=T), format="%Y-%m-%d") ## migration ends when dist to home no longer increases after the bird is south of 22 N
}else{end<-as.Date(max(input$Date), format="%Y-%m-%d")}				## for birds that died on migration the migration ends when the bird is dead
mig_summary$mig_days[mig_summary$Name==a]<-end-start

migdata<-input[as.Date(input$Date)>=start,]
migdata<-migdata[as.Date(migdata$Date)<=end,]
all_migdata<-rbind(all_migdata, migdata)
mig_summary$mig_dist[mig_summary$Name==a]<-sum(migdata$step_dist)
mig_summary$mig_speed[mig_summary$Name==a]<-mean(migdata$speed)
mig_summary$start_mig[mig_summary$Name==a]<-as.character(start)
mig_summary$end_mig[mig_summary$Name==a]<-as.character(end)

}

mig_summary$travel_speed<-mig_summary$mig_dist/mig_summary$mig_days


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DETERMINE START DATES OF SPRING MIGRATION FOR INDIVIDUALS MANUALLY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mig_summary$spring_mig<-as.Date(mig_summary$Tag_date, format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==108050]<-as.Date("2013-03-08", format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==118153]<-as.Date("2014-05-01", format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==118152]<-as.Date("2014-05-24", format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==118148]<-as.Date("2014-05-20", format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==118150]<-as.Date("2014-05-08", format="%Y-%m-%d")
mig_summary$spring_mig[mig_summary$Transmitter_ID==118146]<-as.Date("2014-05-15", format="%Y-%m-%d")
nospringmig<-ifelse(mig_summary$spring_mig<as.Date(mig_summary$start_mig, format="%Y-%m-%d"),0,1)
mig_summary$spring_mig[nospringmig==0]<-NA





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ASSIGN SEASONS TO EACH LOCATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
migration$season<-NA

for (a in mig_summary$Name){
start<-as.POSIXct(mig_summary$start_mig[mig_summary$Name==a], format="%Y-%m-%d")
migration$season[migration$Bird_ID==a & migration$Date < start]<-"premig"
migration$season[migration$Bird_ID==a & migration$Date >= start]<-"aut.mig"
end<-as.POSIXct(mig_summary$end_mig[mig_summary$Name==a], format="%Y-%m-%d")
migration$season[migration$Bird_ID==a & migration$Date >= end]<-"winter"
spring<-as.POSIXct(mig_summary$spring_mig[mig_summary$Name==a], format="%Y-%m-%d")
migration$season[migration$Bird_ID==a & migration$Date >= spring]<-"spring.mig"
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE MIGRATION DATA FOR SUCCESSFUL AND UNSUCCESSFUL MIGRANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mig_table<-aggregate(strptime(mig_summary$start_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=mean)
mig_table$start_sd<-aggregate(strptime(mig_summary$start_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$start_min<-aggregate(strptime(mig_summary$start_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$start_max<-aggregate(strptime(mig_summary$start_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$end_mig<-aggregate(strptime(mig_summary$end_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$end_sd<-aggregate(strptime(mig_summary$end_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$end_min<-aggregate(strptime(mig_summary$end_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$end_max<-aggregate(strptime(mig_summary$end_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$mig_days<-aggregate(mig_days~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$days_sd<-aggregate(mig_days~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$days_min<-aggregate(mig_days~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$days_max<-aggregate(mig_days~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$speed<-aggregate(travel_speed~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$speed_sd<-aggregate(travel_speed~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$speed_min<-aggregate(travel_speed~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$speed_max<-aggregate(travel_speed~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$spring_mig<-aggregate(strptime(mig_summary$spring_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$spring_sd<-aggregate(strptime(mig_summary$spring_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$spring_min<-aggregate(strptime(mig_summary$spring_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$spring_max<-aggregate(strptime(mig_summary$spring_mig,"%Y-%m-%d")$yday~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$distance<-aggregate(mig_dist~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$distance_sd<-aggregate(mig_dist~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$distance_min<-aggregate(mig_dist~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$distance_max<-aggregate(mig_dist~mig_completed+Age, data=mig_summary, FUN=max)[,3]

mig_table$long<-aggregate(Tag_long~mig_completed+Age, data=mig_summary, FUN=mean)[,3]
mig_table$long_sd<-aggregate(Tag_long~mig_completed+Age, data=mig_summary, FUN=sd)[,3]
mig_table$long_min<-aggregate(Tag_long~mig_completed+Age, data=mig_summary, FUN=min)[,3]
mig_table$long_max<-aggregate(Tag_long~mig_completed+Age, data=mig_summary, FUN=max)[,3]



dummy<-as.Date("2013-12-31")
for (col in c(3,5,6,7,9,10,19,21,22)){
mig_table[,col]<-mig_table[,col]+dummy
}

mig_table[,c(3,5,6,7,9,10,19,21,22)]<-format(mig_table[,c(3,5,6,7,9,10,19,21,22)], format="%d %b")

mig_table

write.table(mig_table,"A:\\MANUSCRIPTS\\in_prep\\EGVU_migration\\Table2.csv", sep=",", row.names=F)
write.table(mig_table,"C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\EGVU_migration\\Table2.csv", sep=",", row.names=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANALYSE LONGEST DAILY TRAVEL DISTANCES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dailydist<-aggregate(step_dist~Bird_ID+Date, data=migration, FUN=sum)
dailydist<-dailydist[order(dailydist$step_dist, decreasing=T),]
head(dailydist)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANALYSE WHAT VARIABLES DIFFERENTIATE BETWEEN SUCCESSFUL AND UNSUCCESSFUL MIGRANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### DEFINE DATA SET AS 10 days surrounding departure

decis<-read.table("Mig_decision_dates.csv", header=T, sep=",")
decis$start_mig<-as.POSIXct(decis$start_mig, format="%d/%m/%Y")
decis$mig_decision<-as.POSIXct(decis$mig_decision, format="%d/%m/%Y")

migdep<-data.frame()

for (a in decis$Bird_ID){
start<-decis$start_mig[decis$Bird_ID==a]
end<-decis$mig_decision[decis$Bird_ID==a]
x<-migration[migration$Bird_ID==a & migration$Date >= start,]
x<-x[x$Bird_ID==a & x$Date <= end,]
migdep<-rbind(migdep,x)
}

migdep$success<-mig_summary$mig_completed[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]
#migdep$dep.date<-strptime(mig_summary$start_mig[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)],"%Y-%m-%d")$yday				### single departure day
migdep$dep.date<-as.numeric(format(as.POSIXct(mig_summary$start_mig[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]),format="%U"))	### departure week to avoid overfitting by RF
migdep$dep.long<-mig_summary$Tag_long[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]									### single departure longitude
migdep$dep.lat<-mig_summary$Tag_lat[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]									### single departure longitude
migdep$mass<-mig_summary$Mass[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]	
migdep$head<-mig_summary$Head[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]	
migdep$tail<-mig_summary$Tail[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]	
migdep$tarsus<-mig_summary$Tarsus[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]	
migdep$jday<-strptime(migdep$DateTime,"%Y-%m-%d")$yday
migdep$fledge.date<-as.numeric(format(as.POSIXct(mig_summary$Fledge_date[match(migdep$Transmitter_ID, mig_summary$Transmitter_ID)]),format="%U"))	### departure week to avoid overfitting by RF
head(migdep)

### ANALYSIS AT LOCATION LEVEL ###
## previous analysis included mass - this was excluded on 22 Jan because measurements are not representative of 

RFmig<-randomForest(success~speed+fledge.date+dep.date+dep.long+dep.lat+winddir+windspeed+orog.uplift+therm.uplift+pressure+Cloud+tail+head+tarsus, data=migdep, mtry=8, ntree=1500, replace=F, importance=T, na.action=na.omit)
varImpPlot(RFmig)
RFmig


### PLOT VARIABLE IMPORTANCE ####

IMP<-as.data.frame(importance(RFmig)) 			
names(IMP)[1]<-"MSE"
IMP<-IMP[order(IMP$MSE, decreasing=T),]  			## SORTED BY MSE
IMP$rel_imp<-round((IMP$MSE/IMP$MSE[1])*100,2)

pdf("EGVU_mig_success_variable_importance.pdf", width=10, height=7)
par(mar=c(5,7,2,1))
barplot(IMP$rel_imp[8:1], horiz=T, names.arg=row.names(IMP)[8:1], xlim=c(0,100), las=1,xlab="Relative importance (%)", col='darkred',main="Importance of variables for successful migration")
text(100,1,sprintf("Var explained = %s",round(RFmig$rsq[1500],3)),cex=1.5, adj=1)
dev.off()



### PLOT ENVIRONMENT RELATIONSHIPS ####

postscript(file = "A:\\MANUSCRIPTS\\submitted\\EGVU_migration\\Figure2.eps", width = 9, height = 8)

#jpeg(filename = "A:\\MANUSCRIPTS\\submitted\\EGVU_migration\\Figure4mig.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
#jpeg(filename = "C:\\STEFFEN\\MANUSCRIPTS\\submitted\\EGVU_migration\\Figure4mig.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
#pdf("EGVU_mig_variables.pdf", width=10, height=10)

par(mfrow=c(2,2), mar=c(4.2,4,1.2,1), oma=c(1,2,0,0))
partialPlot(RFmig, pred.data=migdep, x.var=dep.long, xlab="origin (longitude of nest)", main="", rug = F, ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1), xlim=c(20,28))
partialPlot(RFmig, pred.data=migdep[!is.na(migdep$tarsus),], x.var=tarsus, main="", rug = F, xlab="tarsus length (cm)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
partialPlot(RFmig, pred.data=migdep, x.var=dep.date, main="", rug = F, xlab="departure week", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
#partialPlot(RFmig, pred.data=migdep[!is.na(migdep$mass),], x.var=mass, main="", rug = F, xlab="body mass (g)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
partialPlot(RFmig, pred.data=migdep[!is.na(migdep$head),], x.var=head, main="", rug = F, xlab="head length (cm)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
#partialPlot(RFmig, pred.data=migdep[!is.na(migdep$tail),], x.var=tail, main="", rug = F, xlab="tail length (cm)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))


#partialPlot(RFmig, pred.data=migdep, x.var=dep.lat, xlab="origin (latitude of nest)", main="", rug = F, ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1), xlim=c(39,44))

#partialPlot(RFmig, pred.data=migdep[!is.na(migdep$winddir),], x.var=winddir, main="", rug = F, xlab="wind direction (degrees)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
#partialPlot(RFmig, pred.data=migdep, x.var=therm.uplift, main="", rug = F, xlab="velocity of thermal uplift (m/s)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))
#partialPlot(RFmig, pred.data=migdep, x.var=fledge.date, main="", rug = F, xlab="fledge week", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F, ylim=c(0,1))

mtext(side=2,"log(probability of successful migration)/2", cex=1.5, outer = T)

dev.off()
hist(migdep$winddir)


### PRODUCE TABLE FOR PAPER ####
table4<-aggregate(winddir~success, data=migdep, FUN=mean)
table4$windspeed<-aggregate(windspeed~success, data=migdep, FUN=mean)[,2]
table4$orog.uplift<-aggregate(orog.uplift~success, data=migdep, FUN=mean)[,2]
table4$therm.uplift<-aggregate(therm.uplift~success, data=migdep, FUN=mean)[,2]
table4$pressure<-aggregate(pressure~success, data=migdep, FUN=mean)[,2]
table4$Cloud<-aggregate(Cloud~success, data=migdep, FUN=mean)[,2]
table4$speed<-aggregate(speed~success, data=migdep, FUN=mean)[,2]
#table4$mass<-aggregate(mass~success, data=migdep, FUN=mean)[,2]
table4$tarsus<-aggregate(tarsus~success, data=migdep, FUN=mean)[,2]
table4$tail<-aggregate(tail~success, data=migdep, FUN=mean)[,2]
table4$head<-aggregate(head~success, data=migdep, FUN=mean)[,2]


table4$dep.date<-aggregate(as.numeric(format(as.POSIXct(mig_summary$start_mig),format="%j"))~mig_completed, data=mig_summary, FUN=mean)[,2]
table4$fledge.date<-aggregate(as.numeric(format(as.POSIXct(mig_summary$Fledge_date),format="%j"))~mig_completed, data=mig_summary, FUN=mean)[,2]
table4$dep.long<-aggregate(Tag_long~mig_completed, data=mig_summary, FUN=mean)[,2]
table4$dep.lat<-aggregate(Tag_lat~mig_completed, data=mig_summary, FUN=mean)[,2]


df.aree<-as.matrix(table4)
df.aree <- as.data.frame(t(df.aree))
colnames(df.aree) <- df.aree[1, ]
df.aree <- df.aree[-1, ]
df.aree$variable <- as.character(row.names(df.aree))
IMP$variable <- as.character(row.names(IMP))
table4<-merge(df.aree,IMP, by='variable')
table4<-table4[order(table4$rel_imp, decreasing=T),]
write.table(table4,"clipboard", sep="\t", row.names=F)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# previous approach of analysis of pre-migration conditions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

premig<-subset(migration, season=="premig")
premig<-subset(premig, Transmitter_ID %in% mig_summary$Transmitter_ID[mig_summary$Age=="juv"])
premig$success<-mig_summary$mig_completed[match(premig$Transmitter_ID, mig_summary$Transmitter_ID)]
#premig$dep.date<-strptime(mig_summary$start_mig[match(premig$Transmitter_ID, mig_summary$Transmitter_ID)],"%Y-%m-%d")$yday				### single departure day
premig$dep.date<-as.numeric(format(as.POSIXct(mig_summary$start_mig[match(premig$Transmitter_ID, mig_summary$Transmitter_ID)]),format="%j"))	### departure week to avoid overfitting by RF
premig$fledge.date<-as.numeric(format(as.POSIXct(mig_summary$Fledge_date[match(premig$Transmitter_ID, mig_summary$Transmitter_ID)]),format="%j"))	### departure week to avoid overfitting by RF
tail(premig)

premig_sum<-aggregate(windspeed~Transmitter_ID, data=premig, FUN=mean)
premig_sum$winddir<-aggregate(winddir~Transmitter_ID, data=premig, FUN=mean)[,2]
premig_sum$orog.uplift<-aggregate(orog.uplift~Transmitter_ID, data=premig, FUN=mean)[,2]
premig_sum$therm.uplift<-aggregate(therm.uplift~Transmitter_ID, data=premig, FUN=mean)[,2]
premig_sum$pressure<-aggregate(pressure~Transmitter_ID, data=premig, FUN=mean)[,2]
premig_sum$Cloud<-aggregate(Cloud~Transmitter_ID, data=premig, FUN=mean)[,2]


### ANALYSIS AT INDIVIDUAL LEVEL - CRAP PERFORMANCE!! ###
mig_summary<-merge(mig_summary,birds[,c(1,2,15,16)], by=c("Transmitter_ID","Name"), all.x=T)
mig_summary$jday<-strptime(mig_summary$start_mig,"%Y-%m-%d")$yday
mig_summary<-merge(mig_summary,premig_sum, by="Transmitter_ID", all.x=T)
write.table(mig_summary,"EGVU_migration_summary.csv", sep=",", row.names=F)
RFmig<-randomForest(as.factor(mig_completed)~jday+Tag_long+winddir+windspeed+orog.uplift+therm.uplift+pressure+Cloud, data=mig_summary[mig_summary$Age=="juv",], mtry=2, ntree=1500, replace=F, importance=T, na.action=na.omit)
varImpPlot(RFmig)
RFmig

head(all_migdata[order(all_migdata$speed, decreasing=T),])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMINE WHETHER WIND WAS STRONGER FOR FAILED MED CROSSING THAN FOR SUCCESSFUL ONE
# NO EFFECT FOUND - NOT INCLUDED IN MS TO REDUCE COMPLEXITY
# 26 Sept 2014
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Wind_U is the ZONAL VELOCITY, i.e. the component of the horizontal wind TOWARDS EAST. 
## Wind_V is the MERIDIONAL VELOCITY, i.e. the component of the horizontal wind TOWARDS NORTH.
## http://mst.nerc.ac.uk/wind_vect_convs.html

failures<-mig_summary[mig_summary$mig_completed==0,]
migfail<-data.frame()

for (a in unique(failures$Name)){
start<-failures$start_mig[failures$Name==a]
x<-migration[migration$Bird_ID==a & migration$Date >= start,]
end<-max(x$Date)
x<-x[x$Date <= end,]
migfail<-rbind(migfail,x)
}
migfail$FAIL<-1


### ADD THE ONE BIRD THAT SURVIVED THE MEDITERRANEAN
start<-mig_summary$start_mig[mig_summary$Name=="Paschalis"]
x<-migration[migration$Bird_ID=="Paschalis" & migration$Date >= start,]
end<-as.POSIXct("2013-09-15")
x<-x[x$Date <= end,]
x$FAIL<-0
migfail<-rbind(migfail,x)

library(lme4)
summary(glmer(FAIL~Wind_V+Wind_U+(1|Bird_ID), data=migfail, family=binomial))
aggregate(Wind_U~FAIL, data=migfail, FUN=mean)
aggregate(Wind_V~FAIL, data=migfail, FUN=mean)
t.test(Wind_V~FAIL, data=migfail)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT WINTERING DATA SET FOR HOME RANGE AND HABITAT ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## need to constrain winter period to [end fall mig - start spring mig]
birds$Stop_date[is.na(birds$Stop_date)]<-Sys.time()
winbirds<-unique(birds$Name[difftime(as.Date(birds$Stop_date),as.Date(birds$Tag_date))>100])
winbirds<-winbirds[-2]			### remove the adult bird
winHRinput<-data.frame()

for (a in winbirds){
input<-migration[migration$Bird_ID == a,]
windata<-input[as.Date(input$Date)>=mig_summary$end_mig[mig_summary$Name==a],]
end<-ifelse(is.na(mig_summary$spring_mig[mig_summary$Name==a]),Sys.time(),mig_summary$spring_mig[mig_summary$Name==a])
windata<-windata[as.Date(windata$Date)<end,]
winHRinput<-rbind(winHRinput, windata)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE AND RE-LOAD WORKSPACE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save.image("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration.RData")
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration.RData")

load("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration.RData")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE SIMPLE 95% KERNEL DENSITIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
winHRinput<-winHRinput[order(winHRinput$Bird_ID,winHRinput$DateTime),]
WHR<-move(x=winHRinput$long, y=winHRinput$lat, time=winHRinput$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"), animal=winHRinput$Bird_ID)
bbox(WHR)
WHR_df <- as(WHR, "data.frame")
WHR<-spTransform(x=WHR, CRSobj="+proj=aeqd", center=TRUE)
winHR<-move2ade(WHR)

### AREA OF HOME RANGE FOR EACH ANIMAL ####
ud<-kernelUD(xy=winHR,h = "href", same4all = F, kern = "bivnorm",extent = 1.0)		### LSCV does not work!
image(ud) ## Note that the contours correspond to values of probability density
HR_summary<-kernel.area(ud, percent = c(50,95),unin = "m", unout = "km2")


### AREA OF WINTER MOVEMENT RANGE FOR EACH ANIMAL ####
HR_summary<-rbind(HR_summary,mcp.area(xy=winHR,percent = c(100),unin = "m", unout = "km2", plotit = TRUE))		### MCP for movement range
HR<-melt(HR_summary)
HR$type<-rep(c('kernel50','kernel95','MCP100'), 9)
names(HR)<-c('Bird_ID','area','type')
aggregate(area~type, HR, FUN=mean)
aggregate(area~type, HR, FUN=sd)
aggregate(area~type, HR, FUN=min)
aggregate(area~type, HR, FUN=max)


## Plotting of the 95 percent home range
ver <- getverticeshr(ud, 95)
plot(ver, col=ver@data$id)
par(new=T)
ver50 <- getverticeshr(ud, 50)
plot(ver50, col=ver@data$id)
axis(1,at=c(100000,300000,500000), labels=c('0', '200', '400 km'), cex=2, cex.lab=1.8, cex.axis=1.8)


### AREA OF HOME RANGE FOR ALL ANIMALS TOGETHER ####
inHRall<-winHR
inHRall@data$id<-"X108050"
udall<-kernelUD(xy=inHRall,h = "href", same4all = F, kern = "bivnorm",extent = 1.0)
image(udall) ## Note that the contours correspond to values of probability density
kernel.area(udall, percent = c(50,95),unin = "m", unout = "km2")


## Plotting of the 95 percent home range
verall <- getverticeshr(udall, 95)
verall50 <- getverticeshr(udall, 50)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOT A MAP OF WINTERING HOME RANGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ver<-spTransform(x=ver, CRSobj=CRS("+proj=longlat +datum=WGS84"), center=TRUE)
ver50<-spTransform(x=ver50, CRSobj=CRS("+proj=longlat +datum=WGS84"), center=TRUE)
verall<-spTransform(x=verall, CRSobj=CRS("+proj=longlat +datum=WGS84"), center=TRUE)
verall50<-spTransform(x=verall50, CRSobj=CRS("+proj=longlat +datum=WGS84"), center=TRUE)

bbox(extent(ver))

verall_df <- fortify(verall)
verall50_df <- fortify(verall50)
ver_df <- fortify(ver)
ver50_df <- fortify(ver50)

#jpeg(filename = "EGVU_win_HR_map.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
pdf("EGVU_win_HR_pop.pdf", width=12, height=11)

### ONE MAP FOR ALL ANIMALS###
jpeg(filename = "A:\\MANUSCRIPTS\\in_prep\\EGVU_migration\\Figure3.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
postscript(file = "A:\\MANUSCRIPTS\\submitted\\EGVU_migration\\Figure3.eps", width = 9, height = 8)

tiff(file = "A:\\MANUSCRIPTS\\in_press\\EGVU_migration\\Figure3.tif", width = 9, height = 8, units = "in", pointsize = 12,
     compression = "none",res=600)
m <- get_map(location=bbox(ver), source="google", zoom=4, col="bw")
ggmap(m)+geom_polygon(data=verall_df, aes(x=long, y=lat, group=piece),colour="black", lwd=1, lty=2, fill=NA)+
geom_polygon(data=verall50_df, aes(x=long, y=lat, group=group),colour="black", lwd=1.5, fill=NA)

dev.off()



#jpeg(filename = "EGVU_win_HR_map.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
pdf("EGVU_win_HR_ind.pdf", width=12, height=11)

### EACH ANIMAL IN A SEPARATE HOME RANGE POLYGON ###
ggmap(m)+geom_polygon(data=ver_df, aes(x=long, y=lat, group=group),colour="blue", fill=NA)+
geom_polygon(data=ver50_df, aes(x=long, y=lat, group=group),colour="red")

dev.off()


library(rgdal)
writeOGR(verall, dsn = '.', layer = 'EGVU_winHomeRange', driver = "ESRI Shapefile")
writeOGR(verall50, dsn = '.', layer = 'EGVU_winHomeRange50', driver = "ESRI Shapefile")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UPLOAD ENV DATA FOR WINTERING BIRDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")

import<-read.table("Greece-6881645433928365752.csv", header=T, sep=",")
import$DateTime<-as.POSIXct(import$timestamp, format="%Y-%m-%d %H:%M:%S", tz='UTC')
import$Bird_ID<-import$individual.local.identifier
import$Transmitter_ID<-import$tag.local.identifier
import$Date<-format(import$DateTime, format="%Y-%m-%d")
import$Time<-format(import$DateTime, format="%H:%M:%S")
import$long<-import$location.long
import$lat<-import$location.lat
import<-import[,c(20:31)]
names(import)[1:5]<-c('habitat','rain','elev','NDVI','GPP')
head(import)
dim(import)
str(import)
str(winHRinput)
WIN_HAB<-merge(winHRinput, import, by=c("Transmitter_ID","Bird_ID","Time","long","lat","DateTime"), all.x=T)
WIN_HAB<-WIN_HAB[,c(1:6,17:21)]
head(WIN_HAB)
dim(WIN_HAB)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLASSIFY LAND COVER DATA INTO COHERENT LANDSCAPE TYPES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WIN_HAB$LC<-ifelse(WIN_HAB$habitat %in% c(11,14,20,30),"cropland",
			ifelse(WIN_HAB$habitat %in% c(60,130,110),"savannah",
			ifelse(WIN_HAB$habitat %in% c(120,140),"grassland",				
			ifelse(WIN_HAB$habitat %in% c(150,200),"desert","other"))))
			#ifelse(WIN_HAB$habitat %in% c(180,190,210),"mesic vegetation",
			#ifelse(WIN_HAB$habitat %in% c(190),"urban",
			#ifelse(WIN_HAB$habitat %in% c(210),"water","other")))))))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OVERLAY LOCATIONS WITH CORE HOME RANGES TO DETERMINE HABITAT IN CORE HOME RANGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CORE HOME RANGE AT POPULATION LEVEL STORED IN THIS spdf
str(ver50)
proj4string(ver50)

### CREATE SPATIAL POINTS OBJECT

WINsp<-SpatialPointsDataFrame(coords=WIN_HAB[,c(4:5)], proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"), data=WIN_HAB, match.ID=T)
head(WINsp)

### OVERLAY POINTS OBJECT WITH POLYGON SHAPE

OVERLAY<-over(WINsp, ver50)		### this does the actual overlay and produces a vector of 1 and NA
WIN_HAB$HR_50<-ifelse(is.na(OVERLAY[,2]),0,1)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OVERLAY LOCATIONS FOR EACH INDIVIDUAL WITH CORE HOME RANGE FOR THAT INDIVIDUAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WIN_HAB<-WIN_HAB[order(WIN_HAB$Bird_ID,WIN_HAB$DateTime),]
WIN_HAB$HR_50ind<-NA

for (a in winbirds){
MOBJ<-move(x=WIN_HAB$long[WIN_HAB$Bird_ID == a], y=WIN_HAB$lat[WIN_HAB$Bird_ID == a], time=WIN_HAB$DateTime[WIN_HAB$Bird_ID == a], proj=CRS("+proj=longlat +ellps=WGS84"), animal=WIN_HAB$Bird_ID[WIN_HAB$Bird_ID == a])
WHR<-spTransform(x=MOBJ, CRSobj="+proj=aeqd", center=TRUE)
winHR<-move2ade(MOBJ)
ud<-kernelUD(xy=winHR,h = "href", same4all = F, kern = "bivnorm",extent = 1.0)		### LSCV does not work!
ver50 <- getverticeshr(ud, 50)
ver50<-spTransform(x=ver50, CRSobj=CRS("+proj=longlat +datum=WGS84"), center=TRUE)

### CREATE SPATIAL POINTS OBJECT

WINsp<-SpatialPointsDataFrame(coords=WIN_HAB[WIN_HAB$Bird_ID == a,c(4:5)], proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"), data=WIN_HAB[WIN_HAB$Bird_ID == a,], match.ID=T)

### OVERLAY POINTS OBJECT WITH POLYGON SHAPE

OVERLAY<-over(WINsp, ver50)		### this does the actual overlay and produces a vector of 1 and NA
WIN_HAB$HR_50ind[WIN_HAB$Bird_ID == a]<-ifelse(is.na(OVERLAY[,2]),0,1)

}

head(WIN_HAB)
sum(WIN_HAB$HR_50)
dim(WIN_HAB)
str(WIN_HAB)

table(WIN_HAB$LC, WIN_HAB$HR_50ind)





####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~ 		MODELLING HABITAT SEGREGATION BETWEEN CORE HOME RANGE AND REST OF WINTER LOCATIONS 					     ~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

WIN_HAB$LC<-as.factor(WIN_HAB$LC)
WIN_HAB$rain <-as.numeric(as.character(WIN_HAB$rain))
WIN_HAB$elev<-as.numeric(as.character(WIN_HAB$elev))
WIN_HAB$NDVI<-as.numeric(as.character(WIN_HAB$NDVI))
WIN_HAB$GPP<-as.numeric(as.character(WIN_HAB$GPP))



####################### UNIVARIATE ANALYSIS OF EACH VARIABLE ################################################

par(mfrow=c(2,2))
for (i in c(8:11)){
P<-wilcox.test(WIN_HAB[,i]~WIN_HAB$HR_50ind)
pval<-ifelse(round(as.numeric(P$p.value),3)==0,"<0.001",round(as.numeric(P$p.value),3))
title<-paste(names(WIN_HAB)[i],pval, sep="  p = ")
boxplot(WIN_HAB[,i]~WIN_HAB$HR_50ind, main=title)
}



####################### EXAMINE HABITAT DIFFERENCES BETWEEN CORE USE AREAS AND AVAILABLE LOCATIONS #########################

CH<-randomForest(HR_50ind ~ (rain) + (elev) + (NDVI) + (GPP) + (LC), data=WIN_HAB, ntree=1500, mtry=5, replace=T, na.action=na.omit, importance=TRUE)
CH
varImpPlot(CH, main="Variables distinguishing between core and peripheral wintering range")



### PLOT VARIABLE IMPORTANCE ####

IMP<-as.data.frame(importance(CH)) 			## SORTED BY MSE
names(IMP)[1]<-"MSE"
IMP<-IMP[order(IMP$MSE, decreasing=T),]  			## SORTED BY Gini
IMP$rel_imp<-round((IMP$MSE/IMP$MSE[1])*100,2)

pdf("EGVU_habitat_selection.pdf", width=20, height=12)
par(mar=c(5,7,2,1))
barplot(IMP$rel_imp[5:1], horiz=T, names.arg=row.names(IMP)[5:1], xlim=c(0,100), las=1,xlab="Relative importance (%)", col='darkred',main="Importance of variables for habitat choice")
text(100,1,sprintf("Var explained = %s",round(CH$rsq[1500],3)),cex=1.5, adj=1)
dev.off()



### PLOT ENVIRONMENT RELATIONSHIPS ####
source("A:\\RSPB\\Statistics\\partialPlot.r")

jpeg(filename = "A:\\MANUSCRIPTS\\in_prep\\EGVU_migration\\Figure3hab.jpg", width = 9, height = 8, units = "in", pointsize = 12, res=600,quality = 100)
postscript(file = "A:\\MANUSCRIPTS\\submitted\\EGVU_migration\\Figure4.eps", width = 9, height = 8)

par(mfrow=c(2,2), mar=c(4,4,1,0), oma=c(1,2,0,0))
partialPlot(CH, pred.data=WIN_HAB[!is.na(WIN_HAB$elev),], ylim=c(0,0.7), x.var=elev, main="", rug = F, xlab="altitude (m above sea level)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F)
partialPlot(CH, pred.data=WIN_HAB, x.var=LC, xlab="landscape type", main="", ylim=c(0,0.7), rug = F, ylab="", las=1, cex.lab=1.2, cex.axis=1.2)
partialPlot(CH, pred.data=WIN_HAB[!is.na(WIN_HAB$GPP),], x.var=GPP, main="", ylim=c(0,0.7), rug = F, xlab="gross primary productivity", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F)
partialPlot(CH, pred.data=WIN_HAB[!is.na(WIN_HAB$NDVI),], x.var=NDVI, main="", ylim=c(0,0.7), rug = F, xlab="Normalised difference vegetation index (NDVI)", ylab="", las=1, cex.lab=1.2, cex.axis=1.2, frame=F)
mtext(side=2,"log(probability of core use)/2", cex=1.5, outer = T)

dev.off()

partialPlot(CH, pred.data=WIN_HAB[!is.na(WIN_HAB$rain),], x.var=rain)






########## ~~~~~~~~~~~~~~ Conditional RandomForest ~~~~~~~~~~~~~~~~~~~#############################
##### abandoned 24 July 2014 because insufficient memory to process this model

my_cforest_control <- cforest_control(teststat="max", testtype="Univ", mincriterion=qnorm(0.6), savesplitstats=FALSE, ntree=500, mtry=5, replace=F, fraction=0.65)
RF1000<-cforest(HR_50ind ~ (rain) + (elev) + (NDVI) + (GPP) + (LC), data=WIN_HAB, controls=my_cforest_control)
VAR<-varimp(RF1000, mincriterion=0.0)
IMP<-data.frame(variable=names(VAR), Gini=VAR)



WIN_HAB$pred.HR<-predict(RF1000)
cor<-cor.test(WIN_HAB$pred.HR,WIN_HAB$HR_50ind)




##### SUMMARISE HABITAT FOR TABLE 3 ######

head(WIN_HAB)
dim(WIN_HAB)
out<-as.data.frame(table(WIN_HAB$LC))
names(out)<-c("variable","value")
out$use<-out$value/sum(out$value)

out<-as.data.frame(table(WIN_HAB$LC, WIN_HAB$HR_50ind))
names(out)<-c("variable","HR_50ind","value")
out$use[out$HR_50ind==0]<-out$value[out$HR_50ind==0]/sum(out$value[out$HR_50ind==0])
out$use[out$HR_50ind==1]<-out$value[out$HR_50ind==1]/sum(out$value[out$HR_50ind==1])
out<-out[,c(1,2,4)]
write.table(out,"clipboard", sep="\t", row.names=F)

out2<-aggregate(rain~HR_50ind, data=WIN_HAB, FUN=mean)
out2$sd<-aggregate(rain~HR_50ind, data=WIN_HAB, FUN=sd)[,2]
write.table(out2,"clipboard", sep="\t", row.names=F)

out2<-aggregate(elev~HR_50ind, data=WIN_HAB, FUN=mean)
out2$sd<-aggregate(elev~HR_50ind, data=WIN_HAB, FUN=sd)[,2]
write.table(out2,"clipboard", sep="\t", row.names=F)

out2<-aggregate(NDVI~HR_50ind, data=WIN_HAB, FUN=mean)
out2$sd<-aggregate(NDVI~HR_50ind, data=WIN_HAB, FUN=sd)[,2]
write.table(out2,"clipboard", sep="\t", row.names=F)

out2<-aggregate(GPP~HR_50ind, data=WIN_HAB, FUN=mean)
out2$sd<-aggregate(GPP~HR_50ind, data=WIN_HAB, FUN=sd)[,2]
write.table(out2,"clipboard", sep="\t", row.names=F)














#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA INVENTORY TO CHECK WHETHER ALL LOCATIONS ARE IN DATABASE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set working directory
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\Tracking_data")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Raw_Data\\Tracking_data")

pdf("EGVU_sat_track_data.pdf", height=14, width=9)
par(mfrow=c(6,4))
daterange<-seq(min(migration$Date), max(migration$Date),by='1 month')
for (a in anim){
x<-subset(migration,Bird_ID==a)
plottitle<-paste(a,birds$Transmitter_ID[match(a,birds$Name)], sep=" ")
hist(x$Date, breaks="months", xlim=c(min(migration$Date),max(migration$Date)), axes=F, xlab="", ylab="n locations", main=plottitle)
axis(1, at= daterange, labels=format(daterange, "%b %y"), tck=-0.02, mgp=c(3,1,0), cex.axis=1.2)
abline(v=birds$Tag_date[match(a,birds$Name)], col='red')
}
dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DISCONTINUED CODE OF DYN BBMM - did not yield satisfactory home range maps for flyway or wintering area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING A RASTER FOR AUTUMN MIGRATION PATHWAY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")

#### CREATE A RASTER FOR THE MIGRATION FLYWAY EXTENT #####
flyway<-readOGR("EGVU_migration_extent.kml", "EGVU_migration_extent.kml")
#str(flyway)
#extent(flyway)
r <- raster(extent(flyway), ncol=100, nrow=100, crs=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
setValues(r, values=0)
rp <- rasterize(flyway, r, field="Name", mask=F, background=NA)
flyway_raster<- projectRaster(from=rp, crs=proj4string(MOBJ_all))

plot(MOBJ_all)
plot(flyway_raster, add=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BROWNIAN BRIDGE MOVEMENT MODEL FOR AUTUMN MIGRATION PATHWAY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### SUBSET ALL TRACKING DATA TO ONLY USE MIGRATION DATA ########
all_migdata<-all_migdata[!all_migdata$Transmitter_ID %in% c(118149,130675),]
head(all_migdata)

## need to constrain analyses to birds that actually completed migration!
birds$Stop_date[is.na(birds$Stop_date)]<-Sys.time()
winbirds<-unique(birds$Name[difftime(as.Date(birds$Stop_date),as.Date(birds$Tag_date))>100])

## create output raster that contains 0 values in each cell
flyway_output<- flyway_raster
flyway_output@data@values <- ifelse(is.na(flyway_output@data@values),NA,0)



par(mfrow=c(4,3))

for (a in winbirds){

input<-all_migdata[all_migdata$Bird_ID == a,]

MOBJ<-move(x=input$long, y=input$lat, time=input$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"), animal=input$Bird_ID)
MOBJ<-spTransform(x=MOBJ, CRSobj=proj4string(flyway_raster))
#MOBJ_df <- as(MOBJ, "data.frame")
#BB<-brownian.bridge.dyn(MOBJ, dimSize=500, location.error=25,margin=15, window.size=35)
#plot(BB, xlab="", ylab="", mar=c(2,1,0,0), main=paste(a,birds$Name[match(a,birds$Transmitter_ID)], sep=" "))
BB<-brownian.bridge.dyn(MOBJ, raster=flyway_raster, location.error=25,margin=15, window.size=35)
plot(BB, xlab="", ylab="", mar=c(2,1,0,0), main=paste(a,birds$Name[match(a,birds$Transmitter_ID)], sep=" "))
#con<-raster2contour(BB, levels=c(.5,.95))
#plot(con, add=T)
flyway_output@data@values <- flyway_output@data@values+(BB@data@values)
lines(MOBJ, col='lightgray', lwd=0.5)
contour(BB, levels=c(0.95), add=TRUE, lwd=1)
contour(BB, levels=c(0.50), add=TRUE, lty=2, lwd=1)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARIZE UTILISATION DISTRIBUTION ACROSS ALL INDIVIDUALS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## THIS WILL NOT SHOW A USEFUL MAP, BECAUSE THE NARROW FLYWAY IS NOT WHERE BIRDS SPEND TIME< SO THE UD IS VERY LOW
## HEAVY BIAS TOWARDS STOPOVER AREAS


totalUD<-sum(flyway_output@data@values, na.rm=T)			### this should sum to the number of individuals
totalUD<- flyway_output@data@values/totalUD
flyway_output@data@values <- totalUD
#min(flyway_output@data@values, na.rm=T)
flywayUDcon<-rasterToContour(flyway_output, levels = c(0.5,0.95)) 				### does not work when levels are specified: rasterToContour(flyway_output, levels = c(0.5,0.95))
plot(flywayUDcon)

flywayUDcon<-spTransform(x=flywayUDcon, CRSobj=CRS("+proj=longlat +datum=WGS84"))
flywayUDcon_df <- fortify(flywayUDcon)
head(flywayUDcon_df)

m <- get_map(bbox(extent(flywayUDcon)), source="google", zoom=4)
ggmap(m)+geom_polygon(data=flywayUDcon_df[flywayUDcon_df$piece==1,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=flywayUDcon_df[flywayUDcon_df$piece==2,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=flywayUDcon_df[flywayUDcon_df$piece==3,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=flywayUDcon_df[flywayUDcon_df$piece==4,], aes(x=long, y=lat, group=id), col=2)


#ggmap(m)+geom_point(data=all_migdata, aes(x=long, y=lat, group=Bird_ID), col=2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATING A RASTER FOR WINTER HOMERANGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")

#### CREATE A RASTER FOR THE WINTERING EXTENT #####
subsahara<-readOGR("EGVU_winHR_extent.kml", "EGVU_winHR_extent.kml")
r <- raster(extent(subsahara), ncol=50, nrow=50, crs=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
setValues(r, values=0)
rp <- rasterize(subsahara, r, field="Name", mask=F, background=NA)
subsahara_raster<- projectRaster(from=rp, crs=proj4string(MOBJ_all))

plot(MOBJ_all)
plot(subsahara_raster, add=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANALYSE WINTER HOME RANGES FOR INDIVIDUAL ANIMALS THAT ARRIVED IN AFRICA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## need to constrain winter period to [end fall mig - start spring mig]
birds$Stop_date[is.na(birds$Stop_date)]<-Sys.time()
winbirds<-unique(birds$Name[difftime(as.Date(birds$Stop_date),as.Date(birds$Tag_date))>100])
winbirds<-winbirds[-2]			### remove the adult bird
winHRinput<-data.frame()


## create output raster that contains 0 values in each cell
subsahara_output<- subsahara_raster
subsahara_output@data@values <- ifelse(is.na(subsahara_output@data@values),NA,0)





par(mfrow=c(3,3))

for (a in winbirds){

input<-migration[migration$Bird_ID == a,]
windata<-input[as.Date(input$Date)>=mig_summary$end_mig[mig_summary$Name==a],]
end<-ifelse(is.na(mig_summary$spring_mig[mig_summary$Name==a]),Sys.time(),mig_summary$spring_mig[mig_summary$Name==a])
windata<-windata[as.Date(windata$Date)<end,]
winHRinput<-rbind(winHRinput, windata)

MOBJ<-move(x=windata$long, y=windata$lat, time=windata$DateTime, proj=CRS("+proj=longlat +ellps=WGS84"), animal=windata$Bird_ID)
MOBJ<-spTransform(x=MOBJ, CRSobj=proj4string(subsahara_raster))
BB<-brownian.bridge.dyn(MOBJ, raster=subsahara_raster, location.error=25,margin=15, window.size=35)
cont <- getVolumeUD(BB)
cont <- cont<=.95
mig_summary$win_95HR_BB[mig_summary$Name==a]<- sum(values(cont))
MCP<-hrBootstrap(x=MOBJ, rep=25, unin='m', unout='km2', plot=F)
mig_summary$win_MCP[mig_summary$Name==a]<- max(MCP[,3])/100
plot(BB, xlab="", ylab="", mar=c(2,1,0,0), main=paste(a,birds$Transmitter_ID[match(a,birds$Name)], sep=" "))
subsahara_output@data@values <- subsahara_output@data@values+(BB@data@values)
lines(MOBJ, col='lightgray', lwd=0.5)
contour(BB, levels=c(0.95), add=TRUE, lwd=1)
contour(BB, levels=c(0.50), add=TRUE, lty=2, lwd=1)
}


max(BB@data@values, na.rm=T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE MEAN UD ACROSS ALL DYN BBMM MODELS FROM ALL INDIVIDUALS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#totalUD<-sum(subsahara_output@data@values, na.rm=T)			### this should sum to the number of individuals
#totalUD<- subsahara_output@data@values/totalUD				### this should ensure that the grid sums to 1 again (mean across all individuals), but leads to tiny values
#subsahara_output@data@values <- totalUD
#min(subsahara_output@data@values, na.rm=T)
subsaharaUDcon<-rasterToContour(subsahara_output, levels = c(0.1,0.4)) 	### does not work when large levels are specified: if no cell exceeds threshold then error: no contour lines
plot(subsaharaUDcon)

subsaharaUDcon<-spTransform(x=subsaharaUDcon, CRSobj=CRS("+proj=longlat +datum=WGS84"))
subsaharaUDcon_df <- fortify(subsaharaUDcon)
head(subsaharaUDcon_df)

m <- get_map(bbox(extent(subsaharaUDcon)), source="google", zoom=5)
ggmap(m)+geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==1,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==2,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==3,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==4,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==5,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==6,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==7,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==8,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==9,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==10,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==11,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==12,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==13,], aes(x=long, y=lat, group=id), col=2)+
geom_polygon(data=subsaharaUDcon_df[subsaharaUDcon_df$piece==14,], aes(x=long, y=lat, group=id), col=2)


#ggmap(m)+geom_point(data=all_migdata, aes(x=long, y=lat, group=Bird_ID), col=2)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE CORRIDORS FOR ALL ANIMALS [abandoned, as the corridor points appeared in areas that are not migration corridors]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MOBJ_corr<-move(x=migration$long[migration$Bird_ID==anim[2]], y=migration$lat[migration$Bird_ID==anim[2]], time=migration$DateTime[migration$Bird_ID==anim[2]], proj=CRS("+proj=longlat +ellps=WGS84"), animal=migration$Bird_ID[migration$Bird_ID==anim[2]])

### search for animal-defined corridors
corridor(MOBJ_corr, speedProp=.10, circProp=.90, plot=TRUE)

### explore speed and azimuth variance parameters
corridor(MOBJ_corr, speedProp=.85, circProp=.20, plot=TRUE)


