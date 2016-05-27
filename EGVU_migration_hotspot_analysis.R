##########################################################################
# IDENTIFICATION OF MIGRATORY HOTSPOTS FOR SATELLITE TRACKED EGYPTIAN VULTURES
# original script written by Steffen Oppel in September 2014
# modified in February 2016
# update 6 Feb 2016: revised data cleanup
# removed data prep and cleanup - see extra script 'EGVU_migration_data_prep2016.r'
# update 7 Feb 2016: manually identified spring migration events
# update 7 Feb 2016: collated all migration data and ran tripGrid to identify migration hotspots
# update 8 Feb 2016: added home range estimation, revised autumn migration definition for juveniles
# update 9 Feb 2016: added manual autumn migration definition for immature birds, removed home range estimation, added time discretisation and philopatry analysis
# update 10 February: split from script 'EGVU_migration_season_definition.R'
##########################################################################

# Load necessary library
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
library(plotKML)
library(ggplot2)
library(adehabitatHR)
library(reshape)
library(rworldmap)
library(trip)
data(countriesLow)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PREVIOUSLY SAVED WORKSPACE (from 'EGVU_migration_season_definition.R')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set working directory

setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria")

load("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_Feb2016.RData")
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_Feb2016.RData")
load("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria\\EGVU_migration_Feb2016.RData")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT MIGRATION DATA TO TRAJECTORY OBJECT FOR INTERPOLATION - HELPFUL FOR BOTTLENECK IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PREPARE DATA WITH UNIQUE ID
all_migdata$MigID<-paste(all_migdata$Bird_ID,all_migdata$season,all_migdata$year, sep="_")		## this creates a unique identifier for each migration
all_migdata$DateTime<-adjust.duplicateTimes(all_migdata$DateTime, all_migdata$MigID)			## make sure no duplicate time stamps are present
all_migdata<-all_migdata[order(all_migdata$MigID,all_migdata$DateTime),]
head(all_migdata)
str(all_migdata)

### Convert to LTRAJ TO INTERPOLATE DATA
traj<-as.ltraj(xy=all_migdata[,8:9],date=all_migdata[,17],id=all_migdata[,27],infolocs=all_migdata[,c(1,26,25)],typeII=TRUE)

## Rediscretization every 300 seconds
tr <- redisltraj(traj, 300, type="time")

## Convert output into a data frame
str(tr[[1]])
migdat5min<-data.frame()

for (l in 1:length(unique(all_migdata$MigID))){
out<-tr[[l]]
out$MigID<-as.character(attributes(tr[[l]])[4])				#### extracts the MigID from the attribute 'id'
migdat5min<-rbind(migdat5min,out)				#### combines all data
}

head(migdat5min)
dim(migdat5min)
str(migdat5min)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATORY HOTSPOTS AND BOTTLENECKS BY SUMMARISING LOCATIONS PER 100x100km GRID CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A TRIPS OBJECT ###
#EVSP <- SpatialPoints(data.frame(all_migdata$long, all_migdata$lat), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
EVSP <- SpatialPoints(data.frame(migdat5min$x, migdat5min$y), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the 5 min interpolated data
EVSP <- SpatialPointsDataFrame(EVSP, data = migdat5min)		### switch to data = all_migdata if needed
EV_trips<-trip(EVSP, TORnames=c("date","MigID"))			### switch to "DateTime" when using the raw locations


### PLOT THE MIGRATION ROUTES ###
plot(EV_trips,pch=16, cex=0.1, axes = TRUE)
lines(EV_trips)
plot(countriesLow, add=T)


### CREATE A GRID AND COUNT NUMBER OF LOCATIONS IN EACH GRID CELL ###
grd<-makeGridTopology(EV_trips, cellsize = c(100,100), adjust2longlat = TRUE)
trg <- tripGrid(EV_trips, grid=grd,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
spplot(trg)			## plots the trips with a legend
#image(trg)


### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
poly_count<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
poly_count$time<-poly_count$time/(3600*24)								### this converts the seconds into bird days
poly_count<-poly_count[(poly_count$time>5),]					### remove the cells with 0 time

summary(poly_count)
#bbox(EV_trips)
plot(spdf)
summary(spdf)

## EXPORT TO ARCGIS
ras1<- raster(spdf)
spoldf <- rasterToPolygons(ras1, n=4)
writeOGR(spoldf, ".", "EGVU_mig_hotspots_poly", driver="ESRI Shapefile")


### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT MIGRATION SITES ###

MAP <- get_map(bbox(EV_trips), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)

ggmap(MAP)+geom_tile(data=poly_count, aes(x=long,y=lat, fill = time)) +
	scale_fill_gradient(name = 'N bird days', low="yellow", high="red", na.value = 'white', guide = "colourbar", limits=c(5, 31))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATION BOTTLENECKS BY CALCULATING THE PROPORTION OF MIGRATIONS GOING THROUGH EACH CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for each cell we need to list the individual ID's that pass through the cell
# divide that number by the total number of migrations

### OVERALL MAP - [Clementine can do seperate maps for spring and autumn migration]

ras1<- raster(spdf)
spoldf <- rasterToPolygons(ras1, n=4)
proj4string(spoldf)<-proj4string(EVSP)
migrations<-unique(all_migdata$MigID)


bottleneck_plot<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z,n_mig=0)
for (m in migrations){
x<-subset(EVSP, MigID==m)				## improvement: re-discretise time steps to 5 min intervals to not 'miss' any cell.
x<-x[,c(11,3,1,2)]
x$count<-1
proj4string(x)<-proj4string(spoldf)
test<-spoldf %over% x 
test$count[is.na(test$count)]<-0
bottleneck_plot$n_mig<-bottleneck_plot$n_mig+test$count
}

bottleneck_plot$n_mig<-(bottleneck_plot$n_mig/length(migrations))*100
summary(bottleneck_plot)



## EXPORT TO ARCGIS
bottrast<-spdf
bottrast@data$z<-bottleneck_plot$n_mig
ras2<- raster(bottrast)
botpoly <- rasterToPolygons(ras2, n=4)
writeOGR(botpoly, ".", "EGVU_mig_bottlenecks_poly", driver="ESRI Shapefile")


## REDUCE DATA FRAME FOR PLOTTING
bottleneck_plot<-bottleneck_plot[bottleneck_plot$n_mig>25,]



ggmap(MAP)+geom_tile(data=bottleneck_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(20, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE MAPS INTO PDF FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("EGVU_MIGRATION HOTSPOTS.pdf", width=12, height=11)
ggmap(MAP)+geom_tile(data=poly_count, aes(x=long,y=lat, fill = time)) +
	scale_fill_gradient(name = 'N bird days', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(0, 31))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))

ggmap(MAP)+geom_tile(data=bottleneck_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(20, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE WORKSPACE (MIGRATORY BOTTLENECK ANALYSIS AND HOME RANGE ESTIMATION MOVED TO ANOTHER SCRIPT)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set working directory

setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration")
setwd("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria")

save.image("C:\\STEFFEN\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_hotspots.RData")
save.image("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_hotspots.RData")
save.image("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria\\EGVU_migration_hotspots.RData")








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE THE MINIMUM DISTANCE TO THE NEST IN EVERY YEAR AFTER THE TAG YEAR FOR EVERY JUVENILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nest_dist<-data.frame()

head(birds)
juveniles<-birds[(birds$Age=="juv"),c(2,4,5,6,7)]
juveniles<-juveniles[juveniles$Name %in% season_summary$Name,]

for(j in juveniles$Name){
x<-subset(alldata, Bird_ID==j)
x$year<-as.integer(format(x$Date, format="%Y"))
years<-unique(x$year)
years<-years[-1]
for(y in years){

xy<-subset(x, year==y)
if(max(xy$DateTime)>as.POSIXct(sprintf("%s-05-30",y),format="%Y-%m-%d")){
mindist<-xy[xy$home_dist==min(xy$home_dist,na.rm = TRUE),c(1:4,6)]
mindist$year<-y
nest_dist<-rbind(nest_dist,mindist)
}}}		## Close the if statement, year loop, juveniles loop


write.table(nest_dist, "EGVU_philopatry_distance.csv", sep=";", row.names=F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE A PLOT WITH DISTANCE AND THE YEAR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
juveniles<-juveniles[juveniles$Name %in% nest_dist$Bird_ID,]
dim(juveniles)
pdf("EGVU_natal_distance.pdf",width=10, height=7)
par (mfrow=c(2,3), mar=c(3,4,1,0), oma=c(2,2,0,0))

for(j in juveniles$Name){
x<-subset(nest_dist, Bird_ID==j)
x$Age<-x$year-birds$Tag_year[match(j,birds$Name)]
plot(log(home_dist)~Age, type='l', frame=F, data=x, xlim=c(0,4),ylim=log(c(1,5000)), main=j, ylab="",xlab="",axes=F)
axis(2,at=log(c(5,50,500,5000)),labels=c(5,50,500,5000),las=1)
axis(1,at=seq(0,4,1),labels=T,las=1)
}
mtext("Age of bird (year)",1,cex=1.5,outer=T)
mtext("Distance from natal nest (km)",2,cex=1.5,outer=T)
dev.off()




