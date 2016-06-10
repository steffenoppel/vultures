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
# update 28 Apr 2016: create hotspots and bottlenecks maps for autumun and spring separately
# update 27 May 2016 : merge script for spring and autumn bottlenecks to create map with >50% of migrations in BOTH seasons and >70% of migrations either spring or autumn
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
load("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Bulgaria\\Analysis\\Migration\\EGVU_migration_May2016.RData")
load("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria\\EGVU_migration_May2016.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEEP ONLY DATA FOR SPRING MIGRATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spring_migdata<-all_migdata[!(all_migdata$season=="autumn"),]
dim(spring_migdata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT MIGRATION DATA TO TRAJECTORY OBJECT FOR INTERPOLATION - HELPFUL FOR BOTTLENECK IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PREPARE DATA WITH UNIQUE ID
spring_migdata$MigID<-paste(spring_migdata$Bird_ID,spring_migdata$season,spring_migdata$year, sep="_")		## this creates a unique identifier for each migration
spring_migdata$DateTime<-adjust.duplicateTimes(spring_migdata$DateTime, spring_migdata$MigID)			## make sure no duplicate time stamps are present
spring_migdata<-spring_migdata[order(spring_migdata$MigID,spring_migdata$DateTime),]
head(spring_migdata)
str(spring_migdata)

### Convert to LTRAJ TO INTERPOLATE DATA
traj<-as.ltraj(xy=spring_migdata[,8:9],date=spring_migdata[,17],id=spring_migdata[,27],infolocs=spring_migdata[,c(1,26,25)],typeII=TRUE)

## Rediscretization every 300 seconds
tr <- redisltraj(traj, 300, type="time")

## Convert output into a data frame
str(tr[[1]])
spring_migdat5min<-data.frame()

for (l in 1:length(unique(spring_migdata$MigID))){
out<-tr[[l]]
out$MigID<-as.character(attributes(tr[[l]])[4])				#### extracts the MigID from the attribute 'id'
spring_migdat5min<-rbind(spring_migdat5min,out)				#### combines all data
}

head(spring_migdat5min)
dim(spring_migdat5min)
str(spring_migdat5min)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATORY BOTTLENECKS BY SUMMARISING LOCATIONS PER 100x100km GRID CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A TRIPS OBJECT ###
#EVSP <- SpatialPoints(data.frame(spring_migdata$long, spring_migdata$lat), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
EVSP <- SpatialPoints(data.frame(spring_migdat5min$x, spring_migdat5min$y), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the 5 min interpolated data
EVSP <- SpatialPointsDataFrame(EVSP, data = spring_migdat5min)		### switch to data = spring_migdata if needed
EV_spring_trips<-trip(EVSP, TORnames=c("date","MigID"))			### switch to "DateTime" when using the raw locations


### PLOT THE MIGRATION ROUTES ###
plot(EV_spring_trips,pch=16, cex=0.1, axes = TRUE)
lines(EV_spring_trips)
plot(countriesLow, add=T)


### CREATE A GRID AND COUNT NUMBER OF LOCATIONS IN EACH GRID CELL ###
grd<-makeGridTopology(EV_spring_trips, cellsize = c(100,100), adjust2longlat = TRUE)			### CREATE THIS USING all_migdata!
trg <- tripGrid(EV_spring_trips, grid=grd,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
spplot(trg)			## plots the trips with a legend
#image(trg)


### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
spring_poly_count<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
spring_poly_count$time<-spring_poly_count$time/(3600*24)								### this converts the seconds into bird days
spring_poly_count<-spring_poly_count[(spring_poly_count$time>5),]					### remove the cells with 0 time

summary(spring_poly_count)
#bbox(EV_spring_trips)


### EXPORT DATA ###
## Write data into csv table ##

#write.table(spring_poly_count,"EGVU_mig_hotspots_spring.csv", sep=";", row.names=F)			## Export table to use it in GIS - gives only the center of the grid


## Write data into shapefile ##
## First create SpatialPointsDataFrame and define coordinates and coordinates reference system (lat/long WGS84) 

#shape<-SpatialPointsDataFrame(data.frame(spring_poly_count$long, spring_poly_count$lat), spring_poly_count, coords.nrs = numeric(0), proj4string = CRS("+proj=longlat + datum=wgs84"), match.ID = TRUE)
#writeOGR(shape, ".", "spring_hotspots", driver="ESRI Shapefile")

ras1<-raster(spdf)		## converts the SpatialPixelDataFrame into a raster
spoldf<-rasterToPolygons(ras1, n=4, dissolve=T) 		## converts the raster into quadratic polygons

## Export to ArcGis ##

#writeOGR(spoldf, ".", "EGVU_mig_hotspots_spring", driver="ESRI Shapefile")

summary(EV_spring_trips)
summary(trg)

### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT MIGRATION SITES ###

MAP <- get_map(bbox(EV_spring_trips), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)

ggmap(MAP)+geom_tile(data=spring_poly_count, aes(x=long,y=lat, fill = time)) +
	scale_fill_gradient(name = 'N bird days', low="yellow", high="red", na.value = 'white', guide = "colourbar", limits=c(5, 31))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATION BOTTLENECKS IN SPRING BY CALCULATING THE PROPORTION OF MIGRATIONS GOING THROUGH EACH CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for each cell we need to list the individual ID's that pass through the cell
# divide that number by the total number of migrations

### OVERALL MAP -

ras<- raster(spdf)		# converts the SpatialPixelDataFrame into a raster
spoldf <- rasterToPolygons(ras1, n=4) # converts the raster into quadratic polygons
proj4string(spoldf)<-proj4string(EVSP)
spring_migrations<-unique(spring_migdata$MigID)

spring_bottleneck_plot<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z,n_mig=0)
for (m in spring_migrations){
x<-subset(EVSP, MigID==m)				## improvement: re-discretise time steps to 5 min intervals to not 'miss' any cell.
x<-x[,c(11,3,1,2)]
x$count<-1
proj4string(x)<-proj4string(spoldf)
test<-spoldf %over% x 
test$count[is.na(test$count)]<-0
spring_bottleneck_plot$n_mig<-spring_bottleneck_plot$n_mig+test$count
}


### EXPORT DATA ###
#write.table(spring_bottleneck_plot,"EGVU_mig_spring_bottlenecks.csv", sep=",", row.names=F)			## Export table to use it in GIS

## Write data into shapefile ##
## This write only the central point of the 100kmx100km grid

shape_spring<-SpatialPointsDataFrame(data.frame(spring_bottleneck_plot$long, spring_bottleneck_plot$lat), spring_bottleneck_plot, coords.nrs = numeric(0), proj4string = CRS("+proj=longlat + datum=wgs84"), match.ID = TRUE)
#writeOGR(shape_spring, ".", "spring_bottlenecks", driver="ESRI Shapefile")

bottrast<-spdf
bottrast@data$z<-spring_bottleneck_plot$n_mig

## Export to ARCGIS as polygons ##
#writeOGR(spoldf, ".", "EGVU_mig_spring_bottlenecks_poly", driver="ESRI Shapefile")

ras2<-raster(bottrast)
botpoly<-rasterToPolygons(ras2, n=4)
#writeOGR(botpoly, ".", "EGVU_mig_spring_bottlenecks_poly", driver="ESRI Shapefile")

### PLOT AND GET THE MAP ###
spring_bottleneck_plot$n_mig<-(spring_bottleneck_plot$n_mig/length(spring_migrations))*100
summary(spring_bottleneck_plot)

spring_bottleneck_plot<-spring_bottleneck_plot[spring_bottleneck_plot$n_mig>25,]
##spring_bottleneck_plot<-spring_bottleneck_plot[spring_bottleneck_plot$n_mig>20,]		## improvement : get a more sensitive map


ggmap(MAP)+geom_tile(data=spring_bottleneck_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all spring migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(20, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEEP ONLY DATA FOR AUTUMN MIGRATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
autumn_migdata<-all_migdata[!(all_migdata$season=="spring"),]
dim(autumn_migdata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVERT MIGRATION DATA TO TRAJECTORY OBJECT FOR INTERPOLATION - HELPFUL FOR BOTTLENECK IDENTIFICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PREPARE DATA WITH UNIQUE ID
autumn_migdata$MigID<-paste(autumn_migdata$Bird_ID,autumn_migdata$season,autumn_migdata$year, sep="_")		## this creates a unique identifier for each migration
autumn_migdata$DateTime<-adjust.duplicateTimes(autumn_migdata$DateTime, autumn_migdata$MigID)			## make sure no duplicate time stamps are present
autumn_migdata<-autumn_migdata[order(autumn_migdata$MigID,autumn_migdata$DateTime),]
head(autumn_migdata)
str(autumn_migdata)

### Convert to LTRAJ TO INTERPOLATE DATA
traj<-as.ltraj(xy=autumn_migdata[,8:9],date=autumn_migdata[,17],id=autumn_migdata[,27],infolocs=autumn_migdata[,c(1,26,25)],typeII=TRUE)

## Rediscretization every 300 seconds
tr <- redisltraj(traj, 300, type="time")

## Convert output into a data frame
str(tr[[1]])
autumn_migdat5min<-data.frame()

for (l in 1:length(unique(autumn_migdata$MigID))){
out<-tr[[l]]
out$MigID<-as.character(attributes(tr[[l]])[4])				#### extracts the MigID from the attribute 'id'
autumn_migdat5min<-rbind(autumn_migdat5min,out)				#### combines all data
}

head(autumn_migdat5min)
dim(autumn_migdat5min)
str(autumn_migdat5min)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATORY BOTTLENECKS BY SUMMARISING LOCATIONS PER 100x100km GRID CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### CREATE A TRIPS OBJECT ###
#EVSP <- SpatialPoints(data.frame(autumn_migdata$long, autumn_migdata$lat), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the raw location data
EVSP <- SpatialPoints(data.frame(autumn_migdat5min$x, autumn_migdat5min$y), proj4string=CRS("+proj=longlat + datum=wgs84"))		### these are the 5 min interpolated data
EVSP <- SpatialPointsDataFrame(EVSP, data = autumn_migdat5min)		### switch to data = all_migdata if needed
EV_autumn_trips<-trip(EVSP, TORnames=c("date","MigID"))			### switch to "DateTime" when using the raw locations


### PLOT THE MIGRATION ROUTES ###
plot(EV_autumn_trips,pch=16, cex=0.1, axes = TRUE)
lines(EV_autumn_trips)
plot(countriesLow, add=T)

### CREATE A GRID AND COUNT NUMBER OF LOCATIONS IN EACH GRID CELL ###
### INSTEAD OF CREATING A NEW GRID, USE THE SAME AS FOR SPRING MIGRATION
#grd<-makeGridTopology(EV_autumn_trips, cellsize = c(100,100), adjust2longlat = TRUE)
trg <- tripGrid(EV_autumn_trips, grid=grd,method="pixellate")						### this will provide the number of bird seconds spent in each grid cell
spplot(trg)			## plots the trips with a legend
#image(trg)


### CONVERT SPATIAL GRID TO SOMETHING WE CAN PLOT
spdf <- SpatialPixelsDataFrame(points=trg, data=trg@data)
autumn_poly_count<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z)
autumn_poly_count$time<-autumn_poly_count$time/(3600*24)								### this converts the seconds into bird days
autumn_poly_count<-autumn_poly_count[(autumn_poly_count$time>5),]					### remove the cells with 0 time

summary(autumn_poly_count)
#bbox(EV_autumn_trips)


### EXPORT DATA ###
## Write data into csv table ##

#write.table(autumn_poly_count,"EGVU_mig_hotspots_autumn.csv", sep=";", row.names=F)			## Export table to use it in GIS - gives only the center of the grid


## Write data into shapefile ##
## First create SpatialPointsDataFrame and define coordinates and coordinates reference system (lat/long WGS84) 

#shape<-SpatialPointsDataFrame(data.frame(autumn_poly_count$long, autumn_poly_count$lat), autumn_poly_count, coords.nrs = numeric(0), proj4string = CRS("+proj=longlat + datum=wgs84"), match.ID = TRUE)
#writeOGR(shape, ".", "autumn_hotspots", driver="ESRI Shapefile")

ras1<-raster(spdf)		## converts the SpatialPixelDataFrame into a raster
spoldf<-rasterToPolygons(ras1, n=4, dissolve=T) 		## converts the raster into quadratic polygons

## Export to ArcGis ##

#writeOGR(spoldf, ".", "EGVU_mig_hotspots_autumn", driver="ESRI Shapefile")

summary(EV_autumn_trips)
summary(trg)

### PRODUCE NICE AND SHINY MAP WITH THE MOST IMPORTANT MIGRATION SITES ###

MAP <- get_map(bbox(EV_autumn_trips), source="google", zoom=4, color = "bw")		### retrieves a map from Google (requires live internet connection)

ggmap(MAP)+geom_tile(data=autumn_poly_count, aes(x=long,y=lat, fill = time)) +
	scale_fill_gradient(name = 'N bird days', low="yellow", high="red", na.value = 'white', guide = "colourbar", limits=c(5, 31))+
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATION BOTTLENECKS IN AUTUMN BY CALCULATING THE PROPORTION OF MIGRATIONS GOING THROUGH EACH CELL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for each cell we need to list the individual ID's that pass through the cell
# divide that number by the total number of migrations

### OVERALL MAP -

ras<- raster(spdf)		# converts the SpatialPixelDataFrame into a raster
spoldf <- rasterToPolygons(ras1, n=4) # converts the raster into quadratic polygons
proj4string(spoldf)<-proj4string(EVSP)
autumn_migrations<-unique(autumn_migdata$MigID)

autumn_bottleneck_plot<-data.frame(lat=spdf@coords[,2],long=spdf@coords[,1],time=spdf@data$z,n_mig=0)
for (m in autumn_migrations){
x<-subset(EVSP, MigID==m)				## improvement: re-discretise time steps to 5 min intervals to not 'miss' any cell.
x<-x[,c(11,3,1,2)]
x$count<-1
proj4string(x)<-proj4string(spoldf)
test<-spoldf %over% x 
test$count[is.na(test$count)]<-0
autumn_bottleneck_plot$n_mig<-autumn_bottleneck_plot$n_mig+test$count
}



### EXPORT DATA ###
#write.table(autumn_bottleneck_plot,"EGVU_mig_autumn_bottlenecks.csv", sep=",", row.names=F)			## Export table to use it in GIS

## Write data into shapefile ##
## This write only the central point of the 100kmx100km grid

#shape_autumn<-SpatialPointsDataFrame(data.frame(autumn_bottleneck_plot$long, autumn_bottleneck_plot$lat), autumn_bottleneck_plot, coords.nrs = numeric(0), proj4string = CRS("+proj=longlat + datum=wgs84"), match.ID = TRUE)
#writeOGR(shape_autumn, ".", "autumn_bottlenecks", driver="ESRI Shapefile")

bottrast<-spdf
bottrast@data$z<-autumn_bottleneck_plot$n_mig

## Export to ARCGIS as polygons ##
#writeOGR(spoldf, ".", "EGVU_mig_autumn_bottlenecks_poly", driver="ESRI Shapefile")

ras2<-raster(bottrast)
botpoly<-rasterToPolygons(ras2, n=4)
#writeOGR(botpoly, ".", "EGVU_mig_autumn_bottlenecks_poly", driver="ESRI Shapefile")

### PLOT AND GET THE MAP ###
autumn_bottleneck_plot$n_mig<-(autumn_bottleneck_plot$n_mig/length(autumn_migrations))*100
summary(autumn_bottleneck_plot)

autumn_bottleneck_plot<-autumn_bottleneck_plot[autumn_bottleneck_plot$n_mig>25,]
##autumn_bottleneck_plot<-autumn_bottleneck_plot[autumn_bottleneck_plot$n_mig>20,]		## improvement : get a more sensitive map


ggmap(MAP)+geom_tile(data=autumn_bottleneck_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all autumn migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(20, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IDENTIFY MIGRATION BOTTLENECKS WITH >50% OF MIGRATIONS IN BOTH AUTUMN AND SPRING AND >70% OF MIGRATIONS IN EITHER SPRING OR AUTUMN  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### DEFINE THE MOST IMPORTANT CELLS FOR EACH SEASON ###
autumn<-autumn_bottleneck_plot[autumn_bottleneck_plot$n_mig>50,]
spring<-spring_bottleneck_plot[spring_bottleneck_plot$n_mig>50,]

#### QUANTITATIVE INDEX FOR ALL CELLS ####
all_botts<-merge(autumn_bottleneck_plot, spring_bottleneck_plot, by=c("lat","long"), all=T)
names(all_botts)[4]<-"n_mig_autumn"	#autumn						### Rename the columns 4 and 6
names(all_botts)[6]<-"n_mig_spring"	#spring
all_botts$index<-(all_botts$n_mig_autumn*13+all_botts$n_mig_spring*9)/(13+9)
all_botts$rank<-rank(-all_botts$index, ties.method = "min")	

## combine the >50% migrations in both seasons
fifty_aut<-autumn_bottleneck_plot[autumn_bottleneck_plot$n_mig>50,]			### Select the cells with >50% of migrations in autumn, then in spring
fifty_spr<-spring_bottleneck_plot[spring_bottleneck_plot$n_mig>50,]
fifty<-merge(fifty_aut, fifty_spr, by=c("lat","long"), all=T)
names(fifty)[4]<-"n_mig_autumn"	#autumn						### Rename the columns 4 and 6
names(fifty)[6]<-"n_mig_spring"	#spring
names(fifty)
main_bottlenecks_plot<-fifty[!(is.na(fifty$n_mig_autumn)),]				### remove the lines that have no autumn mig info (i.e. those with <50% in autumn)
main_bottlenecks_plot<-main_bottlenecks_plot[!(is.na(main_bottlenecks_plot$n_mig_spring)),]				### remove the lines that have no spring mig info (i.e. those with <50% in spring)
main_bottlenecks_plot

ggmap(MAP)+geom_tile(data=main_bottlenecks_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all autumn migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(25, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE MAPS INTO PDF FILE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("EGVU_MIGRATION MAIN BOTTLENECKS.pdf", width=12, height=11)
ggmap(MAP)+geom_tile(data=main_bottlenecks_plot, aes(x=long,y=lat, fill = n_mig)) +
	scale_fill_gradient(name = '% of all spring migrations', low="white", high="red", na.value = 'white', guide = "colourbar", limits=c(25, 75))+ 
	theme(axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())+
	theme(strip.text.y = element_text(size = 20, colour = "black"), strip.text.x = element_text(size = 15, colour = "black"))
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE WORKSPACE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set working directory

setwd("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria")

save.image("C:\\Users\\Clémentine\\Documents\\Stage M2\\EGVU_Bulgaria\\EGVU_migration_main_bottlenecks.RData")
