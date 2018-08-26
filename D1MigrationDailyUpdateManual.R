##########################################################################
# DAILY REPORTING OF MIGRATION DATA
# original script written by Steffen Oppel
# AUGUST 2018
# data collected by Survey123 app are stored online 
# workflow to create HTML output
# modified for manual download as Python not working
##########################################################################

### CreationDate and EditDate are the exact times in UTC when data were submitted to ArcGISOnline (not when they were entered)
### Date is in UTC converted from the phone's time zone at 00:00 on the day that was manually entered
### https://doc.arcgis.com/en/arcgis-online/create-maps/work-with-fields.htm
### because R cannot store different time zones in a single column, I opted for a make-shift solution
### time zone set to Bangkok (i.e. further east than any project site) to ensure that entered date is also date recorded

# Load necessary library and specify 

library(tidyverse)
library(data.table)
library(lubridate)
library(data.table)
library(markdown)
library(rmarkdown)
library(knitr)
library(plotly)


TODAY<-Sys.time()		###-day(1) if you want to do it on the following day
species<-fread("C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\Species_list.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOWNLOAD DATA MANUALLY FROM ARCGISONLINE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download zip file containing 2 csv files:
# Go to the following web address ( you need your ArcGISOnline password and user ID = the same as for Survey123)
##   http://rspb.maps.arcgis.com/home/item.html?id=9f96dae4476b4b5aa1f9e5c65020613e
# In the top right corer click on the blue 'Download' button

# Save the zip file into this folder:
# C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\D1Migration



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MANIPULATE MIGRATION MONITORING DATA FOR SIMPLE SUMMARY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


############ IMPORT DOWNLOADED CSV FILE ############################
setwd("C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\D1Migration")
outDir<-"C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\D1Migration" # Define the folder where the zip file should be unzipped to 
zipfile<-list.files(pattern=".zip")
unzip(zipfile,exdir=outDir)

files<-list.files(pattern=".csv")

import<-list()
for (f in 1:length(files)){
  import[[f]]<-fread(files[f])
}


### FORMAT MIGRATION SURVEYS

names(import[[1]])

MigSurveys<-import[[1]][,c(2,5,3,4,13)] %>%			### select only columns needed
	mutate(date=mdy_hms(Date)) 	%>%				### extract the date from the date column
	mutate(date=force_tz(date,"UTC")) %>%			### ArcGISOnline stores everything as UTC
	mutate(date=with_tz(date,"Asia/Bangkok")) %>%   ### because we cannot adjust to local time zone, we pick one that is far enough east to not affect dates
  dplyr::select(1,2,6,5) %>%
  filter(date>ymd("2018-08-25"))  ### exclude fake test surveys before 26 Aug 2018


### FORMAT MIGRATION COUNTS

names(import[[2]])

MigCount<-import[[2]][,c(10,3,4)] %>%						### select only columns needed
  filter(ParentGlobalID %in% MigSurveys$GlobalID) %>% ### make sure that fake test surveys are excluded from count data
  setNames(c('GlobalID','Species','N'))

OUT<-merge(MigSurveys,MigCount, by="GlobalID", all=T)
	



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARISE DATA PER DAY AND CREATE TABLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUMMARY<- OUT %>% group_by(date,Species) %>%
  summarise(DailyTotal=sum(N)) %>%
  filter(format(date,'%Y-%m-%d')==format(TODAY,'%Y-%m-%d')) %>%
  mutate(Species=species$EnglishName[match(Species,species$SpeciesCode)])


TOTALSUMMARY<- OUT %>% group_by(Species) %>%
  summarise(AnnualTotal=sum(N)) %>%
  mutate(Species=species$EnglishName[match(Species,species$SpeciesCode)])




##################################################################
### PRODUCE OUTPUT HTML WITH KEY TABLES ###
##################################################################

### create HTML report for overall summary report
#Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

rmarkdown::render('C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\D1Migration\\DailyMigrationSummary.Rmd',
                  output_file = "EGVU_DailyMigrationSummary.html",
                  output_dir = 'C:\\STEFFEN\\RSPB\\Bulgaria\\Raw_Data\\NewLIFE_DATA\\Survey123Downloads\\D1Migration')


