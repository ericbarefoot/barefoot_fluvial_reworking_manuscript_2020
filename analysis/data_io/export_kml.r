library(ggmap)
library(maptools)
library(rgdal)



# load libraries

library(tidyverse)
library(here)
library(broom)
library(lme4)
library(devtools)  # for the next thing
library(colorspace)
library(purrr)
library(rsample)

load_all('/Users/ericfoot/Dropbox/programs/git_repos/paleohydror')

library(conover.test)

if (interactive()) {
    require(colorout)
}

format_function = function(pval, signif = 0.001) {
  if (pval > signif) {
    prep = paste('p =', signif(pval, 3))
  } else {
    prep = paste('p \\ll', signif)
  }
  return(prep)
}

# load data

# chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
# chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
# barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

field_depths = field %>% filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>% mutate(meas_a = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b), formation = as.factor(formation))

depths_export = field %>%
filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>%
mutate(meas_a = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b))

grain_size_export = field %>%
filter(meas_type == 'sample' & !(structure %in% c('outcrop', 'formation', 'bar_drape')) & meas_a < 2000) %>%
mutate(D50 = meas_a * 1e-6)

coordinates(depths_export) = c('loc_x', 'loc_y')
proj4string(depths_export)<-CRS("+proj=longlat +datum=WGS84")
coordinates(grain_size_export) = c('loc_x', 'loc_y')
proj4string(grain_size_export)<-CRS("+proj=longlat +datum=WGS84")

flight_data = read_csv(here('data', 'raw_data', 'model_data', 'drone_flight_location_data.csv'))

interpreted_photopanels = flight_data %>% filter(interpreted == 'yes')

coordinates(flight_data) = c('lon', 'lat')
proj4string(flight_data)<-CRS("+proj=longlat +datum=WGS84")
coordinates(interpreted_photopanels) = c('lon', 'lat')
proj4string(interpreted_photopanels)<-CRS("+proj=longlat +datum=WGS84")

writeOGR(depths_export, dsn="barefoot_flow_depths.kml", layer= "flow_depths", driver="KML")
writeOGR(grain_size_export, dsn="barefoot_grain_size_samples.kml", layer= "grain_size_samples", driver="KML")
writeOGR(flight_data, dsn="barefoot_all_3d_outcrops.kml", layer= "drone_flights", driver="KML")
writeOGR(interpreted_photopanels, dsn="barefoot_interpreted_photopanels.kml", layer= "interpreted_panels", driver="KML")
# writeOGR(grain_size_export, dsn="barefoot_example_photos.kml", layer= "interpreted_panels", driver="KML")


data(crime)                                 #Load the Crime Data Set
crime<-crime[!is.na(crime$lat),]            #Get rid of the NA's in the crime data set
crime<-crime[1:100,]                        #Get only the first 100 records

coordinates(crime)<-c("lon","lat")          #Build a SpatialPointsData Frame
 proj4string(crime)<-CRS("+proj=longlat +datum=WGS84")

##The two lines below convert the month and day columns to character dat
##(both of these line are originally 'factor' data, which is not compatible)
crime$month<-as.character(crime$month)
crime$day<-as.character(crime$day)

##The one line of code below writes our SpatialPointsDataFrame to a KML File
writeOGR(crime, dsn="crime.kml", layer= "crime", driver="KML")
