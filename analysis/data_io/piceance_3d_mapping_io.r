#!/usr/bin/env Rscript
# script for importing shapefile data derived from 3d models 
# and putting it into a data table format to join with field data.
#	Eric Barefoot
#	Sep 2019

#	load required packages

require(tidyverse)
require(here)
require(devtools)	# for the next thing
require(paleohydror)
require(sf)

# read in shapefile data

# fields it will come with are:
# 	FID 					This is a unique id for the GIS, will be reassigned in here
# 	interpretations		Denotes full or partial preservation or truncation or type of scour

# the rest of the necessary information will be read from an auxiliary file as part of this i/o script.
# this information includes the date, a referance waypoint number, as listed below
# 	structure		is it a bar? an intra-belt scour? an avulsion scour?
# 	wp_num			waypoint using the convention in other piceance data
# 	subset			in this particular case, we will use to group bars together
# 	date				what day??
	
model_data_file = here('data','raw_data','model_data')

dirs_to_iterate = list.dirs(model_data_file, recursive = FALSE)

files = vector('list', length = length(dirs_to_iterate))
joined_tab = vector('list', length = length(dirs_to_iterate))
attribute_tab = vector('list', length = length(dirs_to_iterate))
site_attributes = vector('list', length = length(dirs_to_iterate))
j = 1

for (i in dirs_to_iterate) {
	filelist = list.files(i)
	if (length(filelist > 0))
	{
		site_match = grep('*site_vals.csv', filelist)
		sa = read_csv(file.path(i, filelist[site_match]), col_types = "Dccc")
		match = grep('*.shp', filelist)
		files[[j]] = st_read(file.path(i, filelist[match]))
		files[[j]] = st_transform(files[[j]], 2152) # transform into UTM coordinates for 12 N
		geom_tab = st_coordinates(files[[j]]) %>% as.data.frame()
		att_tab = as.data.frame(files[[j]]) %>% select(-geometry)
		att_tab$FID = att_tab$FID + 1
		attribute_tab[[j]] = att_tab %>% rename(id = 'FID') %>% select(-LAYER) %>% 
		mutate(id = id + (j * 1000)) %>% 
		mutate(wp_num = sa$wp_num, date = sa$date, formation = sa$formation, field_site = sa$field_site)
		joined_tab[[j]] = inner_join(att_tab, geom_tab, by = c('FID' = 'L1'))
		joined_tab[[j]] = joined_tab[[j]] %>% rename(id = 'FID', interpretations = 'LAYER')
		joined_tab[[j]]$id = joined_tab[[j]]$id + (j * 1000)
		j = j + 1
	}
	else 
	{
		print(paste('there is no data in', i, '.'))
	}
}

bound_data = bind_rows(joined_tab) %>% as_tibble() %>%  
mutate(id = as.factor(id), interpretations = as.character(interpretations))

# bound_data = bound_data 

no_geom = bind_rows(attribute_tab) %>% as_tibble() %>% 
mutate(id = as.factor(id))

# There is some concern about how to handle it if there's different id numbers for each of the shapes across all the files, for now, that's being handled by having all of them with j*1000 added to the id number. This will not be sustainable long-term. 

#  %>% rowid_to_column(var = "Uid")
# 
# do.call(dplyr::recode, c(list(bound_data), lookup))
# 
# 
# bound_data = bound_data %>% recode(id, no_geom$Uid)

# calculate the vertical relief on a structure.

barThicknessMeas = bound_data %>% group_by(id, interpretations) %>% 
summarize(meas_type = 'thickness', meas_a = diff(range(Z))) 

barFaceLengthMeas = bound_data %>% group_by(id) %>% 
summarize(meas_type = 'face_length', meas_a = barfaceDistFun(X,Y,Z)) 

# TODO: find a way to calculate the centroid of every id, and put it into loc_x loc_y and loc_z but in WGS84

measurements = bind_rows(barThicknessMeas, barFaceLengthMeas)

model_data = left_join(measurements, no_geom, by = 'id') %>% as_tibble() %>% select(-id)

# wrangle shapefile data so that it matches field_data format with whatever is required

model_data %>% mutate(units_a = 'm')

save(model_data, file = here('data','derived_data', 'piceance_3d_model_data.rda'))
