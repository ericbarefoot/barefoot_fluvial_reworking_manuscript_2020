#!/usr/bin/env Rscript

#	i/o and processing for field data collected in the Piceance basin
# and also importing shapefile data and putting it into a similar data table format.
#	Eric Barefoot
#	Sep 2019

#	load required packages

require(tidyverse)
require(here)
require(devtools)	# for the next thing
require(paleohydror)
require(sf)
# require(rgdal)

#	set directories

#	load in field data

fieldDataColumnTypes = cols(
	object_number = col_character(),
	unique_id = col_character(),
	date = col_date(format = "%Y-%m-%d"),
	wp_num = col_character(),
	subset = col_character(),
	set_id = col_character(),
	loc_x = col_double(),
	loc_y = col_double(),
	loc_z = col_double(),
	projection = col_character(),
	formation = col_character(),
	field_site = col_character(),
	structure = col_character(),
	meas_type = col_character(),
	meas_a = col_character(),
	units_a = col_character(),
	meas_b = col_double(),
	units_b = col_character(),
	textures = col_character(),
	composition = col_character(),
	samp_ind = col_integer(),
	samp_code = col_character(),
	samp_description = col_character(),
	image_file = col_character(),
	image_code = col_character(),
	observations = col_character(),
	interpretations = col_character(),
	notes = col_character()
)

data_file = here('data', 'raw_data', 'field_data', 'piceance_paleohydro_data - data_table.csv')
n_rows = length(count.fields(data_file))
col_names = names(read_csv(data_file, n_max = 0))
# field_dat = suppressMessages()
field_dat = read_csv(data_file, col_names = col_names, col_types = fieldDataColumnTypes, skip = 3)

#	process and recompile in useful format

#	find rows with only NA values - store as variable
empty = field_dat %>% is.na %>% `!` %>% rowSums

field_dat$empty = empty		#	append to dataframe

field_d = field_dat %>% # clean and filter
	filter(empty != 1) %>%	#	filter out the master row and any empty rows
	select(-empty) %>% #	remove dummy NA tracker variable
	mutate(unique_id = 1:length(date))			#	add unique ID code for each row

az_data = field_d %>%
	select(c(unique_id, meas_a)) %>%
	filter(substr(meas_a,1,1) == 'N' | substr(meas_a,1,1) == 'S') %>%
	rowwise() %>%
	mutate(az = azStringConvert(meas_a))

azimuthReplace = as.numeric(replace(field_d$meas_a, az_data$unique_id, az_data$az))

field_dd = field_d %>% mutate(meas_a = azimuthReplace)

inch_data = field_dd %>%
	select(c(unique_id, meas_a, units_a)) %>%
	filter(grepl('in', units_a)) %>%
	rowwise() %>%
	mutate(meas_a = meas_a * 0.0254) %>%
	mutate(units_a = 'm')

inchReplace = replace(field_dd$meas_a, inch_data$unique_id, inch_data$meas_a)
field_dd = field_dd %>% mutate(meas_a = inchReplace)
inchUnitsReplace = replace(field_dd$units_a, inch_data$unique_id, inch_data$units_a)
field_dd = field_dd %>% mutate(units_a = inchUnitsReplace)

cm_data = field_dd %>%
	select(c(unique_id, meas_a, units_a)) %>%
	filter(grepl('cm', units_a)) %>%
	rowwise() %>%
	mutate(meas_a = meas_a / 100) %>%
	mutate(units_a = 'm')

cmReplace = replace(field_dd$meas_a, cm_data$unique_id, cm_data$meas_a)
field_dd = field_dd %>% mutate(meas_a = cmReplace)
cmUnitsReplace = replace(field_dd$units_a, cm_data$unique_id, cm_data$units_a)
field_dd = field_dd %>% mutate(units_a = cmUnitsReplace)

cm_data_b = field_dd %>%
	select(c(unique_id, meas_b, units_b)) %>%
	filter(grepl('cm', units_b)) %>%
	rowwise() %>%
	mutate(meas_b = meas_b / 100) %>%
	mutate(units_b = 'm')

cmReplaceB = replace(field_dd$meas_b, cm_data_b$unique_id, cm_data_b$meas_b)
field_dd = field_dd %>% mutate(meas_b = cmReplaceB)
cmUnitsReplaceB = replace(field_dd$units_b, cm_data_b$unique_id, cm_data_b$units_b)
field_dd = field_dd %>% mutate(units_b = cmUnitsReplaceB)

angle_conv = field_dd %>%
filter(meas_type == 'thickness' & units_b == 'deg') %>%
select(c(unique_id, meas_a:units_b)) %>%
mutate(thick = meas_a * sin(meas_b * pi / 180))

thicknessReplace = replace(field_dd$meas_a, angle_conv$unique_id, angle_conv$thick)
field_ddd = field_dd %>% mutate(meas_a = thicknessReplace)
angleDelete = replace(field_dd$meas_b, angle_conv$unique_id, NA)
field_ddd = field_ddd %>% mutate(meas_b = angleDelete)
unitsDelete = replace(field_dd$units_b, angle_conv$unique_id, NA)
field_ddd = field_ddd %>% mutate(units_b = unitsDelete)

field_data = field_ddd

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
j = 1

for (i in dirs_to_iterate) {
	filelist = list.files(i)
	if (length(filelist > 0))
	{
		match = grep('*.shp',filelist)
		files[[j]] = st_read(file.path(i, filelist[match]))
		files[[j]] = st_transform(files[[j]], 2152) # transform into UTM coordinates for 12 N
		geom_tab = st_coordinates(files[[j]]) %>% as.data.frame()
		att_tab = as.data.frame(files[[j]]) %>% select(-geometry)
		att_tab$FID = att_tab$FID + 1
		attribute_tab[[j]] = att_tab %>% rename(id = 'FID', interpretations = 'LAYER') %>% mutate(id = id + (j * 1000))
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
mutate(id = as.factor(id))

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

# this function calculates the distance along the mapped linear feature.

barfaceDistFun = function(X, Y, Z) {
	
	x = cbind(X,Y,Z)
	
	dists = dist(x)
	
	disc = (1 - 4 * -2 * length(dists))

	if (disc > 0) {
		size = (-1 + sqrt(disc)) / 2
	} else {
		error('no solutions!')
	}

	face = numeric(size)

	j = 0

	for (i in 1:size) {
		face[i] = i + j
		j = j + size - i
	}

	barface_length = sum(dists[face])
	
	return(barface_length)
}

barFaceLengthMeas = bound_data %>% group_by(id) %>% 
summarize(meas_type = 'face_length', meas_a = barfaceDistFun(X,Y,Z)) 

measurements = bind_rows(barThicknessMeas, barFaceLengthMeas)

model_data = inner_join(no_geom, measurements, by = 'id')

# wrangle shapefile data so that it matches field_data format



# join up with the original data

comb_data = bind_rows(field_data, model_data)

#	export as saved data file

save(comb_data, file = here('data','derived_data', 'piceance_paleoslope_data.rda'))
