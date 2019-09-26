#!/usr/bin/env Rscript
#	i/o and processing for field data collected in the Piceance basin
#	Eric Barefoot
#	Sep 2019

#	load required packages

require(tidyverse)
require(here)
require(devtools)	# for the next thing
require(paleohydror)

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

saveRDS(field_data, file = here('data','derived_data', 'piceance_field_data.rds'))
