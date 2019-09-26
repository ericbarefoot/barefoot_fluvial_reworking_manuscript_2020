#!/usr/bin/env Rscript
#	i/o and processing for field data collected in the Piceance basin
#	Eric Barefoot
#	Sep 2019

#	load required packages

require(tidyverse)
require(here)

# build datasets TODO: put in logic to control this.

# source(here('analysis','data_io','piceance_3d_mapping_io.r'))
# source(here('analysis','data_io','piceance_field_io.r'))

# load data

model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))

combined = bind_rows(model, field)

saveRDS(combined, file = here('data','derived_data','piceance_field_model_data.rds'))
