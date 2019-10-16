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

# combine them together

combined = bind_rows(field, model)

# save new data file

saveRDS(combined, file = here('data','derived_data','piceance_field_model_data.rds'))


combined %>% filter(interpretations == 'full' | (structure == 'bar' & meas_type == 'thickness')) %>%
ggplot() + geom_histogram(aes(x = meas_a, fill = formation)) + facet_wrap(vars(formation), ncol = 1)


library(broom)
library(purrr)
library(rsample)


just_bars_models = model %>% filter(interpretations %in% c('full','partial','truncated'))

model_boots = bootstraps(just_bars_models, times = 200)

calcFullFrac = function(split)
{
  analysis(split) %>% group_by(formation, interpretations) %>% summarize(count = n())
}

boot_counts = model_boots %>% mutate(stats = map(splits, calcFullFrac))

percent_preserved = boot_counts %>% unnest(stats) %>% group_by(id, formation) %>% summarize(perc_full = count[interpretations == 'full'] / sum(count))

percent_preserved %>% ggplot() + geom_step(aes(x = perc_full, y = ..density.., color = formation), stat = 'bin') 

geom_histogram(aes(x = perc_full, fill = formation)) + facet_wrap(vars(formation), ncol = 1)
