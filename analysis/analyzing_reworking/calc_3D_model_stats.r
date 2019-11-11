#!/usr/bin/env Rscript
#	preliminary analysis of data collected from 3D models from Piceance basin
#	Eric Barefoot
#	Oct 2019

# import libraries

library(broom)
library(purrr)
library(rsample)

# load data (3D model, field, and combined)

model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))

# first, lets just look at the bars, and calculate the relative abundance of fully preserved barforms.

just_bars_models = model %>% filter(interpretations %in% c('full','partial','truncated'))

model_boots = just_bars_models %>% nest(-formation) %>% mutate(boots = map(data, ~ bootstraps(., times = 4))) %>% unnest(boots)
 
calcFullFrac = function(split)
{
  analysis(split) %>% group_by(interpretations) %>% summarize(count = n())
}

boot_counts = model_boots %>% mutate(stats = map(splits, calcFullFrac))

percent_preserved = boot_counts %>% unnest(stats) %>% group_by(id, formation) %>% summarize(perc_full = count[interpretations == 'full'] / sum(count))

saveRDS(percent_preserved, file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

# percent_preserved %>% ggplot() + geom_step(aes(x = perc_full, y = ..density.., color = formation), stat = 'bin') 

# geom_histogram(aes(x = perc_full, fill = formation)) + facet_wrap(vars(formation), ncol = 1)
