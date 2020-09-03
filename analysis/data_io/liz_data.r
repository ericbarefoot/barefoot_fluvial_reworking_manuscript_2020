#!/usr/bin/env Rscript
#	i/o and processing for field data collected in the Piceance basin
#	Eric Barefoot
#	Sep 2019

#	load required packages

require(tidyverse)
require(here)
require(devtools)	# for the next thing
require(paleohydror)
library(broom)
library(purrr)
library(rsample)

field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))

avulsions = field %>% filter(structure == 'avulsion')

avulsion_boots = avulsions %>% select(formation, interpretations) %>% nest(data = c(interpretations)) %>% mutate(boots = map(data, ~ bootstraps(., times = 100))) %>% unnest(boots)

calcFullFrac = function(split)
{
  analysis(split) %>% group_by(interpretations) %>% summarize(count = n())
}

boot_counts = avulsion_boots %>% mutate(stats = map(splits, calcFullFrac))

percent_transitional = boot_counts %>% unnest(stats) %>% group_by(id, formation) %>% summarize(perc_trans = count[interpretations == 'transitional'] / sum(count))

avuls_summary = percent_transitional %>% group_by(formation) %>%
summarize(meanAvul = mean(perc_trans), sdAvul = sd(perc_trans))

save(avuls_summary,
file = here('data', 'data_summaries', 'avuls_summary.rda'))

saveRDS(percent_transitional, file = here('data', 'derived_data', 'piceance_avulsion_style.rds'))
