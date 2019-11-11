# script to test for differences in reworking between molina and other members
# Eric Barefoot
# Nov 2019

# load libraries

library(tidyverse)
library(here)
library(broom)
library(lme4)
library(devtools)  # for the next thing
library(colorspace)
library(paleohydror)
library(purrr)
library(rsample)

library(conover.test)

if (interactive()) {
    require(colorout)
}

# load data

# chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
# chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
# barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

field_depths = field %>% filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>% mutate(meas_a = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b), formation = as.factor(formation)) 

# quick hypothesis test since these are normally distributed.

depths = field_depths$meas_a
formations = field_depths$formation

fit = conover.test(depths, formations, method = 'bonferroni', table = 'false', kw = 'false')

p_comp_AG_M = fit$P[1]
p_comp_AG_S = fit$P[2]
p_comp_S_M = fit$P[3]


model_depths = model %>% filter(interpretations %in% c('full'))

all_depths = bind_rows(field_depths, model_depths, .id = 'data_source') %>% mutate(data_source = recode(.$data_source, `1` = 'field', `2` = 'model'))

all_depths = all_depths %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 3, formation == 'atwell_gulch' ~ 2, formation == 'molina' ~ 1, formation == 'shire' ~ 0))

bedloadSamples = field %>% 
filter(meas_type == 'sample' & !(structure %in% c('outcrop', 'formation', 'bar_drape')) & meas_a < 2000) %>% 
mutate(D50 = meas_a * 1e-6) %>% 
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3, 
  formation == 'atwell_gulch' ~ 2, 
  formation == 'molina' ~ 1, 
  formation == 'shire' ~ 0))

samples_for_slope = bedloadSamples %>% select(set_id, D50, textures, composition, samp_ind, samp_code, samp_description)
depths_for_slope = all_depths %>% filter(data_source == 'field') %>% select(-c(textures, composition, samp_ind, samp_code, samp_description))

paleohydroset = inner_join(samples_for_slope, depths_for_slope, by = 'set_id') %>% 
mutate(S = trampush_slp(D50, meas_a))

# slope tests

slopes = paleohydroset$S
formations = paleohydroset$formation

fit_slope = conover.test(slopes, formations, method = 'bonferroni', table = 'false', kw = 'false')

p_comp_AG_M_slope = fit_slope$P[1]
p_comp_AG_S_slope = fit_slope$P[2]
p_comp_S_M_slope = fit_slope$P[3]

# save(p_comp_reworking, file = here('data', 'data_summaries', 'paleo_geom.rda'))
