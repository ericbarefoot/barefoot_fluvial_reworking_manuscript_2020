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

if (interactive()) {
    require(colorout)
}

# load data

chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

# bootstrap sampling already done in calc_3D_model_stats.r

# but maybe we want fewer resamples. modified script to only take four resamples. 

source(here('analysis', 'analyzing_reworking', 'calc_3D_model_stats.r'))

barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

barpres = barpres %>% mutate(formation = as.factor(formation))

# quick hypothesis test since these are normally distributed.

fit = aov(perc_full ~ formation, data = barpres)

p_comp_reworking = summary(fit)[[1]][["Pr(>F)"]][1]

# save(p_comp_reworking, file = here('data', 'data_summaries', 'reworking_pvalue.rda'))
