#!/usr/bin/env Rscript
#	Making figures for group meeting presentation 18 Oct 2019
#	Eric Barefoot
#	Oct 2019

# load packages

library(tidyverse)
library(here)
library(broom)
library(lme4)
if (interactive()) {
    require(colorout)
}

# load data

chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))

# for the section talking about how paleo geometries suck
# TODO: scavange plots from before about this, and incorporate new data.



# for the section introducing and explaining the basis of Ellen's work

# replica of her plot
chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_point(aes(color = frac_input, shape = model)) 

# plot showing vertical groupings of input ratio experiments
chamb_data %>% filter(model == 'comp_avulsions') %>% 
ggplot(., aes(x = frac_input, y = frac_elems)) + 
geom_point() 

# plot showing vertical groupings of input ratio experiments - colorcoded
chamb_data %>% 
ggplot(., aes(x = frac_input, y = frac_elems)) + 
geom_point(aes(color = exp_set)) 

# plot showing spread of all the groups
chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_point(aes(color = exp_set, shape = model)) 

# plot showing trends of all the groups
chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_smooth(aes(color = exp_set), method = 'lm', se = FALSE) +
geom_point(aes(color = exp_set, shape = model))

# plot showing the fixed effect fit
chamb_data_pred = chamb_data %>% mutate(pred_supply = predict(chamb_model))

ggplot(data = chamb_data_pred) + 
geom_line(aes(x = frac_elems, y = pred_supply, group = frac_input), color = 'red3') + 
geom_point(aes(x = frac_elems, y = frac_supply, color = frac_input, shape = model))

# for the section showing how to derive estimates of reworking

# import data for preservation from 3D models and plot the predicted values from the statistical model from Ellen's paper

barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

model = rep('comp_avulsions', nrow(barpres))
exp_set = rep('Ec', nrow(barpres))
frac_input = rnorm(nrow(barpres), 0.30, 0.0425)

barpres = bind_cols(barpres, tibble(model, exp_set, frac_input))

barpres = barpres %>% rename(frac_elems = 'perc_full') %>% ungroup()

barpres_pred = barpres %>% mutate(frac_supply = predict(chamb_model, newdata = .))

barpres_pred_means = barpres_pred %>% group_by(formation) %>% summarize(mean_supply = mean(frac_supply), mean_elems = mean(frac_elems), mean_input = mean(frac_input))

ggplot() + 
geom_point(aes(x = frac_elems, y = frac_supply, shape = model), data = chamb_data) + 
geom_point(aes(x = frac_elems, y = frac_supply, color = formation), data = barpres_pred) + 
geom_hline(aes(yintercept = mean_supply, color = formation), data = barpres_pred_means) +
geom_point(aes(x = mean_elems, y = mean_supply, fill = formation), size = 5, shape = 'circle filled', color = 'black', data = barpres_pred_means)

ggplot() + 
geom_step(aes(x = frac_supply, color = formation), data = barpres_pred, bins = 20, stat = 'bin')
