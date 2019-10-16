#!/usr/bin/env Rscript
#	Building a statistical model for predicting sediment reworking from fraction of preserved bar structures.
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

# for helpfulness later - plot up the data with factors, etc
# factors are that the model results were generated by either:
  # (1) compensational avulsions or
  # (2) random avulsions.
# The other factors shown are sets of experiments with the same input overbank/channel material ratio ("input ratio")
# This is a strong controlling factor, and is also mapped to a frac_input continuous variable. 

# chamb_data %>%
# ggplot(., aes(x = frac_elems, y = frac_supply)) + 
# geom_point(aes(color = exp_set, shape = model)) 

# We are interested in finding the best predictive model for these data. 
# Let's try out some. 

# null
# fixed-effects model with
# fixed effects on % fully preserved structures

fixed_model_null = lm(frac_supply ~ frac_elems, data = coded_data)

# fixed-effects model with
# fixed effects on % fully preserved structures
# fixed effects on input ratio

fixed_model_input = lm(frac_supply ~ frac_elems + frac_input, data = coded_data)

# mixed-effects model with
# fixed effects on % fully preserved structures
# random effects on experiment set (proxy for input ratio) mean

mixed_model_int_expSet = lmer(frac_supply ~ frac_elems + (1|exp_set), data = coded_data)

# mixed-effects model with
# fixed effects on % fully preserved structures
# random effects on experiment set (proxy for input ratio) mean
# random effects on avulsion type mean

mixed_model_int_expSet_avulType = lmer(frac_supply ~ frac_elems + (1|exp_set) + (1|model), data = coded_data)

# mixed-effects model with 
# fixed effects on % fully preserved structures
# random effects on experiment set (proxy for input ratio) mean and slope

mixed_model_slp_expSet = lmer(frac_supply ~ frac_elems + (1 + frac_elems|exp_set), data = coded_data, REML = FALSE, control = lmerControl(optimizer = 'Nelder_Mead'))

# mixed-effects model with
# fixed effects on % fully preserved structures
# random effects on experiment set (proxy for input ratio) mean and slope
# random effects on avulsion type mean and slope

mixed_model_slp_expSet_avulType = lmer(frac_supply ~ frac_elems + (1+frac_elems|exp_set) + (1+frac_elems|model), data = coded_data, REML = FALSE)

# mixed-effects model with 
# fixed effects on % fully preserved structures
# fixed effects on input ratio
# random effects on experiment set (proxy for input ratio) mean and slope
# random effects on avulsion type mean and slope

mixed_model_input_slp_expSet_avulType = lmer(frac_supply ~ frac_elems + frac_input + (1+frac_elems|model), data = coded_data, REML = FALSE)

# So after doing all of those models, we need to evaluate them. 

models_aic = AIC(fixed_model_null, fixed_model_input, mixed_model_int_expSet, mixed_model_int_expSet_avulType, mixed_model_slp_expSet, mixed_model_slp_expSet_avulType, mixed_model_input_slp_expSet_avulType)

# based on this, it looks like the best model is really the simplest fixed-effects model with the input ratio included. 
# That effect seems to be pretty important! Since we have NO way of constraining that in the ancient, the best route forward is determine the sensitivity of the analysis to overbank vs channel deposition and work that in somehow.

# for fun, plot this model on the data.

# chamb_data_pred = chamb_data %>% mutate(pred_supply = predict(fixed_model_input))
# 
# ggplot(data = chamb_data_pred) + 
# geom_point(aes(x = frac_elems, y = frac_supply, color = frac_input, shape = model)) + 
# geom_point(aes(x = frac_elems, y = pred_supply), color = 'red3')

# export the best model as a derived data product. 

saveRDS(fixed_model_input, file = here('data','derived_data', 'chamberlin_stat_model.rds'))