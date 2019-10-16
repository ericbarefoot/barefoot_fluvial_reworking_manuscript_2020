#!/usr/bin/env Rscript
#	i/o and processing for data digitized from Chamberlin 2019
#	Eric Barefoot
#	Oct 2019

#	load required packages

require(tidyverse)
require(here)
if (interactive()) {
    require(colorout)
}

# read in data

chamb_data = read_csv(file = here('data','raw_data','chamberlin_digitized.csv'))
chamb_greyscale = read_csv(file = here('data','raw_data','chamberlin_colorscale.csv'))

# # plot to demonstrate linearity of color mapping
# chamb_greyscale %>% ggplot(aes(x = greyscale, y = frac_input)) + geom_point() + geom_smooth(method = 'lm', se = FALSE)

color_model = lm(frac_input~greyscale, data = chamb_greyscale)

chamb_data = chamb_data %>% mutate(frac_input = predict(color_model, .))

# simple plot up of basic dataset.
# chamb_data %>%
# ggplot(., aes(x = frac_elems, y = frac_supply)) + 
# geom_point(aes(color = frac_input, shape = model))

# it looks like ellen did runs each with a constant value of overbank sediment input, and varied some other boundary condition to vary preserved elements. 

comp_coded = chamb_data %>% filter(model == 'comp_avulsions') %>% mutate(exp_set = 
case_when(
  between(frac_input, 0.55, 0.6) ~ 'Ac',
  between(frac_input, 0.46, 0.53) ~ 'Bc',
  between(frac_input, 0.39, 0.45) ~ 'Cc',
  between(frac_input, 0.3, 0.35) ~ 'Dc',
  between(frac_input, 0.2, 0.25) ~ 'Ec'
))

rand_coded = chamb_data %>% filter(model == 'random_avulsions') %>% mutate(exp_set = 
case_when(
  between(frac_input, 0.48, 0.53) ~ 'Ar',
  between(frac_input, 0.39, 0.45) ~ 'Br',
  between(frac_input, 0.32, 0.38) ~ 'Cr',
  between(frac_input, 0.25, 0.3) ~ 'Dr',
  between(frac_input, 0.1, 0.25) ~ 'Er'
))

# some plots to show this procedure and goals.
# 
# comp_coded %>% 
# ggplot(., aes(x = frac_input, y = frac_elems)) + 
# geom_point(aes(color = exp_set)) 
# 
# rand_coded %>% 
# ggplot(., aes(x = frac_input, y = frac_elems)) + 
# geom_point(aes(color = exp_set)) 
# 
# comp_coded %>% filter(!is.na(exp_set)) %>% 
# ggplot(., aes(x = frac_elems, y = frac_supply)) + 
# geom_smooth(aes(color = exp_set), method = 'lm', se = FALSE) + 
# geom_point(aes(color = exp_set)) 
# 
# rand_coded %>% filter(!is.na(exp_set)) %>% 
# ggplot(., aes(x = frac_elems, y = frac_supply)) + 
# geom_smooth(aes(color = exp_set), method = 'lm', se = FALSE) + 
# geom_point(aes(color = exp_set)) 

coded_data = bind_rows(comp_coded, rand_coded)
coded_data = coded_data %>% filter_all(all_vars(!is.na(.))) %>%
mutate(frac_elems = frac_elems / 100)

saveRDS(coded_data, file = here('data','derived_data', 'chamberlin_model_data.rds'))
