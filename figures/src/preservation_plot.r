# function(s) to make a plot of data from 3D model mapping data.
# specifically the data on preservation relating it to bypass
# Eric Barefoot
# Sep 2019

# load packages

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

figdir = here('figures', 'outputs', 'manuscript_figures')

# load data

chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

cpal = hcl(h = c(80, 120, 240), c = rep(100, 3), l = c(85, 65, 20))

model = rep('comp_avulsions', nrow(barpres))
exp_set = rep('Ec', nrow(barpres))
frac_input = rnorm(nrow(barpres), 0.30, 0.0425)

barpres = bind_cols(barpres, tibble(model, exp_set, frac_input))

barpres = barpres %>% rename(frac_elems = 'perc_full') %>% ungroup()

barpres_pred = barpres %>% mutate(frac_supply = predict(chamb_model, newdata = .)) %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3,
  formation == 'atwell_gulch' ~ 2,
  formation == 'molina' ~ 1,
  formation == 'shire' ~ 0))

barpres_pred_means = barpres_pred %>% group_by(formation) %>% summarize(mean_supply = mean(frac_supply), mean_elems = mean(frac_elems), mean_input = mean(frac_input)) %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3,
  formation == 'atwell_gulch' ~ 2,
  formation == 'molina' ~ 1,
  formation == 'shire' ~ 0))

legend_ord = levels(with(barpres_pred, reorder(formation, formID)))
legend_ord = levels(with(barpres_pred_means, reorder(formation, formID)))


petm_bar_preservation = ggplot() +
geom_step(aes(x = frac_elems, y = ..density.., color = reorder(formation, formID)), data = barpres_pred, bins = 20, stat = 'bin', size = 1.5) +
theme_minimal() + # xlim(0,1) + ylim(0.2,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = 'Fraction of fully preserved bars', y = 'Probability Density', color = 'Member') +
facet_wrap(vars(formID), ncol = 1) +
theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.position = 'none')


petm_sediment_retention = ggplot() +
geom_step(aes(x = frac_supply, y = ..density.., color = reorder(formation, formID)), data = barpres_pred, bins = 20, stat = 'bin', size = 1.5) +
theme_minimal() + # xlim(0,1) + ylim(0.2,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = 'Fraction of sediment supply retained', y = NULL, color = 'Member') +
facet_wrap(vars(formID), ncol = 1) +
theme(strip.background = element_blank(), strip.text.x = element_blank())

w1 = 2.6
w2 = 3.7
h1 = 3.1

ggsave(plot = petm_bar_preservation, filename = "petm_bar_preservation.png", path = figdir, width = w1, height = h1, units = "in")

ggsave(plot = petm_sediment_retention, filename = "petm_sediment_retention.png", path = figdir, width = w2, height = h1, units = "in")
