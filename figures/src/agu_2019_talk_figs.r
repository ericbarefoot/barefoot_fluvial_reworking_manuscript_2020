#!/usr/bin/env Rscript
#	Making figures for AGU Fall Meeting presentation 9 Dec 2019
#	Eric Barefoot
#	Dec 2019

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
library(segmented)

if (interactive()) {
    require(colorout)
}

# establish directories

figdir = here('figures', 'outputs', 'agu_2019_talk')

# load data

chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))
strat_model_data = readRDS(here('data', 'derived_data', 'strat_model_results_20191202.rds'))

# Define a color palette for the stratigraphy.

cpal = hcl(h = c(80, 120, 240), c = rep(100, 3), l = c(85, 65, 20))

# Get data from field organized into tables for plotting

field_depths = field %>% filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>% mutate(meas_a = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b)) 

model_depths = model %>% filter(interpretations %in% c('full'))

all_depths = bind_rows(field_depths, model_depths, .id = 'data_source') %>% mutate(data_source = recode(.$data_source, `1` = 'field', `2` = 'model'))

all_depths = all_depths %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 3, formation == 'atwell_gulch' ~ 2, formation == 'molina' ~ 1, formation == 'shire' ~ 0))

legend_ord = levels(with(all_depths, reorder(formation, formID)))

formation_labeller = as_labeller(c(`0` = 'Shire', `1` = 'Molina', `2` = 'Atwell Gulch'))
source_labeller = as_labeller(c('field' = 'Field Data', 'model' = '3D Model Data'))

bedloadSamples = field %>% 
filter(meas_type == 'sample' & !(structure %in% c('outcrop', 'formation', 'bar_drape')) & meas_a < 2000) %>% 
mutate(D50 = meas_a * 1e-6) %>% 
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3, 
  formation == 'atwell_gulch' ~ 2, 
  formation == 'molina' ~ 1, 
  formation == 'shire' ~ 0))

legend_ord = levels(with(bedloadSamples, reorder(formation, formID)))

samples_for_slope = bedloadSamples %>% select(set_id, D50, textures, composition, samp_ind, samp_code, samp_description)
depths_for_slope = all_depths %>% filter(data_source == 'field') %>% select(-c(textures, composition, samp_ind, samp_code, samp_description))

paleohydroset = inner_join(samples_for_slope, depths_for_slope, by = 'set_id') %>% 
mutate(S = trampush_slp(D50, meas_a))

# modify and model data from suite of model runs.

strat_model_dataFilt = strat_model_data %>% 
filter(between(IR, 0.3, 0.55)) %>% 
filter(most_common_avulsion == 'compensational') %>% 
mutate(datatype = 'Model Run')

strat_model_dataFilt2 = strat_model_dataFilt %>% mutate(presInv = 1 / (pres_percents-1.01), n2gInv = 1/n2g)

n2gMod = lm(n2gInv ~ presInv + latmob + IR, data = strat_model_dataFilt2, weight = 1-pres_percents)

n2gModSimple = lm(n2gInv ~ presInv, data = strat_model_dataFilt2, weight = 1-pres_percents)

reworkingMod = lm(reworking ~ pres_percents + latmob + IR, data = strat_model_dataFilt2)
reworkingModSeg = segmented(reworkingMod, seg.Z = ~pres_percents)

reworkingModSimple = lm(reworking ~ pres_percents, data = strat_model_dataFilt2)
reworkingModSegSimple = segmented(reworkingModSimple, seg.Z = ~pres_percents)

strat_model_dataFilt3 = strat_model_dataFilt2 %>% mutate(n2gPred = 1/predict(n2gModSimple)) %>% 
mutate(reworkingPred = predict(reworkingModSegSimple))

avul = rep('compensational', nrow(barpres))
IR = rnorm(nrow(barpres), 0.4, 0.05)
latmob = rnorm(nrow(barpres), 8, 2)

barpres2 = bind_cols(barpres, tibble(avul, latmob, IR))

barpres3 = barpres2 %>% rename(pres_percents = 'perc_full') %>%
mutate(formation = recode(formation, shire = 'atwell_gulch', atwell_gulch = 'shire'), presInv = 1 / (pres_percents-1.01)) %>% ungroup()

barpres_pred = barpres3 %>% mutate(n2g = 1/predict(n2gMod, newdata = .), reworking = predict(reworkingModSeg, newdata = .)) %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3, 
  formation == 'atwell_gulch' ~ 2, 
  formation == 'molina' ~ 1, 
  formation == 'shire' ~ 0))

barpres_pred_means = barpres_pred %>% group_by(formation) %>% summarize(mean_reworking = mean(reworking), mean_elems = mean(pres_percents)) %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3, 
  formation == 'atwell_gulch' ~ 2, 
  formation == 'molina' ~ 1, 
  formation == 'shire' ~ 0))

d = barpres_pred_means

ann_data = tibble(x = c(rep(0,3), d$mean_elems), y = c(d$mean_reworking, rep(0.2,3)), xend = rep(d$mean_elems,2), yend = rep(d$mean_reworking,2), formation = rep(d$formation,2), formID = rep(d$formID,2))

#######################################################################

start_theme = theme_get()

pres_theme = theme_update(#line = element_line(color = '#bdbdbd'),
panel.background = element_rect(fill = '#bdbdbd', color = '#bdbdbd'), 
legend.key = element_rect(fill = '#bdbdbd', color = '#bdbdbd'), 
panel.grid = element_line(color = '#a79c99'), 
strip.background = element_blank(), 
strip.text.x = element_blank(), 
plot.background = element_rect(fill = "transparent", color = 'transparent'), 
legend.background = element_rect(fill = "transparent"))

depth_freq_field = ggplot() +
stat_bin(aes(x = meas_a, y = ..density.., fill = reorder(formation, formID)), geom = 'bar', position = "identity", bins = 20, data = filter(all_depths, data_source == 'field'), color = 'black', size = 0.5) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
labs(x = 'Flow Depth (m)', y = 'Probability Density') + 
facet_wrap(vars(formID), ncol = 1, labeller = formation_labeller)

slope_freq = paleohydroset %>% ggplot() + stat_bin(aes(x = S, y = ..density.., fill = reorder(formation, formID)), geom = 'bar', position = 'identity', bins = 8, color = 'black', size = 0.5) +
scale_x_log10(breaks = c(1e-04, 2e-04, 5e-04, 0.001, 0.003), labels = c(expression(10^-4), expression(2%*%10^-4), expression(5%*%10^-4), expression(10^-3), expression(3%*%10^-3))) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
labs(x = 'Slope (-)', y = 'Probability Density') + 
facet_wrap(vars(formID), ncol = 1)

plot_width = 5

ggsave(plot = slope_freq, filename = "petm_slope_piceance_frequency.png", path = figdir, width = plot_width, height = plot_width, units = "in", bg = 'transparent')

ggsave(plot = depth_freq_field, filename = "petm_depth_piceance_frequency.png", path = figdir, width = plot_width, height = plot_width, units = "in", bg = 'transparent')

strat_model_runs = ggplot(aes(x = pres_percents, y = reworking), data = strat_model_dataFilt3) + 
geom_point() + 
xlim(0,1) + ylim(0.4,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = '% fully preserved bars', y = '% undisturbed stratigraphy') 

ggsave(plot = strat_model_runs, filename = "strat_model_runs.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

strat_model_runs_trend = strat_model_runs + 
geom_line(aes(x = pres_percents, y = reworkingPred), data = strat_model_dataFilt3, color = 'red3', size = 1.5)

ggsave(plot = strat_model_runs_trend, filename = "strat_model_runs_trend.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

strat_model_runs_overlay = strat_model_runs +
geom_point(aes(x = pres_percents, y = reworking, color = formation), data = barpres_pred) + 
theme(legend.position="none")

ggsave(plot = strat_model_runs_overlay, filename = "strat_model_runs_overlay.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

strat_model_runs_overlay_legend = strat_model_runs + 
geom_point(aes(x = pres_percents, y = reworking, color = formation), data = barpres_pred) 

ggsave(plot = strat_model_runs_overlay_legend, filename = "strat_model_runs_overlay_legend.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

bootstrap_plot = ggplot(aes(y = ..density.., fill = reorder(formation, formID)), data = barpres_pred) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
facet_wrap(vars(formID), ncol = 1) 

petm_bar_preservation = bootstrap_plot + 
stat_bin(aes(x = pres_percents, color = reorder(formation, formID)), bins = 20, geom = 'bar', size = 0.5, alpha = 0.4) + 
geom_vline(aes(xintercept = mean_elems, color = reorder(formation, formID)), data = d, size = 2) +
labs(x = '% fully preserved bars', y = 'Probability Density', color = 'Member')

ggsave(plot = petm_bar_preservation, filename = "petm_bar_preservation.png", path = figdir, width = plot_width*1.3, height = plot_width, units = "in", bg = 'transparent')

petm_sediment_retention = bootstrap_plot + 
stat_bin(aes(x = reworking, color = reorder(formation, formID)), data = barpres_pred, bins = 20, geom = 'bar', size = 0.5, alpha = 0.4) + 
geom_vline(aes(xintercept = mean_reworking, color = reorder(formation, formID)), data = d, size = 2) +
labs(x = '% undisturbed stratigraphy', y = 'Probability Density', color = 'Member')

ggsave(plot = last_plot(), filename = "petm_sediment_retention.png", path = figdir, width = plot_width*1.3, height = plot_width, units = "in", bg = 'transparent')

strat_model_runs_n2g = ggplot(aes(x = pres_percents, y = n2g), data = strat_model_dataFilt3) + 
geom_point() + 
xlim(0,1) + ylim(0,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = '% fully preserved bars', y = '% Sand') 

ggsave(plot = strat_model_runs_n2g, filename = "strat_model_runs_n2g.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

strat_model_runs_n2g_trend = strat_model_runs_n2g + 
geom_line(aes(x = pres_percents, y = n2gPred), data = strat_model_dataFilt3, color = 'red3', size = 1.5)

ggsave(plot = strat_model_runs_n2g_trend, filename = "strat_model_runs_n2g_trend.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

############### BONUS SLIDES ###############


strat_model_runs_n2g_overlay = strat_model_runs_n2g + 
geom_point(aes(x = pres_percents, y = n2g, color = formation), data = barpres_pred) + 
theme(legend.position="none")

ggsave(plot = strat_model_runs_n2g_overlay, filename = "strat_model_runs_n2g_overlay.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')

stratFilt = strat_model_data %>% 
filter(between(IR, 0.3, 0.55))

reworking_by_avulsion_type = ggplot(aes(x = pres_percents, y = reworking, color = most_common_avulsion), data = stratFilt) + 
geom_point() + 
geom_line(aes(x = pres_percents, y = reworkingPred), data = strat_model_dataFilt3, color = '#009e2f', size = 1.5) + 
xlim(0,1) + ylim(0.4,1) +
geom_abline() + 
labs(x = '% fully preserved bars', y = '% undisturbed stratigraphy', color = 'Avulsion Rule') 

ggsave(plot = reworking_by_avulsion_type, filename = "reworking_by_avulsion_type.png", path = figdir, width = 7, height = 5, units = "in", bg = 'transparent')
