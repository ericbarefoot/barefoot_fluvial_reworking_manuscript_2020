#!/usr/bin/env Rscript
#	Making figures for group meeting presentation 18 Oct 2019
#	Eric Barefoot
#	Oct 2019

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

# establish directories

figdir = here('figures', 'outputs', 'group_meeting_20191018_interactions')

# load data

chamb_data = readRDS(file = here('data','derived_data', 'chamberlin_model_data.rds'))
chamb_model = readRDS(file = here('data','derived_data', 'chamberlin_stat_model.rds'))
model = readRDS(file = here('data','derived_data','piceance_3d_model_data.rds'))
field = readRDS(file = here('data','derived_data','piceance_field_data.rds'))
combined = readRDS(file = here('data', 'derived_data', 'piceance_field_model_data.rds'))
barpres = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

# for the section talking about how paleo geometries suck

cpal = hcl(h = c(80, 120, 240), c = rep(100, 3), l = c(85, 65, 20))

# mean_depth = function(split)
# {
#   analysis(split) %>% summarize(mean_depth = mean(meas_a, na.rm = TRUE))
# }
# 
field_depths = field %>% filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>% mutate(meas_a = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b)) 
# 
# field_depths_bootstraps = field_depths %>% nest(data = -formation) %>% mutate(boots = map(data, ~ bootstraps(., times = 200))) %>% unnest(boots)
# 
# depth_estimates_field = field_depths_bootstraps %>% mutate(depth = map(splits, mean_depth)) %>% 
# unnest(depth) %>% select(formation, mean_depth) %>% rownames_to_column(var = 'id')
# 
model_depths = model %>% filter(interpretations %in% c('full'))
# 
# model_boots = model_depths %>% nest(data = -formation) %>% mutate(boots = map(data, ~ bootstraps(., times = 200))) %>% unnest(boots)
# 
# depth_estimates_model = model_boots %>% mutate(depth = map(splits, mean_depth)) %>% unnest(depth) %>% select(formation, mean_depth) %>% rownames_to_column(var = 'id')
# 
# depth_estimates = inner_join(depth_estimates_field, depth_estimates_model, by = c('formation', 'id'), suffix = c('_field', "_model")) %>% select(-id) %>% gather(-formation, key = 'data_source', value = 'depth') %>% mutate(data_source = recode(.$data_source, mean_depth_field = 'field', mean_depth_model = 'model'))
# 
# ggplot(depth_estimates) + geom_step(aes(x = depth, y = ..density.., color = formation), stat = 'bin', bins = 20) + facet_wrap(vars(data_source), ncol = 1)

all_depths = bind_rows(field_depths, model_depths, .id = 'data_source') %>% mutate(data_source = recode(.$data_source, `1` = 'field', `2` = 'model'))

all_depths = all_depths %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 3, formation == 'atwell_gulch' ~ 2, formation == 'molina' ~ 1, formation == 'shire' ~ 0))

legend_ord = levels(with(all_depths, reorder(formation, formID)))

formation_labeller = as_labeller(c(`0` = 'Shire', `1` = 'Molina', `2` = 'Atwell Gulch'))
source_labeller = as_labeller(c('field' = 'Field Data', 'model' = '3D Model Data'))

depth_freq_field = ggplot() +
stat_bin(aes(x = meas_a, y = ..density.., color = reorder(formation, formID)), geom = 'step', position = "identity", bins = 20, size = 1.5, data = filter(all_depths, data_source == 'field')) +
theme_minimal() + 
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
labs(x = 'Flow Depth (m)', y = 'Probability Density') + facet_wrap(vars(formID), ncol = 1, labeller = formation_labeller) + 
theme(strip.background = element_blank(), strip.text.x = element_blank())

bedloadSamples = field %>% 
filter(meas_type == 'sample' & !(structure %in% c('outcrop', 'formation', 'bar_drape')) & meas_a < 2000) %>% 
mutate(D50 = meas_a * 1e-6) %>% 
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3, 
  formation == 'atwell_gulch' ~ 2, 
  formation == 'molina' ~ 1, 
  formation == 'shire' ~ 0))

legend_ord = levels(with(bedloadSamples, reorder(formation, formID)))

grainsize_freq = bedloadSamples %>% ggplot() + 
stat_bin(aes(x = D50 * 1e6, y = ..density.., color = reorder(formation, formID)), geom = 'step', position = "identity", binwidth = 65, size = 1.5) + 
theme_minimal() + 
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
labs(x = expression(paste('Bedload Sample ', D[50], ' (', mu, 'm)')), y = 'Probability Density') + 
facet_wrap(vars(formID), ncol = 1, labeller = formation_labeller) + 
theme(strip.background = element_blank(), strip.text.x = element_blank())

samples_for_slope = bedloadSamples %>% select(set_id, D50, textures, composition, samp_ind, samp_code, samp_description)
depths_for_slope = all_depths %>% filter(data_source == 'field') %>% select(-c(textures, composition, samp_ind, samp_code, samp_description))

paleohydroset = inner_join(samples_for_slope, depths_for_slope, by = 'set_id') %>% 
mutate(S = trampush_slp(D50, meas_a))

slope_freq = paleohydroset %>% ggplot() + stat_bin(aes(x = S, y = ..density.., color = reorder(formation, formID)), geom = 'step', position = 'identity', bins = 8, size = 1.5) +
scale_x_log10(breaks = c(1e-04, 2e-04, 3e-04, 5e-04, 0.001, 0.002, 0.003), labels = c(expression(10^-4), expression(2%*%10^-4), expression(3%*%10^-4), expression(5%*%10^-4), expression(10^-3), expression(2%*%10^-3), expression(3%*%10^-3))) +
theme_minimal() + 
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
labs(x = 'Slope (-)', y = 'Probability Density') + 
facet_wrap(vars(formID), ncol = 1)+ 
theme(strip.background = element_blank(), strip.text.x = element_blank())

ggsave(plot = slope_freq, filename = "petm_slope_piceance_frequency.png", path = figdir, width = 7, height = 7, units = "in")

ggsave(plot = grainsize_freq, filename = "petm_grainsize_piceance_frequency.png", path = figdir, width = 7, height = 7, units = "in")

ggsave(plot = depth_freq_field, filename = "petm_depth_piceance_frequency.png", path = figdir, width = 7, height = 7, units = "in")

# for the section introducing and explaining the basis of Ellen's work

avulsion_labeller = as_labeller(c('comp_avulsions' = 'Compensational', 'random_avulsions' = 'Random'))

# replica of her plot
chamb_replot = chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_point(aes(color = frac_input, shape = model)) + 
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Overbank \ncontribution', shape = 'Avulsion process') + 
theme_minimal() + xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller)
# theme(legend.position = 'bottom')

ggsave(plot = chamb_replot, filename = "chamb_replot.png", path = figdir, width = 7, height = 5, units = "in")

# replica of her plot
chamb_replot_trend = chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_smooth(se = FALSE, method = 'lm', color = 'red3') +
geom_point(aes(color = frac_input, shape = model)) + 
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Overbank \ncontribution', shape = 'Avulsion process') + 
theme_minimal() + xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller)
# theme(legend.position = 'bottom')

ggsave(plot = chamb_replot_trend, filename = "chamb_replot_trend.png", path = figdir, width = 7, height = 5, units = "in")

# plot showing vertical groupings of input ratio experiments
input_range = chamb_data %>% #filter(model == 'comp_avulsions') %>% 
ggplot(., aes(x = frac_input, y = frac_elems, shape = model)) + 
geom_point() +
labs(x = 'Overbank contribution', y = 'Fraction of fully preserved elements', shape = 'Avulsion process') + 
theme_minimal() +
scale_shape(labels = avulsion_labeller)

ggsave(plot = input_range, filename = "chamb_input.png", path = figdir, width = 7, height = 5, units = "in")

# plot showing vertical groupings of input ratio experiments - colorcoded
input_groups = chamb_data %>% # filter(model == 'comp_avulsions') %>% 
ggplot(., aes(x = frac_input, y = frac_elems)) + 
geom_point(aes(color = exp_set, shape = model)) + 
labs(x = 'Overbank contribution', y = 'Fraction of fully preserved elements', color = 'Experiment\nset', shape = 'Avulsion process') + 
theme_minimal()+
scale_shape(labels = avulsion_labeller)

ggsave(plot = input_groups, filename = "chamb_input_colored.png", path = figdir, width = 7, height = 5, units = "in")

# plot showing spread of all the groups
plot_exp_sets = chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_point(aes(color = exp_set, shape = model)) + 
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Experiment\nset', shape = 'Avulsion process') + 
theme_minimal() +
xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller)

ggsave(plot = plot_exp_sets, filename = "chamb_replot_colored.png", path = figdir, width = 7, height = 5, units = "in")

# plot showing trends of all the groups
plot_exp_sets_trends = chamb_data %>% 
ggplot(., aes(x = frac_elems, y = frac_supply)) + 
geom_smooth(aes(color = exp_set), method = 'lm', se = FALSE) +
geom_point(aes(color = exp_set, shape = model)) + 
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Experiment\nset', shape = 'Avulsion process') + 
theme_minimal() +
xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller)

ggsave(plot = plot_exp_sets_trends, filename = "chamb_replot_colored_trends.png", path = figdir, width = 7, height = 5, units = "in")

# plot showing the fixed effect fit
chamb_data_pred = chamb_data %>% mutate(pred_supply = predict(chamb_model))

chamb_fixed_effects_plot = ggplot(data = chamb_data_pred) + 
geom_line(aes(x = frac_elems, y = pred_supply, group = frac_input), color = 'red3') + 
geom_point(aes(x = frac_elems, y = frac_supply, color = frac_input, shape = model)) +
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Overbank\nContribution', shape = 'Avulsion process') + 
theme_minimal() +
xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller)

ggsave(plot = chamb_fixed_effects_plot, filename = "chamb_replot_fixed_effects_model.png", path = figdir, width = 7, height = 5, units = "in")

# for the section showing how to derive estimates of reworking

# import data for preservation from 3D models and plot the predicted values from the statistical model from Ellen's paper

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

d = barpres_pred_means

ann_data = tibble(x = c(rep(0,3), d$mean_elems), y = c(d$mean_supply, rep(0.2,3)), xend = rep(d$mean_elems,2), yend = rep(d$mean_supply,2), formation = rep(d$formation,2), formID = rep(d$formID,2))

chamb_piceance_overlay = ggplot() + 
geom_point(aes(x = frac_elems, y = frac_supply, shape = model), data = chamb_data) + 
geom_point(aes(x = frac_elems, y = frac_supply, color = reorder(formation, formID)), data = barpres_pred) + 
scale_x_continuous(expand = expand_scale(mult = c(0, .05))) + scale_y_continuous(expand = expand_scale(mult = c(0, .05))) + 
theme_minimal() +
xlim(0,1) + ylim(0.2,1) +
scale_shape(labels = avulsion_labeller) + 
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
scale_fill_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = 'Fraction of fully preserved elements', y = 'Fraction of sediment supply retained', color = 'Overbank\nContribution', shape = 'Avulsion process') 

ggsave(plot = chamb_piceance_overlay, filename = "chamb_piceance_overlay.png", path = figdir, width = 7, height = 5, units = "in")

chamb_piceance_overlay_ann = chamb_piceance_overlay + 
geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = reorder(formation, formID)), data = ann_data) +
geom_point(aes(x = mean_elems, y = mean_supply, fill = reorder(formation, formID)), size = 5, shape = 'circle filled', color = 'black', data = barpres_pred_means) 

ggsave(plot = chamb_piceance_overlay_ann, filename = "chamb_piceance_overlay_ann.png", path = figdir, width = 7, height = 5, units = "in")

petm_bar_preservation = ggplot() + 
geom_step(aes(x = frac_elems, y = ..density.., color = reorder(formation, formID)), data = barpres_pred, bins = 20, stat = 'bin', size = 1.5) + 
theme_minimal() + # xlim(0,1) + ylim(0.2,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = 'Fraction of fully preserved elements', y = 'Probability Density', color = 'Member') +
facet_wrap(vars(formID), ncol = 1) +
theme(strip.background = element_blank(), strip.text.x = element_blank())

ggsave(plot = last_plot(), filename = "petm_bar_preservation.png", path = figdir, width = 7, height = 7, units = "in")

petm_sediment_retention = ggplot() + 
geom_step(aes(x = frac_supply, y = ..density.., color = reorder(formation, formID)), data = barpres_pred, bins = 20, stat = 'bin', size = 1.5) + 
theme_minimal() + # xlim(0,1) + ylim(0.2,1) +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) +
labs(x = 'Fraction of sediment supply retained', y = 'Probability Density', color = 'Member') +
facet_wrap(vars(formID), ncol = 1) +
theme(strip.background = element_blank(), strip.text.x = element_blank())

ggsave(plot = last_plot(), filename = "petm_sediment_retention.png", path = figdir, width = 7, height = 7, units = "in")

depth_freq_compare_3d = ggplot(all_depths) + geom_step(aes(x = meas_a, y = ..count.., color = reorder(formation, formID)), stat = 'bin', bins = 20, size = 1.5) + 
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + 
facet_grid(vars(formID), vars(data_source), labeller = labeller(formID = formation_labeller, data_source = source_labeller)) + 
labs(x = 'Flow Depth (m)', y = 'Frequency') + 
theme_minimal()

ggsave(plot = depth_freq_compare_3d, filename = "petm_depth_piceance_compare_field_drone.png", path = figdir, width = 12, height = 7, units = "in")
