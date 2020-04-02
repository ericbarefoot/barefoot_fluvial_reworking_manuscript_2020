# function(s) to make a plot of data from 3D model mapping data and field data.
# specifically the data on kinematics vs geometric estimations
# Eric Barefoot
# Mar 2020

# first, load libraries

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
barpres_raw = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

barpres = barpres_raw %>%
rename(frac_elems = 'perc_full') %>% ungroup() %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3,
  formation == 'atwell_gulch' ~ 2,
  formation == 'molina' ~ 1,
  formation == 'shire' ~ 0)
) %>%
select(formID, id, frac_elems) %>%
mutate(depth = NA, slope = NA) %>%
slice(1:round(n() * 0.5))

barpres_means = barpres %>% group_by(formID) %>%
summarize(mean_elems = median(frac_elems)) %>%
mutate(depth = NA, slope = NA)

cpal = hcl(h = c(80, 120, 240), c = rep(100, 3), l = c(85, 65, 20))

field_depths = field %>%
filter(structure %in% c('bar','channel') & meas_type %in% c('thickness', 'dimensions')) %>%
mutate(depth = case_when(meas_type == 'thickness' ~ meas_a, meas_type == 'dimensions' ~ meas_b)) %>%
select(formation, set_id, depth) %>%
mutate(formID = case_when(
  formation == 'ohio_creek' ~ 3,
  formation == 'atwell_gulch' ~ 2,
  formation == 'molina' ~ 1,
  formation == 'shire' ~ 0)
)

bedloadSamples = field %>%
filter(meas_type == 'sample' & !(structure %in% c('outcrop', 'formation', 'bar_drape')) & meas_a < 2000) %>%
transmute(set_id = set_id, D50 = meas_a * 1e-6)

slopes =
inner_join(field_depths, bedloadSamples, by = 'set_id') %>%
group_by(formID, set_id) %>% summarize(depth = median(depth), D50 = median(D50)) %>%
mutate(slope = atan(trampush_slp(D50, depth)) * (90 / pi)) %>%
select(formID, slope) %>%
mutate(id = NA, frac_elems = NA, mean_elems = NA, depth = NA)

depths = field_depths %>% select(formID, depth) %>%
mutate(id = NA, frac_elems = NA, mean_elems = NA, slope = NA)

data_table = bind_rows(depths, slopes, barpres, barpres_means, .id = 'split') %>%
mutate(formID = as.factor(formID)) %>% select(-id) %>%
pivot_longer(c(depth, slope, frac_elems, mean_elems), names_to = 'meas', values_to = 'vals') %>%
mutate(meas = ifelse(meas == 'mean_elems', 'frac_elems', meas)) %>%
mutate(meas = as.factor(
  case_when(
  meas == 'depth' ~ 100,
  meas == 'slope' ~ 200,
  meas == 'frac_elems' ~ 300)
)
)

labels = as_labeller(c(`0` = 'Shire', `1` = 'Molina', `2` = 'Atwell Gulch', `100` = 'Flow Depth (m)', `200` = 'Fluvial Slope (º)', `300` = '% Fully Preserved Bars'))
form_labs = c('Shire','Molina','Atwell Gulch')

data_figure = ggplot() +
# plot depths
stat_bin(data = filter(data_table, split == 1),
aes(x = vals, y = ..count.., fill = formID),
geom = 'bar', position = "identity", bins = 20, color = 'black', size = 0.5) +
# plot slopes
stat_bin(data = filter(data_table, split == 2),
aes(x = vals, y = ..count.., fill = formID),
geom = 'bar', position = "identity", bins = 7, color = 'black', size = 0.5) +
# plot histogram of % bar estimates from bootstrapping
stat_bin(data = filter(data_table, split == 3),
aes(x = vals, y = ..count.., fill = formID, color = formID),
geom = 'bar', position = "identity", bins = 20, size = 0.5, alpha = 0.4) +
# plot mean line
geom_vline(data = filter(data_table, split == 4),
aes(xintercept = vals, color = formID), size = 2) +
# adjust other parameters for the whole plot.
scale_fill_manual(values = cpal, name = 'Stratigraphic Member', labels = form_labs) +
scale_color_manual(values = cpal) +
facet_grid(rows = vars(formID), cols = vars(meas), scales = 'free_x', labeller = labels, switch = 'x') +
labs(x = NULL, y = 'Frequency') +
theme_minimal() +
theme(strip.background = element_blank(), legend.position = 'none', strip.placement = "outside") +
guides(color = 'none')

# save figure

ggsave(plot = data_figure, filename = "petm_data.pdf", path = figdir, width = 6, height = 4, units = "in")
ggsave(plot = data_figure, filename = "petm_data.png", path = figdir, width = 6, height = 4, units = "in")
