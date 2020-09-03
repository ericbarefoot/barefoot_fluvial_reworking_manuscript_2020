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
library(purrr)
library(rsample)

load_all('/Users/ericfoot/Dropbox/programs/git_repos/paleohydror')

library(conover.test)

if (interactive()) {
    require(colorout)
}

format_function = function(pval, signif = 0.001) {
  if (pval > signif) {
    prep = paste('p =', signif(pval, 3))
  } else {
    prep = paste('p \\ll', signif)
  }
  return(prep)
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
# fit = conover.test(depths, formations, method = 'bonferroni')

p_comp_AG_M = format_function(fit$P[1])
p_comp_AG_S = format_function(fit$P[2])
p_comp_S_M = format_function(fit$P[3])

depth_summary = field_depths %>% group_by(formation) %>% summarize(meanDepth = mean(meas_a), sdDepth = sd(meas_a))

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
mutate(S = trampush_slp(D50, meas_a), tau = (1000 * 9.8 * meas_a * S), taustar = tau / (1065 * 9.8 * D50), Rep = Rep(D50, (1650 / 1000), 9.8, 1.004e-6))

# slope tests

slopes = paleohydroset$S
formations = paleohydroset$formation

fit_slope = conover.test(slopes, formations, method = 'bonferroni', table = 'false', kw = 'false')
# fit_slope = conover.test(slopes, formations, method = 'bonferroni')
# fit_tau = conover.test(paleohydroset$taustar, formations, method = 'bonferroni')

p_comp_AG_M_slope = format_function(fit_slope$P[1])
p_comp_AG_S_slope = format_function(fit_slope$P[2])
p_comp_S_M_slope = format_function(fit_slope$P[3])

slope_summary = paleohydroset %>% group_by(formation) %>% summarize(meanSlope = mean(S), sdSlope = sd(S))

tau_summary = paleohydroset %>% group_by(formation) %>% summarize(meanTau = mean(taustar), sdTau = sd(taustar))

geoms = inner_join(slope_summary, depth_summary, by = 'formation')
geoms = inner_join(geoms, tau_summary, by = 'formation')

# calculating abundance of different avulsion types

avuls = field %>% filter(structure == 'avulsion') %>%
# avuls = field %>% filter(set_id == 'EHAJEKA', structure == 'avulsion') %>%
# avuls = field %>% filter(set_id == 'FORMANA', structure == 'avulsion') %>%
select(formation, interpretations) %>%
transmute(formation = case_when(
  formation == 'molina' ~ 'molina',
  TRUE ~ 'bounding'),
  interp = interpretations
)

avultest = chisq.test(table(avuls))

p_comp_avulsion = format_function(avultest$p.value)

# find out how many new data points we have.

ndata = field_depths %>% mutate(mine = ifelse(set_id == 'FORMANA', 'no', 'yes')) %>% group_by(mine) %>% summarize(n = n()) %>% pull(n)

nold = ndata[1]
nnew = ndata[2]

save(nold,
  nnew,
  p_comp_AG_M_slope,
  p_comp_AG_S_slope,
  p_comp_S_M_slope,
  p_comp_AG_M,
  p_comp_AG_S,
  p_comp_S_M,
  p_comp_avulsion,
  avultest,
  geoms,
  file = here('data', 'data_summaries', 'paleo_geom.rda'))
