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
boots = readRDS(file = here('data', 'derived_data', 'piceance_bar_preservation.rds'))

# bootstrap sampling already done in calc_3D_model_stats.r

format_function = function(pval, signif = 0.001) {
  if (pval > signif) {
    prep = paste('p =', signif(pval, 3))
  } else {
    prep = paste('p \\ll', signif)
  }
  return(prep)
}

# but maybe we want fewer resamples. modified script to only take four resamples.

source(here('analysis', 'analyzing_reworking', 'calc_3D_model_stats.r'))

barpres = model %>% filter(interpretations %in% c('full','partial','truncated')) %>%
transmute(formation = case_when(
  formation == 'molina' ~ 'molina',
  TRUE ~ 'bounding'),
  interp = case_when(
  interpretations == 'full' ~ 'full',
  TRUE ~ 'partial')
)

data_tab = table(barpres$interp, barpres$formation)

barpres_test = chisq.test(data_tab)

# quick hypothesis test since these are normally distributed.

# fit = aov(perc_full ~ formation, data = barpres)

# p_comp_reworking = summary(fit)[[1]][["Pr(>F)"]][1]

p_comp_reworking = format_function(barpres_test$p.value)

# if (barpres_test$p.value > 0.001) {
#   p_comp_reworking = paste('p =', barpres_test$p.value)
# } else {
#   p_comp_reworking = paste('p \\ll 0.001')
# }

# save(p_comp_reworking, file = here('data', 'data_summaries', 'reworking_pvalue.rda'))

barfractions = model %>% filter(interpretations %in% c('full','partial','truncated')) %>%
group_by(interpretations, formation) %>% summarize(n = n()) %>%
pivot_wider(values_from = n, names_from = interpretations) %>%
mutate(perc = 100 * full / (full + partial + truncated))

shirePerc = barfractions %>% filter(formation == 'shire') %>% pull(perc)
molinaPerc = barfractions %>% filter(formation == 'molina') %>% pull(perc)
atwellPerc = barfractions %>% filter(formation == 'atwell_gulch') %>% pull(perc)

boot_sds = boots %>% group_by(formation) %>% summarize(sdPres = sd(perc_full) * 100)

estimates = tibble(meanPres = c(shirePerc,molinaPerc,atwellPerc), formation = c('shire','molina','atwell_gulch'))

pres = inner_join(boot_sds, estimates, by = 'formation')

save(p_comp_reworking,
pres,
barpres_test,
file = here('data', 'data_summaries', 'paleo_pres.rda'))
