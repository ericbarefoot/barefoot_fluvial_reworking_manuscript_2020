## scripts and utilities to read in csv from several different runs of strat model
## Eric Barefoot
## Dec 2019

# load libraries

# library(colorout)
library(tidyverse) 
library(here)

# load data

epc_data_1 = read_csv(here('data', 'raw_data', 'strat_model_data', 'out_tab_20191202_1.csv'))
epc_data_2 = read_csv(here('data', 'raw_data', 'strat_model_data', 'out_tab_20191202_2.csv'))
epc_data_3 = read_csv(here('data', 'raw_data', 'strat_model_data', 'out_tab_20191202_3.csv'))
epc_data_4a = read_csv(here('data', 'raw_data', 'strat_model_data', 'out_tab_20191202_4a.csv'))
epc_data_4b = read_csv(here('data', 'raw_data', 'strat_model_data', 'out_tab_20191202_4b.csv'))

epc_data = bind_rows(epc_data_1, epc_data_2, epc_data_3, epc_data_4a, epc_data_4b, .id = 'chunk')

epc_data = epc_data %>% filter(reworking < 1) %>% 
mutate(avul = recode(avul, `0` = 'clustered', `1` = 'compensational', `2` = 'random')) %>% 
mutate(most_common_avulsion = recode(most_common_avulsion, `0` = 'clustered', `1` = 'compensational', `2` = 'random')) %>% 
mutate(pres_percents = pres_percents / 100)

saveRDS(epc_data, here('data', 'derived_data', 'strat_model_results_20191202.rds'))

# This section is for when I want to investigate the variability of the preserved elements estimator. 

epc_pres_elem_1 = read_csv(here('data', 'raw_data', 'strat_model_data', 'elems_tab_20191202_1.csv'))
epc_pres_elem_2 = read_csv(here('data', 'raw_data', 'strat_model_data', 'elems_tab_20191202_2.csv'))
epc_pres_elem_3 = read_csv(here('data', 'raw_data', 'strat_model_data', 'elems_tab_20191202_3.csv'))
epc_pres_elem_4a = read_csv(here('data', 'raw_data', 'strat_model_data', 'elems_tab_20191202_4a.csv'))
epc_pres_elem_4b = read_csv(here('data', 'raw_data', 'strat_model_data', 'elems_tab_20191202_4b.csv'))

epc_pres_elem = bind_rows(epc_pres_elem_1, epc_pres_elem_2, epc_pres_elem_3, epc_pres_elem_4b, epc_pres_elem_4b, .id = 'chunk')

saveRDS(epc_pres_elem, here('data', 'derived_data', 'strat_model_results_preservation_20191202.rds'))
