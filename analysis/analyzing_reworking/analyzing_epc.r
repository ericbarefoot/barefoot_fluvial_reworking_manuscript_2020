## visualizing outputs of strat model
## Eric Barefoot
## Dec 2019

# load libraries

library(colorout)
library(tidyverse) 
library(here)

# load data

epc_data = readRDS(here('data', 'derived_data', 'strat_model_results_20191202.rds'))

ggplot(data = epc_data, aes(x = pres_percents, y = reworking)) + 
geom_abline() + geom_point(aes(color = avul))

epc_data %>%  filter(most_common_avulsion == 'clustered') %>% 
ggplot(aes(x = pres_percents, y = reworking)) + 
geom_abline() + geom_point()

n2g_pres_plot = epc_data %>%  filter(most_common_avulsion == 'clustered') %>% 
filter(IR == 0.4) %>% 
ggplot(data = ., aes(x = pres_percents, y = n2g)) 

n2g_pres_plot + geom_point(aes(color = vertagg_rate)) 
# n2g_pres_plot + geom_point() 

nagg_pres_plot = epc_data %>% filter(most_common_avulsion == 'clustered') %>% 
# filter(latmob == 4) %>% 
ggplot(data = ., aes(x = pres_percents, y = nagg)) 

# nagg_pres_plot + geom_point() 
nagg_pres_plot + geom_point(aes(color = IR)) 

n2g_nagg_plot_log = epc_data %>% # filter(avul == 'random') %>% 
ggplot(data = ., aes(x = log(nagg), y = log(n2g)))

n2g_nagg_plot_log + geom_point(aes(color = avul), alpha = 0.4) + 
geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, aes(color = avul,group = avul))

n2g_nagg_plot = epc_data %>% # filter(avul == 'clustered') %>% 
ggplot(data = ., aes(x = nagg, y = n2g))

n2g_nagg_plot + geom_point(aes(color = avul), alpha = 0.4)

nagg_vagg_plot = epc_data %>% 
# filter(avul == 'compensational') %>% 
# filter(IR == 0.4) %>% 
ggplot(data = ., aes(x = vertagg_rate, y = nagg))

nagg_vagg_plot + geom_jitter(aes(color = IR))

for (i in seq(0.4,0.9,0.1)) {
  n2g_vagg_plot = epc_data %>% 
  filter(avul == 'compensational') %>% 
  filter(IR == round(i,1)) %>% 
  ggplot(data = ., aes(x = log(vertagg_rate), y = log(n2g)))
  # ggplot(data = ., aes(x = vertagg_rate, y = n2g))

  print(n2g_vagg_plot + geom_jitter())
  Sys.sleep(5)
}

n2g_vagg_plot = epc_data %>% 
# filter(avul == 'compensational') %>% 
filter(between(IR, 0.3, 0.65)) %>% 
ggplot(data = ., aes(x = log(vertagg_rate), y = log(n2g)))
# ggplot(data = ., aes(x = vertagg_rate, y = n2g))
n2g_vagg_plot + geom_jitter()
# n2g_vagg_plot + geom_jitter(aes(color = avul))

epc_just_logvert = epc_data %>%
filter(between(IR, 0.3, 0.65)) %>% 
transmute(lvert = log(vertagg_rate), ln2g = log(n2g))

model_lvert_ln2g = lm(ln2g~lvert, data = epc_just_logvert)

library(segmented)

model_lvert_ln2g_seg = segmented(model_lvert_ln2g)

brpoint = model_lvert_ln2g_seg$psi[2]

n2g_vagg_plot_breakpoint = epc_data %>% 
# filter(avul == 'compensational') %>% 
filter(between(IR, 0.3, 0.65)) %>% 
ggplot(data = ., aes(x = log(vertagg_rate), y = log(n2g)))

n2g_vagg_plot_breakpoint + geom_jitter() + 
geom_vline(xintercept = brpoint, color = 'red3')
# n2g_vagg_plot + geom_jitter(aes(color = avul))

epc_filt = epc_data %>% 
# filter(vertagg_rate > exp(brpoint)) %>% 
filter(between(IR, 0.3, 0.55)) %>% 
filter(most_common_avulsion == 'compensational')

ggplot() + 
geom_point(aes(x = pres_percents, y = reworking), data = epc_data) + 
geom_point(aes(x = pres_percents, y = reworking), data = epc_filt, color = 'white', size = 0.5) + 
geom_abline()

n2g_nagg_plot = epc_filt %>% filter(avul == 'random') %>% 
ggplot(data = ., aes(x = nagg, y = n2g))

n2g_nagg_plot + geom_point(aes(color = IR))


n2g_vagg_plot = epc_filt %>% 
# filter(avul == 'compensational') %>% 
# filter(between(IR, 0.3, 0.65)) %>% 
ggplot(data = ., aes(x = log(vertagg_rate), y = log(n2g)))
# ggplot(data = ., aes(x = vertagg_rate, y = n2g))
n2g_vagg_plot + geom_jitter(aes(color = log(IR / vertagg_rate)))

epc_ratio = epc_filt %>% transmute(x = log10(IR/vertagg_rate), y = log10(n2g))

model_irvert_ln2g = lm(y~x, data = epc_ratio)

model_irvert_ln2g_seg = segmented(model_irvert_ln2g)

brpoint = model_irvert_ln2g_seg$psi[2]

n2g_vagg_plot = epc_filt %>% 
# filter(avul == 'compensational') %>% 
# filter(between(IR, 0.3, 0.65)) %>% 
ggplot(data = ., aes(x = log10(IR/vertagg_rate), y = log10(n2g)))
# ggplot(data = ., aes(x = vertagg_rate, y = n2g))
n2g_vagg_plot + geom_jitter() + geom_vline(xintercept = brpoint, colo = 'red')

n2g_pres_plot = epc_filt %>% 
# filter(most_common_avulsion == 'clustered') %>% 
# filter(IR == 0.4) %>% 
ggplot(data = ., aes(y = pres_percents, x = n2g)) 

n2g_pres_plot + geom_point(aes(color = latmob)) +
geom_line(aes(x = x, y = y, group = g), data = f)


epc_filt2 = epc_filt %>% mutate(presInv = 1 / (pres_percents-1.01), n2gInv = 1/n2g)
epc_filt2 %>% ggplot(aes(x = pres_percents, y = n2g)) + geom_point()

n2gMod = lm(n2gInv ~ presInv, data = epc_filt2, weight = 1-pres_percents)
reworkingMod = lm(reworking ~ pres_percents, data = epc_filt2)
reworkingModSeg = segmented(reworkingMod)


epc_filt3 = epc_filt2 %>% mutate(n2gPred = 1/predict(n2gMod)) %>% 
mutate(reworkingPred = predict(reworkingModSeg))

ggplot(aes(x = presInv), data = epc_filt3) + 
geom_point(aes(y = n2gInv)) + 
geom_line(aes(y = 1/n2gPred), color = 'red3')

ggplot(aes(x = pres_percents), data = epc_filt3) + 
geom_point(aes(y = n2g)) + 
geom_line(aes(y = n2gPred), color = 'red3')

ggplot(aes(x = pres_percents), data = epc_filt3) + 
geom_point(aes(y = reworking)) + geom_abline() + 
geom_vline(xintercept = reworkingModSeg$psi[2], color = 'red3') + 
geom_line(aes(y = reworkingPred), color = 'blue3')





nagg_IR_plot = epc_data %>% filter(avul == 'compensational') %>% 
filter(vertagg_rate == 0.305) %>% 
ggplot(data = ., aes(x = IR, y = nagg))

nagg_IR_plot + geom_jitter(aes(color = latmob))

nagg_latmob_plot = epc_data %>% # filter(most_common_avulsion == 'clustered') %>% 
# filter(latmob == 4) %>% 
ggplot(data = ., aes(x = latmob, y = n2g)) 

nagg_latmob_plot + geom_jitter() 

nagg_reworking_plot = epc_data %>% 
# filter(most_common_avulsion == 'clustered') %>% 
# filter(latmob == 4) %>% 
ggplot(data = ., aes(x = reworking, y = nagg)) 

# nagg_pres_plot + geom_point() 
nagg_reworking_plot + geom_point(aes(color = input_n2g)) 

nagg_input_plot = epc_data %>% filter(most_common_avulsion == 'random') %>% 
filter(IR == 0.4) %>% 
ggplot(data = ., aes(x = vertagg_rate, y = input_n2g)) 

nagg_input_plot + geom_jitter(aes(color = IR)) 

n2g_reworking_plot = epc_data %>% 
filter(most_common_avulsion == 'random') %>% 
# filter(latmob == 4) %>% 
ggplot(data = ., aes(x = reworking, y = n2g)) 

n2g_reworking_plot + geom_point() 
n2g_reworking_plot + geom_point(aes(color = nagg)) 

enrichment_plot = epc_data %>% filter(most_common_avulsion == 'clustered') %>% 
# filter(IR == 0.4) %>% 
ggplot(data = ., aes(x = input_n2g, y = n2g)) 

enrichment_plot + geom_jitter(aes(color = avul)) + geom_abline()
# enrichment_plot + geom_jitter() + geom_abline()


reworking_input_plot = epc_data %>%  filter(most_common_avulsion == 'random') %>% 
# filter(IR == 0.4) %>% 
ggplot(data = ., aes(x = input_n2g, y = reworking)) 

reworking_input_plot + geom_jitter(aes(color = n2g)) 
