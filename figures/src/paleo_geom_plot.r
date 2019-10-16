# function(s) to make a plot of data from 3D model mapping data and field data.
# specifically the data on kinematics vs geometric estimations
# Eric Barefoot
# Sep 2019

# first, load libraries

require(tidyverse)  # for ggplot and stuff
require(devtools)  # for the next thing
require(here)  # for directory sorting out.
require(colorspace)
library(gtable)
require(paleohydror)
if (interactive()) {
    require(colorout)
}

figdir = here("figures")

# load in the data

load(here("data", "derived_data", "piceance_paleoslope_data.rda")) # outcrop data

fieldDataGrouped = field_data %>% group_by(formation, set_id) %>% filter(set_id != 'FORMANA')

# make plot with field data of geometry-based paleohydraulic estimates

ggplot() +

# add panels with kinematic estimates.


# first step is to generate a dataset which has all the same attributes but
# estimates the depth based on paleohydro.

# group the data by formation and the set_id, which indicates a set of paleohydro
# measurements for a single location/outcrop in time. can only compare/combine
# measurements from the same set to calculate slope, etc.

fieldDataGrouped = field_data %>% group_by(formation, set_id) %>% filter(set_id != 'FORMANA')

fieldDataGrouped_BZF = field_data %>% group_by(formation, set_id)

# calculate the depth based on a x-set to depth calculation.  May fine-tune this
# later, but for now chooses the most ham-fisted option.  TODO: provide for
# selecting different methods based on field observations of certain kinds, or
# the availability of certain data.


flowDepthSetXset = fieldDataGrouped %>% filter(structure == "cross_set" & meas_type ==
    "thickness") %>% summarize(xsetDepth = suppressMessages(xset2depthSimple(meas_a)))

# calculate the depth based on channel scour fill thicknesses.  TODO: estimate
# some compaction estimation.

flowDepthSetChannel = fieldDataGrouped %>% filter(structure == "channel" & meas_type ==
    "thickness") %>% summarize(channelDepth = mean(meas_a))

# calculate the depth based on point bar thicknesses.

flowDepthSetBarThick = fieldDataGrouped %>% filter(structure == "bar" & meas_type ==
    "thickness") %>% summarize(barDepthThick = mean(meas_a))

# calculate the depth based on point bar y dimensions as measured from outcrop.

flowDepthSetBarDims = fieldDataGrouped %>% filter(structure == "bar" & meas_type ==
    "dimensions") %>% summarize(barDepthDims = mean(meas_b))

# merge these all back together into one tibble.

flowDepthsSetAll = Reduce(function(...) {
    full_join(..., by = c(formation = "formation", set_id = "set_id"))
}, list(flowDepthSetChannel, flowDepthSetBarThick, flowDepthSetBarDims))

# reshape for easier summary.

flowDepthsSetAll = flowDepthsSetAll %>% gather(key = "measType", value = "depth",
    c(-formation, -set_id))

# assign weights to each of the techniques for estimating depth. x-set
# measurements, though the most heavily weighted, have significant uncertainties,
# although simply by virtue of the fact that I have a lot of them and they are at
# all the sites, I am going to consider them the most important ones.

# They are, in order: xset thicknesses, channel thicknesses, bar thicknesses, and
# bar ydimensions (thicknesses) TODO: figure out how to weight bar measuments by
# grade. Future stuff.

weightByType = tibble(measType = unique(flowDepthsSetAll$measType), weight = c(10, 10, 10))

# join the weight vector to the table.

flowDepthsSetAll = left_join(flowDepthsSetAll, weightByType, by = "measType")

# group variables by formation, then by set id.

flowDepthSet = group_by(flowDepthsSetAll, formation, set_id)

# take a weighted mean of the different techniques by set, then by formation.

flowDepthSetMean = flowDepthSet %>% summarize(depth = weighted.mean(x = depth, w = weight,
    na.rm = T))
flowDepthFormMean = flowDepthSetMean %>% summarize(depth = mean(depth, na.rm = T))

# exploratory plots of the distributions.

# Now that we have depth, the next thing we really need is grain size. I don't
# have any data on this yet, so I'll just build the framework for now. The data
# for each sample is going to come in as a separate dataset. I can provide here a
# list of sample names.

sampleData = fieldDataGrouped %>% filter(meas_type == 'sample')

bedloadSamples = sampleData %>% filter(structure != 'outcrop' & structure != 'formation' & structure != 'bar_drape' & meas_a < 2000)

bedloadSamples = bedloadSamples %>% mutate(D50 = meas_a * 1e-6)

# sampCodes = as.vector(na.omit(unique(fieldDataGrouped$samp_code)))
#
# sample_info = fieldDataGrouped %>% filter(!is.na(samp_code)) %>% select(samp_code,
#     samp_description, textures, observations) %>% distinct()
#
# D50 = rep(NA, nrow(sample_info))
# D95 = rep(NA, nrow(sample_info))
#
# for (i in 1:nrow(sample_info)) {
#     str1 = sample_info$samp_description[i]
#     str2 = sample_info$textures[i]
#     if (grepl("bedload", str1)) {
#         if (grepl("coarse", str2)) {
#             D = rnorm(1, 5e-04, 1e-04)
#         } else if (grepl("fine", str2)) {
#             D = rnorm(1, 2e-04, 1e-04)
#         } else {
#             # message('this doesnt have a coarse-fine choice, going in between')
#             D = rnorm(1, 3e-04, 1e-04)
#         }
#         D50[i] = D
#     } else if (grepl("suspended", str1)) {
#         D = rnorm(1, 7.5e-05, 1e-05)
#         D95[i] = D
#     } else {
#
#     }
# }
#
# sizeGuess = tibble(D50 = D50, D95 = D95, samp_code = sampCodes)

# sampSizes = left_join(sample_info, sizeGuess, by = "samp_code")

# flowDepthSetGrainSize = left_join(flowDepthSetMean, ungroup(sampSizes) %>%
# select(-(samp_description:observations), -formation), by = "set_id")


flowDepthSetGrainSize = inner_join(flowDepthSetMean, bedloadSamples %>% ungroup() %>% group_by(set_id) %>% summarize(D50 = mean(D50)), by = "set_id")


# slopeShearEstimates = flowDepthSetGrainSize %>% mutate(S = trampush_slp(D = D50, H = depth)) %>% mutate(tau = 1000 * 9.8 * depth * S)
#
# # flowDepthSetGrainSize$D50[is.na(flowDepthSetGrainSize$samp_code)] = 3e-04
# # flowDepthSetGrainSize$samp_code[is.na(flowDepthSetGrainSize$samp_code)] = "20171007S6MMPB6014S006"
#
# S = numeric(length(unique(flowDepthSetGrainSize$set_id)))
# F = numeric(length(unique(flowDepthSetGrainSize$set_id)))
# j = 1
# for (i in sort(unique(flowDepthSetGrainSize$set_id))) {
#     dat = flowDepthSetGrainSize %>% filter(set_id == i)
#     Dbed = mean(dat$D50, na.rm = TRUE)
#     # Dsusp = mean(dat$D95, na.rm = TRUE)
#     H = mean(dat$depth, na.rm = TRUE)
#     S[j] = trampush_slp(D = Dbed, H = H)
#     if (is.na(Dsusp)) {
#     } else {
#         S[j] = trampush_slp(Dbed = Dbed, H = H)
#         F[j] = NA
#     }
#     j = j + 1
# }
#
# tibb_slopeF = tibble(set_id = sort(unique(flowDepthSetGrainSize$set_id)), S, F)
#
# paleoHydroSetEst = left_join(flowDepthSetGrainSize, tibb_slopeF, by = "set_id") %>%
#     mutate(tau = 1000 * 9.8 * depth * S)


paleoHydroSetEst = flowDepthSetGrainSize %>% mutate(S = trampush_slp(D = D50, H = depth)) %>% mutate(tau = 1000 * 9.8 * depth * S)

paleoHydroSetEst = ungroup(flowDepthSetGrainSize) %>% group_by(formation) %>% mutate(S = trampush_slp(D = D50, H = depth)) %>% mutate(tau = 1000 * 9.8 * depth * S)

paleoHydroFormEst = paleoHydroSetEst %>% summarize(h = median(depth), depthsd = sd(depth), slope = median(S), slopesd = sd(S), tau = median(tau), Dbed = median(D50, na.rm = TRUE), D50sd = sd(D50), n = n())

paleoHydroFormEst = paleoHydroSetEst %>% summarize(depth = median(depth), slope = median(S), tau = median(tau), D50 = median(D50, na.rm = TRUE), n = n())

cpal = hcl(h = c(80, 120, 240), c = rep(100, 3), l = c(85, 65, 20))

paleoHydroSetEst = paleoHydroSetEst %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 3, formation == 'atwell_gulch' ~ 2, formation == 'molina' ~ 1, formation == 'shire' ~ 0))

# df = field_data %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 3, formation == 'atwell_gulch' ~ 2, formation == 'molina' ~ 1, formation == 'shire' ~ 0))

legend_ord = levels(with(paleoHydroSetEst, reorder(formation, formID)))
#
# df %>% group_by(formID, loc_x, loc_y, formation) %>%
# summarize(n = n()) %>% filter(formation != 'ohio_creek') %>%
# ggplot() + aes(x = loc_x, y = loc_y, color = reorder(formation, formID)) +
# scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch','Ohio Creek'), breaks = legend_ord) +
# geom_point() + labs(x = 'Longitude', y = 'Latitude') + theme_minimal()

# slope_plot

slope_freq = paleoHydroSetEst %>% ggplot() + stat_bin(aes(x = S, y = ..density.., color = reorder(formation, formID)), geom = 'line', position = 'identity', bins = 8, size = 1.5) +
scale_x_log10(breaks = c(1e-04, 2e-04, 3e-04, 5e-04, 7e-04, 0.001, 0.002, 0.003), labels = c(expression(10^-4), expression(2%*%10^-4), expression(3%*%10^-4), expression(5%*%10^-4), expression(7%*%10^-4), expression(10^-3), expression(2%*%10^-3), expression(3%*%10^-3))) +
theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = 'Slope (-)', y = 'Probability Density')

# slope_freq

# depth plot

depth_freq = paleoHydroSetEst  %>% ggplot() + stat_bin(aes(x = depth, y = ..density.., color = reorder(formation, formID)), geom = 'line', position = "identity", binwidth = 0.7, size = 1.5) + theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = 'Flow Depth (m)', y = 'Probability Density')

# depth_freq

# grain size plot

grainsize_freq = paleoHydroSetEst %>% ggplot() + stat_bin(aes(x = D50 * 1e6, y = ..density.., color = reorder(formation, formID)), geom = 'line', position = "identity", binwidth = 65, size = 1.5) + theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = expression(paste('Bedload Sample ', D[50], ' (', mu, 'm)')), y = 'Probability Density')

# grainsize_freq

# slope_plot

slope_curves = paleoHydroSetEst %>% ggplot() + stat_density(aes(x = S, color = reorder(formation, formID)), geom = 'line', position = "identity", bw = 0.15) +
scale_x_log10(breaks = c(1e-04, 2e-04, 3e-04, 5e-04, 7e-04, 0.001, 0.002, 0.003), labels = c(expression(10^-4), expression(2%*%10^-4), expression(3%*%10^-4), expression(5%*%10^-4), expression(7%*%10^-4), '0.001', '0.002', '0.003')) + theme_minimal() +
scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = 'Slope (-)', y = 'Probability Density')

# slope_curves
# depth plot

depth_curves = paleoHydroSetEst  %>% ggplot() + stat_density(aes(x = depth, color = reorder(formation, formID)), geom = 'line', position = "identity", bw = 0.7) + theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = 'Flow Depth (m)', y = 'Probability Density')

# grain size plot

grainsize_curves = paleoHydroSetEst %>% ggplot() + stat_density(aes(x = D50 * 1e6, color = reorder(formation, formID)), geom = 'line', position = "identity", bw = 50) + theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = expression(paste('Sample ', D[50], ' (', mu, 'm)')), y = 'Probability Density')

# ggsave(plot = slope_curves, filename = "petm_slope_piceance_density_estimates.pdf", path = figdir, width = 7, height = 4, units = "in")
#
# ggsave(plot = grainsize_curves, filename = "petm_grainsize_piceance_density_estimates.pdf", path = figdir, width = 7, height = 4, units = "in")
#
# ggsave(plot = depth_curves, filename = "petm_depth_piceance_density_estimates.pdf", path = figdir, width = 7, height = 4, units = "in")

ggsave(plot = slope_freq, filename = "petm_slope_piceance_frequency.png",
    path = figdir, width = 7, height = 4, units = "in")

ggsave(plot = grainsize_freq, filename = "petm_grainsize_piceance_frequency.png",
    path = figdir, width = 7, height = 4, units = "in")

ggsave(plot = depth_freq, filename = "petm_depth_piceance_frequency.png",
    path = figdir, width = 7, height = 4, units = "in")

save(file = here("data", "derived_data", "paleoHydroFormationEstimates.rda"), paleoHydroFormEst)
save(file = here("data", "derived_data", "paleoHydroSetEstimates.rda"), paleoHydroSetEst)


# formNames = list('Ohio Creek','Shire','Molina','Atwell Gulch')
#
# formlabels = function(var,v) {
#   return(formNames[[v]])
# }

SizeplotDepth = paleoHydroSetEst %>% ggplot() + theme_minimal() + scale_color_manual(values = cpal, name = 'Member', labels = c('Shire','Molina','Atwell Gulch'), breaks = legend_ord) + labs(x = expression(paste('Grain Size (', mu, 'm)')), y = 'Flow Depth (m)') + facet_grid(rows = vars(formation)) +
geom_smooth(method = 'lm', aes(x = D50 * 1e6, y = depth, color = reorder(formation, formID)), fill = hcl(c = 0, l = 85), level = 0.99) + theme(legend.position = 'bottom') + geom_point(aes(x = D50 * 1e6, y = depth, color = reorder(formation, formID)), position = "identity")

ggsave(plot = SizeplotDepth, filename = "petm_depth_grainsize_piceance.pdf",
    path = figdir, width = 4, height = 7, units = "in")

# and once we have a list of the names, which will be the data frame names, we
# just need to calculate the D50 of each, and the Dmaxs of each, and compile them
# here.

######################################################################################################################## Some code extracting the D50 from each sample code, ending in a predictably
######################################################################################################################## named dataset grainSizeSamples

# left_join(fieldDataGrouped, grainSizeSamples, by = 'set_id')

######################################################################################################################## Some code generating some fake grainsize data based on notes.  coarse/medium
######################################################################################################################## sand = 500 micron medium/fine sand = 200 micron medium silt = 20 micron


paleoHydroFormEst = paleoHydroFormEst %>% mutate(formID = case_when(formation == 'ohio_creek' ~ 0, formation == 'atwell_gulch' ~ 1, formation == 'molina' ~ 2, formation == 'shire' ~ 3))


stratdat = paleoHydroFormEst %>% mutate(D50 = D50 * 10^6, slope = slope * 10^4) %>% gather(var, data, -c(formation,formID))

stratplot = stratdat %>% filter(var != 'n') %>% ggplot() + aes(y = formID, x = data) + geom_path() + geom_point() + facet_grid(cols = vars(var), scales = 'free_x') + theme_minimal()

stratplot
