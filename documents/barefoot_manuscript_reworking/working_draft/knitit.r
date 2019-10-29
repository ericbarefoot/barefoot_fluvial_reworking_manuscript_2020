#!/usr/bin/env Rscript

# generic script to take an input file and knit it to a output file.
# Stolen from the internet by Eric Barefoot
# Oct 2019

library(knitr)
# render_html()
# source("hooks.R") # mods to defaults
inFile = commandArgs(trailingOnly=TRUE)[1]
outFile = commandArgs(trailingOnly=TRUE)[2]
knit(inFile,output=outFile)
